"""Verify and analyze the preregistered MP2 comparison without running new energies."""
import argparse, os, sys, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
ap=argparse.ArgumentParser()
ap.add_argument("--packages",required=True,type=Path)
args=ap.parse_args()
sys.path.insert(0,str(args.packages.resolve()))
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
os.environ["MPLCONFIGDIR"]=str(ROOT/"plot-cache")
import numpy as np
from pyscf import gto
from pyscf.hessian import thermo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
read=lambda f:json.loads(f.read_text(encoding="utf-8-sig"))
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
def write(name,data):
    (ROOT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")
p=read(ROOT/"protocol.json")
receipt=read(ROOT/"execution-receipt.json")
start=read(ROOT/"run-start.json")
prereg=read(ROOT/"preregistration-receipt.json")
ph=sha(ROOT/"protocol.json")
assert receipt["status"]=="completed"
assert ph==receipt["protocol_sha256"]==start["protocol_sha256"]==prereg["protocol_sha256"]
assert sha(ROOT/"run_comparison.py")==receipt["runner_sha256"]==start["runner_sha256"]
assert sha(ROOT/"implementation-receipt.json")==start["implementation_receipt_sha256"]
assert sha(ROOT/"preregistration-receipt.json")==start["preregistration_receipt_sha256"]
assert datetime.fromisoformat(prereg["checksum_corrected_before_calculations_utc"]) < datetime.fromisoformat(start["started_at_utc"])
for name,digest in p["source_inputs_sha256"].items():
    assert sha(ROOT/"source"/name)==digest
records={}
for row in receipt["evaluations"]:
    f=ROOT/row["file"]
    assert sha(f)==row["sha256"]
    rec=read(f)
    assert rec["protocol_sha256"]==ph and rec["id"]==row["id"]
    records[row["id"]]=rec
assert len(records)==receipt["energy_gradient_evaluation_count"]
opt=read(ROOT/"optimizations.json")
raw=read(ROOT/"hessians.json")
xyz=np.array(opt[0]["geometry_bohr"])
mol=gto.M(atom=[(s,tuple(r)) for s,r in zip(["O","H","H"],xyz)],unit="Bohr",basis="cc-pvtz",verbose=0)
checks=[]
def check(name,value,limit,scope="global"):
    checks.append({"name":name,"value":float(value),"limit":float(limit),"passes":bool(value<=limit),"scope":scope})
v=p["numerical_validity"]
check("all recorded RHF solves converged",sum(not r["RHF_converged"] for r in records.values()),0)
check("stationarity at both optimized geometries",max(o["max_cartesian_gradient_Eh_per_bohr"] for o in opt),v["maximum_cartesian_gradient_Eh_per_bohr"])
check("starting-geometry agreement: bond length",abs(opt[0]["q"][0]-opt[1]["q"][0]),v["maximum_start_geometry_difference_bohr"])
check("starting-geometry agreement: angle",abs(opt[0]["q"][1]-opt[1]["q"][1]),v["maximum_start_angle_difference_rad"])
check("starting-geometry agreement: energy",abs(opt[0]["energy_Eh"]-opt[1]["energy_Eh"]),v["maximum_start_energy_difference_Eh"])
C=5140.487143715828
rows=[]
def modes(H,masses):
    sm=np.sqrt(np.repeat(masses,3))
    mw=H/sm[:,None]/sm[None,:]
    com=np.sum(xyz*masses[:,None],axis=0)/sum(masses)
    tr=[]
    for axis in np.eye(3):
        tr.append((np.tile(axis,(3,1))*np.sqrt(masses)[:,None]).reshape(-1))
        tr.append((np.cross(np.tile(axis,(3,1)),xyz-com)*np.sqrt(masses)[:,None]).reshape(-1))
    U,S,_=np.linalg.svd(np.array(tr).T,full_matrices=True)
    assert sum(S>1e-10)==6
    Z=U[:,6:]
    vals,vec=np.linalg.eigh(Z.T@mw@Z)
    frequencies=np.sign(vals)*np.sqrt(abs(vals))*C
    cart=(Z@vec/sm[:,None]).T.reshape(3,3,3)
    bond_units=np.array([(xyz[i]-xyz[0])/np.linalg.norm(xyz[i]-xyz[0]) for i in [1,2]])
    stretches=np.array([[np.dot(bond_units[j],mode[j+1]-mode[0]) for j in range(2)] for mode in cart])
    anti_score=(stretches[:,0]-stretches[:,1])**2/(2*np.sum(stretches**2,axis=1))
    a=int(np.argmax(anti_score))
    s=max([i for i in range(3) if i!=a],key=lambda i:frequencies[i])
    b=next(i for i in range(3) if i not in [a,s])
    identified=bool(anti_score[a]>.99 and stretches[a,0]*stretches[a,1]<0 and stretches[s,0]*stretches[s,1]>0)
    native=thermo.harmonic_analysis(mol,H.reshape(3,3,3,3).transpose(0,2,1,3),mass=masses)
    nf=np.asarray(native["freq_wavenumber"])
    return {"frequencies_sorted_cm1":frequencies.tolist(),
            "bend_cm1":float(frequencies[b]),"symmetric_stretch_cm1":float(frequencies[s]),
            "antisymmetric_stretch_cm1":float(frequencies[a]),
            "gap_cm1":float(frequencies[a]-frequencies[s]),
            "internal_eigenvalues":vals.tolist(),"mode_bond_changes":stretches.tolist(),
            "antisymmetric_scores":anti_score.tolist(),"modes_identified":identified,
            "independent_pyscf_frequencies_cm1":np.real(nf).tolist(),
            "independent_max_frequency_difference_cm1":float(np.max(abs(np.sort(np.real(nf))-frequencies))),
            "native_imaginary_modes":int(native["freq_error"])}
for item in raw:
    h=item["step_bohr"]
    reconstructed=np.zeros((9,9))
    for pair in item["evaluation_references"]:
        axis=pair["axis"]
        pp,mm=records[pair["plus"]],records[pair["minus"]]
        d=np.zeros((3,3));d.reshape(-1)[axis]=h
        assert np.max(abs(np.array(pp["geometry_bohr"])-(xyz+d)))<1e-14
        assert np.max(abs(np.array(mm["geometry_bohr"])-(xyz-d)))<1e-14
        reconstructed[:,axis]=(np.array(pp["gradient_Eh_per_bohr"]).reshape(-1)-np.array(mm["gradient_Eh_per_bohr"]).reshape(-1))/(2*h)
    assert np.array_equal(reconstructed,np.array(item["raw_hessian_Eh_per_bohr2"]))
    H=(reconstructed+reconstructed.T)/2
    assert np.array_equal(H,np.array(item["symmetrized_hessian_Eh_per_bohr2"]))
    check(f"step {h}: Hessian antisymmetry",np.max(abs(reconstructed-reconstructed.T)),v["maximum_hessian_antisymmetry_Eh_per_bohr2"])
    row={"step_bohr":h,"isotopes":{}}
    for iso,key in [("H2O","H1"),("D2O","D2")]:
        masses=np.array([p["mass_u"]["O16"]]+[p["mass_u"][key]]*2)
        out=modes(H,masses)
        row["isotopes"][iso]=out
        check(f"step {h}: nonpositive internal modes",sum(x<=0 for x in out["internal_eigenvalues"]),0,iso)
        check(f"step {h}: independent PySCF imaginary modes",out["native_imaginary_modes"],0,iso)
        check(f"step {h}: mode assignment unresolved",not out["modes_identified"],0,iso)
        check(f"step {h}: independent mode calculation agreement",out["independent_max_frequency_difference_cm1"],v["maximum_independent_mode_frequency_difference_cm1"],iso)
    rows.append(row)
for iso in ["H2O","D2O"]:
    aa,bb=[r["isotopes"][iso] for r in rows[-2:]]
    check("last two steps: maximum frequency change",np.max(abs(np.array(aa["frequencies_sorted_cm1"])-np.array(bb["frequencies_sorted_cm1"]))),v["maximum_last_two_steps_frequency_change_cm1"],iso)
    check("last two steps: gap change",abs(aa["gap_cm1"]-bb["gap_cm1"]),v["maximum_last_two_steps_gap_change_cm1"],iso)
primary=next(r for r in rows if r["step_bohr"]==p["hessian"]["primary_reported_step_bohr"])["isotopes"]
valid=all(c["passes"] for c in checks if c["scope"] in ["global","H2O"])
test=p["primary_gap_test"]
error=abs(primary["H2O"]["gap_cm1"]-test["reference_gap_cm1"])
status=("PASS" if error<=test["absolute_gap_error_limit_cm1"] else "FAIL") if valid else "INCONCLUSIVE"
old=read(ROOT/"source/stretch-gap-results.json")
baseline={r["isotope"]:r for r in old["results"]}
result={"label":"Separate preregistered MP2/cc-pVTZ comparison","completed_at_utc":datetime.now(timezone.utc).isoformat(),
        "protocol_sha256":ph,"analyzer_sha256":sha(Path(__file__)),"run_started_at_utc":start["started_at_utc"],
        "registered_at_utc":p["registered_before_any_MP2_energy_or_optimization_utc"],
        "versions":receipt["versions"],"energy_gradient_evaluations":len(records),
        "primary_method":"all-electron RMP2/cc-pVTZ, no density fitting, unscaled harmonic frequencies",
        "primary_reported_step_bohr":p["hessian"]["primary_reported_step_bohr"],
        "primary_H2O_outcome":status,"H2O_absolute_gap_error_cm1":error,
        "H2O_gap_reference_cm1":test["reference_gap_cm1"],"H2O_absolute_gap_error_limit_cm1":test["absolute_gap_error_limit_cm1"],
        "H2O_numerical_validity_passes":valid,
        "numerical_checks":checks,"all_numerical_checks_pass":all(c["passes"] for c in checks),
        "primary_results":primary,"all_steps":rows,
        "optimized_geometries":[{k:r[k] for k in ["start","r_angstrom","angle_deg","optimizer_success","optimizer_message","max_cartesian_gradient_Eh_per_bohr"]} for r in opt],
        "historical_xTB_diagnostic":baseline,"historical_xTB_harmonic_diagnostic":old["h2o_harmonic_comparison"],
        "historical_results_regraded":False,"D2O_status":"diagnostic only; harmonic prediction compared with a fundamental reference is not a like-for-like test",
        "physical_parameters_tuned":False,"lean_run":False,"published_to_github":False}
write("comparison-results.json",result)
h2,d2=primary["H2O"],primary["D2O"]
fig,axs=plt.subplots(1,2,figsize=(13,5),gridspec_kw={"width_ratios":[1.25,1]},layout="constrained")
a,b=axs
sets=[("Published xTB",3643.3059614053795,3651.4359004635226,"#607D8B"),
      ("New MP2",h2["symmetric_stretch_cm1"],h2["antisymmetric_stretch_cm1"],"#166C85"),
      ("Harmonic reference",3832.,3943.,"#242424")]
for y,(label,sym,anti,color) in enumerate(sets):
    a.plot([sym,anti],[y,y],color=color,lw=3,marker="o")
    a.annotate(f"{anti-sym:.2f}",((sym+anti)/2,y),xytext=(0,10),textcoords="offset points",ha="center",color=color)
a.set_yticks(range(3),[s[0] for s in sets]);a.invert_yaxis()
a.set_ylim(2.5,-.6);a.set_xlabel("Stretch-band position (cm$^{-1}$)")
a.set_title("Band positions and their gap",loc="left",weight="bold")
a.grid(axis="x",alpha=.2)
b.axvspan(101,121,color="#E1EFE8",label="Predeclared pass interval")
b.axvline(111,color="#242424",ls="--",lw=1,label="Harmonic reference")
gap_old=float(baseline["H2O"]["predicted_gap_cm1"])
for y,(label,gap,color) in enumerate([("Published xTB",gap_old,"#607D8B"),("New MP2",h2["gap_cm1"],"#166C85")]):
    b.scatter([gap],[y],s=75,color=color,zorder=3)
    b.annotate(f"{gap:.2f}",(gap,y),xytext=(0,12),textcoords="offset points",ha="center")
b.set_yticks([0,1],["Published xTB\nDiagnostic","New MP2\n"+status])
b.set_ylim(1.8,-.8);b.set_xlim(0,max(130,h2["gap_cm1"]+15))
b.set_xlabel("H2O stretch gap (cm$^{-1}$)");b.set_title("Rule fixed before the run: 101-121",loc="left",weight="bold")
b.grid(axis="x",alpha=.2);b.legend(loc="lower left",fontsize=8)
fig.suptitle("A separate MP2 comparison of water's stretch gap",fontsize=16,weight="bold")
fig.savefig("stretch-gap-comparison.png",dpi=170)
fig.savefig("stretch-gap-comparison.svg")
plt.close(fig)
passed=sum(c["passes"] for c in checks)
settling={iso:next(c["value"] for c in checks if c["scope"]==iso and c["name"]=="last two steps: gap change") for iso in ["H2O","D2O"]}
report=f"""# A separate MP2 comparison of the stretch gap

**{status}:** MP2/cc-pVTZ gives an H2O harmonic stretch gap of **{h2["gap_cm1"]:.3f} cm^-1**. Before any molecular calculation, we chose a pass rule of **101-121 cm^-1**: within 10 cm^-1 of the same 111 cm^-1 harmonic reference used in the earlier xTB diagnostic. The error here is {error:.3f} cm^-1. This is a chosen project rule, not an established accuracy standard.

![Comparison of stretch bands and the predeclared gap interval](stretch-gap-comparison.png)

| H2O harmonic quantity (cm^-1) | Published GFN2-xTB | New MP2/cc-pVTZ | Harmonic reference |
|---|---:|---:|---:|
| Symmetric stretch | 3643.306 | {h2["symmetric_stretch_cm1"]:.3f} | 3832 |
| Antisymmetric stretch | 3651.436 | {h2["antisymmetric_stretch_cm1"]:.3f} | 3943 |
| Gap | {gap_old:.3f} | {h2["gap_cm1"]:.3f} | 111 |

Each method uses its own optimized gas-phase geometry and unscaled harmonic frequencies. The reference comes from the [NIST CCCBDB H2O table](https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0). The measured fundamental gap of about 99 belongs to a different comparison; it was not used to grade this harmonic prediction.

MP2 is a different electronic-structure method that adds a correlation correction to a Hartree-Fock calculation. We used the fixed cc-pVTZ basis, all electrons and no density fitting. No coefficient or frequency scale was fitted. [PySCF method documentation](https://pyscf.org/user/mp.html).

The [protocol](protocol.json) was recorded at {p["registered_before_any_MP2_energy_or_optimization_utc"]}. The first calculation began at {start["started_at_utc"]}. An initial preflight caught a Windows line-ending mismatch in the receipt's checksum before any molecular calculation. The original receipt and log are retained; the corrected receipt hashes the unchanged saved protocol. No method, reference or pass criterion changed.

We completed {len(records)} energy-and-gradient evaluations, optimizing from two starting geometries and forming Hessians at displacements of 0.005, 0.0025 and 0.00125 bohr. The finest step was designated for reporting before running. **{passed} of {len(checks)} numerical controls passed.** Halving the final displacement changed the H2O gap by {settling["H2O"]:.6f} cm^-1 and the D2O gap by {settling["D2O"]:.6f} cm^-1. The controls also cover geometry agreement, stationarity, Hessian symmetry, positive internal frequencies, mode identification and agreement with PySCF's separate harmonic-analysis routine. Full values and optimizer statuses are in [the results](comparison-results.json).

For D2O, the new harmonic prediction gives stretches of {d2["symmetric_stretch_cm1"]:.3f} and {d2["antisymmetric_stretch_cm1"]:.3f} cm^-1, a gap of **{d2["gap_cm1"]:.3f} cm^-1**, compared with the saved xTB gap of 54.668. This remains diagnostic: the checked [NIST D2O table](https://cccbdb.nist.gov/exp2x.asp?casno=7789200&charge=0) lists fundamental values but no harmonic values. We did not grade a harmonic prediction against its approximately 117 cm^-1 fundamental gap.

This comparison tests the specified MP2 method and finite basis for this molecule. It does not isolate which xTB term is physically inaccurate, establish basis-set convergence, or validate a liquid-water simulation. Those are separate questions.

The original xTB position PASS under the 5% rule and its gap diagnostics (8 against 99 for H2O, 55 against 117 for D2O, and H2O harmonic 8 against 111) are unchanged. The new gap rule was not applied retrospectively to rename the published xTB outcome. No Lean check was run. These results are a local addition and have not been published to GitHub.

Reproduction: use the [locked packages](requirements-lock.txt), [runner](run_comparison.py) and [analyzer](analyze_comparison.py). Both accept --packages pointing to the isolated package directory. The runner refuses to overwrite an existing run. The [execution receipt](execution-receipt.json) hashes every saved energy/gradient record; the analyzer reconstructs each Hessian from those records and verifies the protocol, runner and implementation hashes before computing the result.

Generated by Codex under the user's direction.
"""
(ROOT/"README.md").write_text(report,encoding="utf-8")
print(json.dumps({"H2O_outcome":status,"H2O_gap_cm1":h2["gap_cm1"],"H2O_error_cm1":error,
                  "D2O_gap_cm1":d2["gap_cm1"],"numerical_controls_passed":passed,"numerical_controls_total":len(checks),
                  "energy_gradient_evaluations":len(records),"failed_controls":[c for c in checks if not c["passes"]]}))
