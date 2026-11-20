"""Check the MP2 gradient response equation numerically; keep the primary run intact."""
import argparse,os,sys,json,hashlib,time,contextlib,io
from pathlib import Path
from datetime import datetime,timezone
ap=argparse.ArgumentParser();ap.add_argument("--packages",required=True,type=Path);args=ap.parse_args()
sys.path.insert(0,str(args.packages.resolve()))
ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent
os.chdir(ROOT)
for key in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"]:os.environ[key]="2"
(ROOT/"scratch").mkdir(exist_ok=True)
os.environ["PYSCF_TMPDIR"]=str(ROOT/"scratch")
import numpy as np
from scipy.sparse.linalg import gmres,LinearOperator
from pyscf import gto,scf,mp,lib
from pyscf.scf import cphf
from pyscf.hessian import thermo
lib.num_threads(2)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding="utf-8-sig"))
now=lambda:datetime.now(timezone.utc).isoformat()
def write(name,obj):(ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n",encoding="utf-8")
p=read(ROOT/"protocol.json");base=read(PARENT/"protocol.json");prior=read(PARENT/"comparison-results.json")
ph=sha(ROOT/"protocol.json")
assert ph==read(ROOT/"preregistration-receipt.json")["protocol_sha256"]
for name,field in [("protocol.json","parent_protocol_sha256"),("comparison-results.json","parent_results_sha256"),("hessians.json","parent_hessians_sha256"),("optimizations.json","parent_optimizations_sha256")]:
    assert sha(PARENT/name)==p[field]
assert not (ROOT/"run-start.json").exists()
started=now()
write("run-start.json",{"started_at_utc":started,"protocol_sha256":ph,"runner_sha256":sha(Path(__file__))})
assert datetime.fromisoformat(started)>datetime.fromisoformat(p["written_before_followup_calculations_utc"])
response_records=[]
def checked_response(fvind,mo_energy,mo_occ,h1,s1=None,max_cycle=50,tol=1e-9,hermi=False,verbose=None,level_shift=0):
    assert s1 is None and level_shift==0 and not hermi and h1.ndim==2
    gap=mo_energy[mo_occ==0,None]-mo_energy[None,mo_occ>0]
    assert h1.shape==gap.shape and np.min(gap)>0
    shape=h1.shape;n=h1.size
    def op(x):
        x=np.asarray(x).reshape(shape)
        return (x+np.asarray(fvind(x)).reshape(shape)/gap).reshape(-1)
    A=LinearOperator((n,n),matvec=op,dtype=np.float64)
    b=(-h1/gap).reshape(-1)
    calls=[]
    x,info=gmres(A,b,rtol=p["solver"]["rtol"],atol=p["solver"]["atol"],
                 restart=min(p["solver"]["restart"],n),maxiter=p["solver"]["maxiter"],
                 callback=lambda value:calls.append(float(value)),callback_type="pr_norm")
    residual=float(np.linalg.norm(op(x)-b))
    response_records.append({"dimension":n,"gmres_info":int(info),"iterations":len(calls),
                            "true_preconditioned_residual_norm":residual})
    assert info==0 and residual<=p["validation_limits"]["maximum_true_preconditioned_response_residual"]
    return x.reshape(shape),None
cphf.solve=checked_response
began=time.monotonic()
xyz=np.array(read(PARENT/"optimizations.json")[0]["geometry_bohr"])
h=p["step_bohr"]
cases=[("baseline",xyz)]
for j in range(9):
    d=np.zeros(9);d[j]=h
    cases.extend([(f"axis-{j}-plus",xyz+d.reshape(3,3)),(f"axis-{j}-minus",xyz-d.reshape(3,3))])
(ROOT/"evaluations").mkdir()
records={}
receipts=[]
for tag,geom in cases:
    assert time.monotonic()-began<p["bounded_work"]["maximum_wall_seconds"]
    output=io.StringIO()
    previous=len(response_records)
    with contextlib.redirect_stdout(output),contextlib.redirect_stderr(output):
        mol=gto.M(atom=[(s,tuple(r)) for s,r in zip(["O","H","H"],geom)],unit="Bohr",
                  basis="cc-pvtz",charge=0,spin=0,symmetry=False,verbose=0,max_memory=base["resources"]["memory_MB"])
        mf=scf.RHF(mol);mf.chkfile=None
        mf.conv_tol=base["scf"]["conv_tol_Eh"];mf.conv_tol_grad=base["scf"]["conv_tol_grad"]
        mf.direct_scf_tol=base["scf"]["direct_scf_tol"];mf.max_cycle=base["scf"]["max_cycle"]
        mf.kernel();assert mf.converged
        model=mp.MP2(mf,frozen=0);model.kernel()
        grad=model.nuc_grad_method().kernel()
    assert len(response_records)==previous+1
    rec={"tag":tag,"completed_at_utc":now(),"protocol_sha256":ph,"geometry_bohr":geom.tolist(),
         "energy_Eh":float(model.e_tot),"gradient_Eh_per_bohr":grad.tolist(),
         "RHF_converged":bool(mf.converged),"response":response_records[-1],"captured_messages":output.getvalue()}
    records[tag]=rec
    f=ROOT/"evaluations"/(tag+".json")
    f.write_text(json.dumps(rec,indent=2)+"\n",encoding="utf-8")
    receipts.append({"file":f.relative_to(ROOT).as_posix(),"sha256":sha(f)})
    print(json.dumps({"completed":tag,"response_residual":response_records[-1]["true_preconditioned_residual_norm"]}),flush=True)
assert len(records)==p["bounded_work"]["maximum_gradient_evaluations"]
H=np.zeros((9,9))
for j in range(9):
    H[:,j]=(np.array(records[f"axis-{j}-plus"]["gradient_Eh_per_bohr"]).reshape(-1)-np.array(records[f"axis-{j}-minus"]["gradient_Eh_per_bohr"]).reshape(-1))/(2*h)
rawH=H.copy();H=(H+H.T)/2
write("hessian.json",{"step_bohr":h,"raw_hessian_Eh_per_bohr2":rawH.tolist(),"symmetrized_hessian_Eh_per_bohr2":H.tolist()})
mol=gto.M(atom=[(s,tuple(r)) for s,r in zip(["O","H","H"],xyz)],unit="Bohr",basis="cc-pvtz",verbose=0)
checks=[]
def check(name,value,limit):checks.append({"name":name,"value":float(value),"limit":limit,"passes":bool(value<=limit)})
lim=p["validation_limits"]
check("maximum true response residual",max(x["true_preconditioned_residual_norm"] for x in response_records),lim["maximum_true_preconditioned_response_residual"])
check("baseline stationarity",np.max(abs(np.array(records["baseline"]["gradient_Eh_per_bohr"]))),lim["maximum_primary_geometry_gradient_Eh_per_bohr"])
out={}
for iso,key in [("H2O","H1"),("D2O","D2")]:
    mass=np.array([base["mass_u"]["O16"]]+[base["mass_u"][key]]*2)
    x=thermo.harmonic_analysis(mol,H.reshape(3,3,3,3).transpose(0,2,1,3),mass=mass)
    assert x["freq_error"]==0
    f=np.sort(np.real(x["freq_wavenumber"]))
    old=prior["primary_results"][iso]
    assert old["frequencies_sorted_cm1"][1]==old["symmetric_stretch_cm1"]
    assert old["frequencies_sorted_cm1"][2]==old["antisymmetric_stretch_cm1"]
    gap=float(f[2]-f[1]);delta=abs(gap-old["gap_cm1"])
    check(iso+" frequency change",np.max(abs(f-np.array(old["frequencies_sorted_cm1"]))),lim["maximum_absolute_frequency_change_cm1"])
    check(iso+" gap change",delta,lim["maximum_absolute_gap_change_cm1"])
    out[iso]={"frequencies_cm1":f.tolist(),"gap_cm1":gap,"primary_gap_cm1":old["gap_cm1"],"gap_change_cm1":delta}
result={"label":"Separate response-solver verification","protocol_sha256":ph,"runner_sha256":sha(Path(__file__)),
        "started_at_utc":started,"finished_at_utc":now(),"elapsed_seconds":time.monotonic()-began,
        "evaluation_count":len(records),"evaluations":receipts,"numerical_checks":checks,
        "all_validation_checks_pass":all(c["passes"] for c in checks),"results":out,
        "original_primary_result_preserved":True,"physical_method_and_gap_rule_unchanged":True,
        "H2O_original_outcome":prior["primary_H2O_outcome"],
        "H2O_followup_gap_inside_original_interval":101<=out["H2O"]["gap_cm1"]<=121}
write("verification-results.json",result)
assert sha(PARENT/"comparison-results.json")==p["parent_results_sha256"]
print(json.dumps({"validation_passes":result["all_validation_checks_pass"],"results":out,"controls":checks}),flush=True)
