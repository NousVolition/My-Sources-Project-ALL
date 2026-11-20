"""Run the preregistered MP2 comparison; never edit the frozen protocol."""
import argparse, os, sys, json, hashlib, time, inspect
from pathlib import Path
from datetime import datetime, timezone
ap = argparse.ArgumentParser()
ap.add_argument("--packages", required=True, type=Path)
args = ap.parse_args()
sys.path.insert(0, str(args.packages.resolve()))
ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
for key in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[key] = "2"
(ROOT / "scratch").mkdir(exist_ok=True)
os.environ["PYSCF_TMPDIR"] = str(ROOT / "scratch")
import numpy as np
import scipy
from scipy.optimize import minimize
import pyscf
from pyscf import gto, scf, mp, lib
from pyscf.grad import mp2 as mp2_grad
from pyscf.scf import cphf
from pyscf.hessian import thermo
from importlib import metadata

lib.num_threads(2)
BOHR_A = 0.529177210903
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
read = lambda p: json.loads(p.read_text(encoding="utf-8"))
now = lambda: datetime.now(timezone.utc).isoformat()
def write(name, data):
    (ROOT / name).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
p = read(ROOT / "protocol.json")
prereg = read(ROOT / "preregistration-receipt.json")
ph = sha(ROOT / "protocol.json")
assert ph == prereg["protocol_sha256"]
assert not (ROOT / "run-start.json").exists(), "A run already started; inspect its records rather than rerun."
for name, digest in p["source_inputs_sha256"].items():
    assert sha(ROOT / "source" / name) == digest
versions = {name:metadata.version(name) for name in ["pyscf","numpy","scipy","h5py","matplotlib"]}
source_hashes = {}
for mod in [mp.mp2, mp2_grad, scf.hf, cphf, thermo, gto.mole]:
    f = Path(inspect.getsourcefile(mod))
    source_hashes[f.relative_to(Path(pyscf.__file__).parent).as_posix()] = sha(f)
libraries = {}
libroot = Path(pyscf.__file__).parent / "lib"
for f in libroot.rglob("*.so*"):
    if f.is_file():
        libraries[f.relative_to(libroot).as_posix()] = sha(f)
implementation = {
    "recorded_before_first_calculation_utc":now(), "protocol_sha256":ph,
    "runner_sha256":sha(Path(__file__)), "versions":versions,
    "python_version":sys.version, "pyscf_source_hashes":source_hashes,
    "pyscf_library_hashes":libraries,
    "basis":{s:gto.basis.load("cc-pvtz",s) for s in ["H","O"]},
    "gradient_response":{"solver":"PySCF default CPHF","tolerance":inspect.signature(cphf.solve).parameters["tol"].default,
                         "maximum_cycles_in_MP2_gradient":30},
    "BOHR_A":BOHR_A, "new_calculations_executed":False,
}
write("implementation-receipt.json",implementation)
dists = sorted((d.metadata["Name"],d.version) for d in metadata.distributions(path=[str(args.packages.resolve())]))
(ROOT / "requirements-lock.txt").write_text("\n".join(f"{n}=={v}" for n,v in dists)+"\n",encoding="utf-8")
start = now()
assert datetime.fromisoformat(start) > datetime.fromisoformat(p["registered_before_any_MP2_energy_or_optimization_utc"])
write("run-start.json",{"started_at_utc":start,"protocol_sha256":ph,
    "preregistration_receipt_sha256":sha(ROOT/"preregistration-receipt.json"),
    "implementation_receipt_sha256":sha(ROOT/"implementation-receipt.json"),
    "runner_sha256":sha(Path(__file__)),"new_method":"MP2/cc-pVTZ"})
(ROOT / "evaluations").mkdir()
began = time.monotonic()
records, cache = [], {}
def coords(q):
    r,t = q
    return np.array([[0.,0.,0.],[r*np.sin(t/2),0.,r*np.cos(t/2)],[-r*np.sin(t/2),0.,r*np.cos(t/2)]])
def evaluate(xyz, stage):
    key = np.ascontiguousarray(xyz).tobytes()
    if key in cache:
        return cache[key]
    assert len(records) < p["resources"]["maximum_energy_gradient_evaluations"]
    assert time.monotonic()-began < p["resources"]["maximum_calculation_wall_seconds"]
    assert sha(ROOT/"protocol.json")==ph
    tag = f"eval-{len(records)+1:04d}"
    t = time.monotonic()
    mol = gto.M(atom=[(s,tuple(x)) for s,x in zip(["O","H","H"],xyz)],
                unit="Bohr", basis="cc-pvtz", charge=0, spin=0, symmetry=False,
                verbose=0, max_memory=p["resources"]["memory_MB"])
    mf = scf.RHF(mol)
    mf.conv_tol = p["scf"]["conv_tol_Eh"]
    mf.conv_tol_grad = p["scf"]["conv_tol_grad"]
    mf.direct_scf_tol = p["scf"]["direct_scf_tol"]
    mf.max_cycle = p["scf"]["max_cycle"]
    mf.chkfile = None
    mf.kernel()
    assert mf.converged, f"RHF did not converge in {tag}"
    model = mp.MP2(mf, frozen=0)
    assert model.__class__.__module__=="pyscf.mp.mp2"
    model.kernel()
    grad = model.nuc_grad_method().kernel()
    assert np.isfinite(grad).all() and np.isfinite(model.e_tot)
    result = {"id":tag,"stage":stage,"completed_at_utc":now(),"protocol_sha256":ph,
              "geometry_bohr":xyz.tolist(),"energy_Eh":float(model.e_tot),
              "RHF_energy_Eh":float(mf.e_tot),"MP2_correlation_energy_Eh":float(model.e_corr),
              "gradient_Eh_per_bohr":grad.tolist(),"RHF_converged":bool(mf.converged),
              "basis_functions":mol.nao_nr(),"occupied_orbitals":int((mf.mo_occ>0).sum()),
              "MO_gap_Eh":float(mf.mo_energy[5]-mf.mo_energy[4]),
              "elapsed_seconds":time.monotonic()-t}
    path = ROOT/"evaluations"/(tag+".json")
    path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    records.append({"id":tag,"file":path.relative_to(ROOT).as_posix(),"sha256":sha(path)})
    cache[key]=result
    print(json.dumps({"completed":tag,"stage":stage,"seconds":round(result["elapsed_seconds"],3)}),flush=True)
    return result

optimizations = []
try:
    for index, initial in enumerate(p["geometry_optimization"]["starts"]):
        q0 = np.array([initial["r_angstrom"]/BOHR_A,np.deg2rad(initial["angle_deg"])])
        history=[]
        def objective(q):
            r,t=q
            assert .5/BOHR_A<r<1.5/BOHR_A and .5<t<3., "Optimizer left water search region."
            result=evaluate(coords(q),f"geometry-start-{index}")
            g=np.array(result["gradient_Eh_per_bohr"])
            dr=np.array([[0,0,0],[np.sin(t/2),0,np.cos(t/2)],[-np.sin(t/2),0,np.cos(t/2)]])
            dt=np.array([[0,0,0],[r/2*np.cos(t/2),0,-r/2*np.sin(t/2)],[-r/2*np.cos(t/2),0,-r/2*np.sin(t/2)]])
            dq=np.array([np.sum(g*dr),np.sum(g*dt)])
            history.append({"q":q.tolist(),"evaluation":result["id"],"internal_gradient":dq.tolist()})
            return result["energy_Eh"],dq
        opt=minimize(objective,q0,jac=True,method="BFGS",
                     options={"gtol":p["geometry_optimization"]["internal_gradient_gtol"],
                              "maxiter":p["geometry_optimization"]["max_iterations_per_start"]})
        final=evaluate(coords(opt.x),f"geometry-start-{index}-final")
        row={"start":index,"q":opt.x.tolist(),"r_angstrom":float(opt.x[0]*BOHR_A),
             "angle_deg":float(np.rad2deg(opt.x[1])),"geometry_bohr":coords(opt.x).tolist(),
             "optimizer_success":bool(opt.success),"optimizer_message":str(opt.message),
             "optimizer_iterations":int(opt.nit),"energy_Eh":final["energy_Eh"],
             "max_cartesian_gradient_Eh_per_bohr":float(np.max(abs(np.array(final["gradient_Eh_per_bohr"])))),
             "final_evaluation":final["id"],"history":history}
        optimizations.append(row)
        write("optimizations.json",optimizations)
        assert row["max_cartesian_gradient_Eh_per_bohr"]<=p["numerical_validity"]["maximum_cartesian_gradient_Eh_per_bohr"], "Geometry is not stationary enough."
        print(json.dumps({"geometry_finished":index,"r_angstrom":row["r_angstrom"],"angle_deg":row["angle_deg"],"max_gradient":row["max_cartesian_gradient_Eh_per_bohr"]}),flush=True)
    xyz=np.array(optimizations[0]["geometry_bohr"])
    hessians=[]
    for h in p["hessian"]["steps_bohr"]:
        hess=np.zeros((9,9))
        refs=[]
        for axis in range(9):
            delta=np.zeros(9);delta[axis]=h
            plus=evaluate(xyz+delta.reshape(3,3),f"hessian-{h}-axis-{axis}-plus")
            minus=evaluate(xyz-delta.reshape(3,3),f"hessian-{h}-axis-{axis}-minus")
            hess[:,axis]=(np.array(plus["gradient_Eh_per_bohr"]).reshape(-1)-np.array(minus["gradient_Eh_per_bohr"]).reshape(-1))/(2*h)
            refs.append({"axis":axis,"plus":plus["id"],"minus":minus["id"]})
        hessians.append({"step_bohr":h,"raw_hessian_Eh_per_bohr2":hess.tolist(),
                         "symmetrized_hessian_Eh_per_bohr2":((hess+hess.T)/2).tolist(),
                         "maximum_antisymmetry_Eh_per_bohr2":float(np.max(abs(hess-hess.T))),
                         "evaluation_references":refs})
        write("hessians.json",hessians)
        print(json.dumps({"hessian_complete":h,"max_antisymmetry":hessians[-1]["maximum_antisymmetry_Eh_per_bohr2"]}),flush=True)
    status="completed"
except BaseException as exc:
    status="incomplete"
    write("execution-error.json",{"type":type(exc).__name__,"message":str(exc),"at_utc":now(),"completed_energy_gradient_evaluations":len(records)})
    raise
finally:
    write("execution-receipt.json",{"status":status,"started_at_utc":start,"finished_at_utc":now(),
         "protocol_sha256":ph,"runner_sha256":sha(Path(__file__)),"evaluations":records,
         "energy_gradient_evaluation_count":len(records),"versions":versions,"elapsed_seconds":time.monotonic()-began,
         "published_xTB_files_written":False,"alternative_method_tuned":False})
