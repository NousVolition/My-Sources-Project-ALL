"""Small independent ensemble using the existing matched-stretch NS solver."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from numerics import Flow, L


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sample(u, positions):
    """Periodic trilinear u AND exact spatial Jacobian of that interpolant.

    Gradient convention J[i,j]=du_i/dx_j. Coordinates are unwrapped [-3,3).
    Cell faces have a one-sided gradient; measure-zero ties use floor().
    """
    n = u.shape[1]
    q = (positions + L/2) % L * n/L
    a = np.floor(q).astype(int)
    b = q-a
    vel = np.zeros((len(q), 3))
    jac = np.zeros((len(q), 3, 3))
    for bits in np.ndindex(2,2,2):
        s = np.array(bits)
        ix = (a+s) % n
        f = np.where(s, b, 1-b)
        values = u[:,ix[:,0],ix[:,1],ix[:,2]].T
        vel += values * np.prod(f,axis=1)[:,None]
        for j in range(3):
            weight = (2*s[j]-1)*n/L*np.prod(f[:,[k for k in range(3) if k != j]],axis=1)
            jac[:,:,j] += values*weight[:,None]
    return vel, jac


def initial(flow, seed):
    """Grid-independent smooth random Fourier polynomial; normalization by energy."""
    rng = np.random.default_rng(seed)
    axis = np.arange(flow.n)*flow.dx-L/2
    xyz = np.meshgrid(axis,axis,axis,indexing="ij",sparse=True)
    modes = [(1,1,0),(1,0,2),(0,2,1),(2,-1,1),(1,2,-2),(2,2,1),
             (1,0,0),(0,1,0),(0,0,1),(1,-1,1),(3,1,0),(0,1,3)]
    raw = np.zeros((3,flow.n,flow.n,flow.n))
    for mode in modes:
        phase = sum(2*np.pi/L*m*x for m,x in zip(mode,xyz))+rng.uniform(0,2*np.pi)
        raw += rng.normal(size=3)[:,None,None,None]*np.cos(phase)/np.linalg.norm(mode)
    h = flow.project(flow.hat(raw))
    h[:,0,0,0] = 0
    h /= np.sqrt(flow.inner(h,h)/L**3)
    return h


def run(seed, config, destination):
    destination = Path(destination)
    destination.mkdir(parents=True,exist_ok=True)
    name = f"s{seed}_n{config['n']}_nu{config['nu']}_dt{config['dt']}"
    path = destination/(name+".npz")
    meta = destination/(name+".json")
    source = {p.name: sha(p) for p in [Path(__file__),ROOT.parent/"numerics.py",ROOT/"protocol.json"]}
    if path.exists() and meta.exists():
        old=json.loads(meta.read_text())
        if old["sources"] != source or old["config"] != config or old["sha256"] != sha(path):
            raise ValueError("Stale or modified simulation: " + str(path))
        return path
    f = Flow(config["n"],config["nu"],workers=1)
    h = initial(f,seed)
    # Random labels have no numerical or physical role. Same seeding across variants.
    x = np.random.default_rng(seed+80000).uniform(-L/2,L/2,(config["markers"],3))
    M = np.broadcast_to(np.eye(3),(len(x),3,3)).copy()
    dt = config["dt"]
    nstep = round(config["end"]/dt)
    stride = round(config["save_dt"]/dt)
    assert np.isclose(nstep*dt,config["end"]) and np.isclose(stride*dt,config["save_dt"])
    arrays = {k:[] for k in ["time","positions","velocity","gradient","tangent"]}
    diag=[]
    start=time.perf_counter()
    energy0=f.inner(h,h)/2
    max_cfl=0.
    integrated_diss=0.
    for step in range(nstep+1):
        u=f.real(h)
        vel,J=sample(u,x)
        if step % stride == 0:
            for key,value in zip(arrays,[step*dt,x.copy(),vel,J,M.copy()]):
                arrays[key].append(value)
            wh=f.curl(h)
            divergence=f.real(1j*sum(k*c for k,c in zip(f.k,h)))
            energy=f.inner(h,h)/2
            diag.append({"t":step*dt,"energy":energy,"energy_budget_residual":(energy-energy0+integrated_diss)/energy0,
                         "divergence_max":float(abs(divergence).max()),
                         "high_band_energy_fraction":f.inner(h*f.high,h*f.high)/max(2*energy,1e-30),
                         "max_tangent_condition":float(np.linalg.cond(M).max()),
                         "max_det_error":float(abs(np.linalg.det(M)-1).max())})
        if step == nstep: break
        # Identical Heun fluid update to Flow.step; joint marker/tangent Heun stages.
        r1,a,speed1=f.rhs(h)
        hp=h+dt*r1
        xp=x+dt*vel
        dM=J@M
        Mp=M+dt*dM
        vp,Jp=sample(f.real(hp),xp)
        r2,b,speed2=f.rhs(hp)
        h=f.project(h+dt/2*(r1+r2))
        x=x+dt/2*(vel+vp)
        M=M+dt/2*(dM+Jp@Mp)
        integrated_diss+=dt/2*(a[1]+b[1])
        max_cfl=max(max_cfl,max(speed1,speed2)*dt/f.dx)
        if max_cfl > .5 or not np.isfinite(M).all() or not np.isfinite(h).all():
            raise FloatingPointError("Pilot numerical safety check failed")
    np.savez_compressed(path,**arrays)
    save_json(meta,{"name":name,"independent_group":seed,"seed":seed,"config":config,
                    "sources":source,"sha256":sha(path),"seconds":time.perf_counter()-start,
                    "max_cfl":max_cfl,"diagnostics":diag,"Re_Urms_L_over_nu":L/config["nu"]})
    print(name+f" complete ({time.perf_counter()-start:.1f}s)",flush=True)
    return path


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",type=Path,default=ROOT/"data")
    ap.add_argument("--pilot",action="store_true")
    ap.add_argument("--sensitivity",action="store_true")
    args=ap.parse_args()
    p=json.loads((ROOT/"protocol.json").read_text())
    seeds=p["train_seeds"]+p["validation_seeds"]+p["test_seeds"]
    for seed in (seeds[:2] if args.pilot else seeds): run(seed,p["base"],args.data)
    if args.sensitivity:
        for seed in p["refinement_seeds"]:
            for change in [{"n":38},{"dt":.0025},{"nu":.01},{"nu":.04}]:
                run(seed,{**p["base"],**change},args.data)


if __name__ == "__main__": main()
