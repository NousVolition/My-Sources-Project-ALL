from pathlib import Path
import argparse, concurrent.futures, hashlib, json, platform, sys, time
import numpy as np
from solver import Solver

ROOT=Path(__file__).resolve().parent
def dump(p,x): p.write_text(json.dumps(x,indent=2),encoding='utf-8')
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def name(t):return f"{t['family']}-n{t['n']}-nu{t['nu']:g}-e{t['epsilon']:g}-s{t['seed']}-dt{t['dt']:g}"
def tasks():
    out=[]
    def add(f,n,nu,eps=0,seed=1,dt=.01):out.append(dict(family=f,n=n,nu=nu,epsilon=eps,seed=seed,dt=dt,end=2. if f=='kida' else 4.))
    for n in [32,48,64,96]:
        for nu in [.02,.005,.001]:add('kida',n,nu)
    for nu in [.02,.001]:add('kida',128,nu)
    for n in [32,64,96]:
        for nu in [.01,.001]:add('taylor_green',n,nu)
    for seed in [1,2,3,4]:add('kida',48,.005,.001,seed)
    for n in [64,96]:
        for eps in [1e-6,.001,.01]:add('kida',n,.005,eps)
    for nu,eps in [(.005,0),(.001,0),(.005,.001)]:add('kida',64,nu,eps,dt=.005)
    return out

def run(t):
    folder=ROOT/'data'/name(t);folder.mkdir(parents=True,exist_ok=True)
    if (folder/'result.json').exists():return name(t)+' cached'
    tic=time.time();s=Solver(t['n'],t['nu']);h=s.initial(t['family'],t['epsilon'],t['seed'])
    e0=.5*s.inner(h,h);clock=0.;loss=0.;steps=0;cflmax=0;rows=[];spectra=[];slices=[];saved=[];status='complete';error=None
    sources={p.name:digest(p) for p in [ROOT/'solver.py',ROOT/'run.py',ROOT/'protocol.json']}
    # Every run stores original settings and source hashes before evolving.
    dump(folder/'settings.json',dict(**t,source_sha256=sources))
    try:
        for target in np.linspace(0,t['end'],round(t['end']/.1)+1):
            while clock < target-1e-12:
                u=s.real(h);umax=float(np.sum(abs(u),axis=0).max())
                dt=min(t['dt'],.35*s.dx/max(umax,1e-20),.4/(t['nu']*float(s.k2[s.keep].max())),target-clock)
                h,dl,cfl=s.step(h,dt);loss+=dl;clock+=dt;steps+=1;cflmax=max(cflmax,cfl)
                if not np.isfinite(h).all():raise FloatingPointError('Nonfinite field')
                if cfl>1:raise FloatingPointError('Gross CFL failure')
            row,sp,sl=s.observe(h,float(target),loss,e0);row['max_cfl_so_far']=cflmax
            row['local_screen_pass']=bool(row['tail_Z']<.01 and abs(row['energy_residual'])<.001 and row['divergence_rms']<1e-9 and cflmax<=.5)
            rows.append(row);spectra.append(sp);slices.append(sl)
            interval=.5 if t['family']=='kida' else 1.
            if abs(target/interval-round(target/interval))<1e-9:
                fn=f"field-{target:.2f}.npz";np.savez_compressed(folder/fn,h=h,t=target)
                saved.append(dict(file=fn,sha256=digest(folder/fn),t=float(target)))
            dump(folder/'progress.json',dict(t=float(target),steps=steps,seconds=time.time()-tic,tail_Z=row['tail_Z'],local_screen_pass=row['local_screen_pass']))
            if abs(row['energy_residual'])>.05:raise FloatingPointError('Energy budget residual exceeds 5 percent')
    except Exception as exc:status='failed';error=repr(exc)
    np.savez_compressed(folder/'diagnostics.npz',spectrum=np.array(spectra),vorticity_slice=np.array(slices),times=[r['t'] for r in rows])
    result=dict(settings=t,status=status,error=error,rows=rows,fields=saved,diagnostics_sha256=digest(folder/'diagnostics.npz'),steps=steps,seconds=time.time()-tic,source_sha256=sources)
    dump(folder/'result.json',result)
    return f"{name(t)} {status}, t={clock:.2f}, {result['seconds']:.1f}s, tail={rows[-1]['tail_Z']:.4g}"

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);p.add_argument('--max-n',type=int,default=128);p.add_argument('--only');p.add_argument('--limit',type=int);a=p.parse_args()
    jobs=[t for t in tasks() if t['n']<=a.max_n and (not a.only or a.only in name(t))]
    if a.limit:jobs=jobs[:a.limit]
    (ROOT/'data').mkdir(exist_ok=True)
    dump(ROOT/'planned_tasks.json',tasks());dump(ROOT/'data'/'environment.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),argv=sys.argv))
    with concurrent.futures.ProcessPoolExecutor(a.workers) as pool:
        for msg in pool.map(run,jobs): print(msg,flush=True)
