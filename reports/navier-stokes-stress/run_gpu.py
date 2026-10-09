"""Run the unchanged planned matrix with the separately validated GPU backend.
CPU results remain in data/; GPU results remain in data-gpu/.
"""
from pathlib import Path
import argparse,json,platform,sys,time
import cupy as cp
import numpy as np
from gpu_solver import GPUSolver
from run import ROOT,tasks,name,dump,digest

def run_gpu(t):
    folder=ROOT/'data-gpu'/name(t);folder.mkdir(parents=True,exist_ok=True)
    if (folder/'result.json').exists():return name(t)+' cached'
    tic=time.time();s=GPUSolver(t['n'],t['nu']);h=s.initial(t['family'],t['epsilon'],t['seed'])
    e0=.5*s.inner(h,h);clock=0.;loss=0.;steps=0;cflmax=0;rows=[];spectra=[];slices=[];saved=[];status='complete';error=None
    sources={p.name:digest(p) for p in [ROOT/'solver.py',ROOT/'run.py',ROOT/'protocol.json',ROOT/'gpu_solver.py',ROOT/'run_gpu.py',ROOT/'gpu-amendment.json']}
    dump(folder/'settings.json',dict(**t,backend='cupy',source_sha256=sources))
    maxk=float(s.k2[s.keep].max())
    try:
        for target in np.linspace(0,t['end'],round(t['end']/.1)+1):
            while clock<target-1e-12:
                u=s.real(h);umax=float(cp.sum(abs(u),axis=0).max())
                dt=min(t['dt'],.35*s.dx/max(umax,1e-20),.4/(t['nu']*maxk),target-clock)
                h,dl,cfl=s.step(h,dt);loss+=dl;clock+=dt;steps+=1;cflmax=max(cflmax,cfl)
                if not bool(cp.isfinite(h).all()):raise FloatingPointError('Nonfinite field')
                if cfl>1:raise FloatingPointError('Gross CFL failure')
            row,sp,sl=s.observe(h,float(target),loss,e0);row['max_cfl_so_far']=cflmax
            row['local_screen_pass']=bool(row['tail_Z']<.01 and abs(row['energy_residual'])<.001 and row['divergence_rms']<1e-9 and cflmax<=.5)
            rows.append(row);spectra.append(sp);slices.append(sl)
            interval=.5 if t['family']=='kida' else 1.
            if abs(target/interval-round(target/interval))<1e-9:
                fn=f'field-{target:.2f}.npz';np.savez_compressed(folder/fn,h=cp.asnumpy(h),t=target)
                saved.append(dict(file=fn,sha256=digest(folder/fn),t=float(target)))
            dump(folder/'progress.json',dict(t=float(target),steps=steps,seconds=time.time()-tic,tail_Z=row['tail_Z'],local_screen_pass=row['local_screen_pass']))
            if abs(row['energy_residual'])>.05:raise FloatingPointError('Energy budget residual exceeds 5 percent')
    except Exception as exc:status='failed';error=repr(exc)
    np.savez_compressed(folder/'diagnostics.npz',spectrum=np.array(spectra),vorticity_slice=np.array(slices),times=[r['t'] for r in rows])
    dump(folder/'result.json',dict(settings=t,backend='cupy',status=status,error=error,rows=rows,fields=saved,diagnostics_sha256=digest(folder/'diagnostics.npz'),steps=steps,seconds=time.time()-tic,source_sha256=sources))
    del h,s;cp.get_default_memory_pool().free_all_blocks()
    return f'{name(t)} {status}, t={clock:.2f}, {time.time()-tic:.1f}s, tail={rows[-1]["tail_Z"]:.4g}'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int);p.add_argument('--only');a=p.parse_args()
    check=json.loads((ROOT/'gpu-validation.json').read_text());assert check['passed'],'GPU verification must pass first'
    data=ROOT/'data-gpu';data.mkdir(exist_ok=True)
    dump(data/'environment.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,cupy=cp.__version__,gpu=check['gpu'],started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),argv=sys.argv))
    jobs=[t for t in tasks() if not a.only or a.only in name(t)]
    if a.limit:jobs=jobs[:a.limit]
    for t in jobs:print(run_gpu(t),flush=True)
