import json,time
import numpy as np
from final_experiment import ROOT,P,DrivenSolver,digest,dump,DIAG
from analysis import read,filename

def run_one(seed,regime,n,dt=None):
    cfg=P['regimes'][regime];dt=cfg['dt'] if dt is None else dt
    source=filename(seed,regime,'coincident',n,dt);a,m=read(source,['modes','start_fields','start_positions','probe_basis','pulse_coefficients','prepared_fields'])
    folder=ROOT/'unforced-readout-data';folder.mkdir(exist_ok=True);stem=f'{regime}-s{seed}-n{n}-dt{dt:g}'
    raw=folder/(stem+'.npz');mp=raw.with_suffix('.json');src={f:digest(ROOT/f) for f in ['safe_unforced_readout.py','unforced-readout-protocol.json','corrected_solver.py','corrected_observation.py','final_experiment.py']}
    if mp.exists():
        old=json.loads(mp.read_text())
        if old['sources']!=src or old['parent_raw_sha256']!=m['raw_sha256'] or old['raw_sha256']!=digest(raw):raise RuntimeError('Unforced control cache mismatch')
        return old
    tic=time.time();s=DrivenSolver(n,cfg['nu'],0,'gpu');xp=s.xp;times=np.arange(0,P['horizon']+1e-9,P['save_dt']);pieces=[];ends=[]
    batch=4 if n>=128 else 8
    for off in range(0,len(a['start_fields']),batch):
        h=s.unpack(a['start_fields'][off:off+batch]);b=len(h);x=xp.asarray(np.tile(a['start_positions'],(b,1,1)));f=xp.asarray(np.tile(np.eye(3),(b,P['centers'],1,1)))
        e0=.5*s.inner(h,h);budget=xp.zeros(b);cfl=xp.zeros(b);t=0.;xx=[];ff=[];dd=[]
        for target in times:
            while t<target-1e-12:
                step=min(dt,target-t);h,x,f,db,cf=s.step(h,x,f,step);budget+=db;cfl=xp.maximum(cfl,cf);t+=step
            xx.append(s.cpu(x));ff.append(s.cpu(f));dd.append(s.diagnostics(h,x,f,budget,e0,cfl))
        pieces.append((xx,ff,dd));ends.append(s.pack(h))
    positions=np.concatenate([r[0] for r in pieces],axis=1);tangent=np.concatenate([r[1] for r in pieces],axis=1);diagnostics=np.concatenate([r[2] for r in pieces],axis=1)
    np.savez_compressed(raw,**a,time=times,positions=positions,tangent=tangent,diagnostics=diagnostics,end_fields=np.concatenate(ends))
    maxima={name:float(np.max(abs(diagnostics[...,j]))) if j==1 else float(np.max(diagnostics[...,j])) for j,name in enumerate(DIAG)}
    result=dict(config=dict(seed=seed,regime=regime,n=n,dt=dt,injection=0),branches=m['branches'],sources=src,parent_raw_sha256=m['raw_sha256'],raw_sha256=digest(raw),maxima=maxima,seconds=time.time()-tic)
    dump(mp,result);print('unforced readout '+stem+f' {result["seconds"]:.1f}s',flush=True);xp.get_default_memory_pool().free_all_blocks();return result
