import json,time,numpy as np
from final_experiment import ROOT,P,DrivenSolver,digest,dump,DIAG
from analysis import read,effective_regimes

def run_one(seed,regime,n,dt=None):
    dt=effective_regimes()[regime]['dt'] if dt is None else dt
    source=ROOT/'common-probe-data'/f'{regime}-coincident-s{seed}-n{n}-dt{dt:g}.npz';a,m=read(source)
    folder=ROOT/'common-unforced-data';folder.mkdir(exist_ok=True);stem=f'{regime}-s{seed}-n{n}-dt{dt:g}'
    raw=folder/(stem+'.npz');mp=raw.with_suffix('.json');src={f:digest(ROOT/f) for f in ['safe_common_unforced.py','unforced-readout-protocol.json','final_experiment.py','corrected_solver.py','corrected_observation.py']}
    if mp.exists():
        old=json.loads(mp.read_text())
        if old['sources']!=src or old['parent_raw_sha256']!=m['raw_sha256'] or old['raw_sha256']!=digest(raw):raise RuntimeError('Common unforced cache mismatch')
        return old
    tic=time.time();s=DrivenSolver(n,P['regimes'][regime]['nu'],0,'gpu');xp=s.xp
    starts=np.stack([h for j in range(len(m['branches'])) for h in [a['base_start_fields'][j],a['start_fields'][j]]])
    labels=[dict(**b,kicked=bool(k)) for b in m['branches'] for k in range(2)];times=a['time'];pieces=[];ends=[]
    batch=8 if n>=128 else 16
    for off in range(0,len(starts),batch):
        h=s.unpack(starts[off:off+batch]);b=len(h);x=xp.asarray(np.stack([a['positions'][0,(off+j)//2] for j in range(b)]));f=xp.asarray(np.tile(np.eye(3),(b,P['centers'],1,1)))
        e0=.5*s.inner(h,h);budget=xp.zeros(b);cfl=xp.zeros(b);t=0.;xx=[];ff=[];dd=[]
        for target in times:
            while t<target-1e-12:
                step=min(dt,target-t);h,x,f,db,cf=s.step(h,x,f,step);budget+=db;cfl=xp.maximum(cfl,cf);t+=step
            xx.append(s.cpu(x));ff.append(s.cpu(f));dd.append(s.diagnostics(h,x,f,budget,e0,cfl))
        pieces.append((xx,ff,dd));ends.append(s.pack(h))
    positions=np.concatenate([r[0] for r in pieces],axis=1);tangent=np.concatenate([r[1] for r in pieces],axis=1);diagnostics=np.concatenate([r[2] for r in pieces],axis=1)
    np.savez_compressed(raw,modes=a['modes'],pulse=a['pulse'],start_fields=starts,end_fields=np.concatenate(ends),time=times,positions=positions,tangent=tangent,diagnostics=diagnostics)
    maxima={name:float(np.max(abs(diagnostics[...,j]))) if j==1 else float(np.max(diagnostics[...,j])) for j,name in enumerate(DIAG)}
    result=dict(config=dict(seed=seed,regime=regime,n=n,dt=dt,injection=0),branches=labels,sources=src,parent_raw_sha256=m['raw_sha256'],raw_sha256=digest(raw),maxima=maxima,seconds=time.time()-tic)
    dump(mp,result);print('common unforced '+stem+f' {result["seconds"]:.1f}s',flush=True);xp.get_default_memory_pool().free_all_blocks();return result
