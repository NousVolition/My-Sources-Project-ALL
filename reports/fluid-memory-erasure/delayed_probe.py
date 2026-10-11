"""Disturb the fluid after a waiting interval; marker resets are observations only."""
import argparse,json,time
from pathlib import Path
import numpy as np
from solver import Solver
from run import ROOT,P,dump,digest,filename,DIAG_NAMES

DESIGN=json.loads((ROOT/'delayed-probe-protocol.json').read_text())


def run_one(seed,n,dt,wait,data,output,backend='gpu'):
    source_stem=filename(seed,n,dt,P['nu'],P['probe_velocity_rms'],max(P['horizons']),backend)
    stem=source_stem+f'-wait{wait:g}';rp=output/(stem+'.npz');mp=rp.with_suffix('.json')
    sources={f:digest(ROOT/f) for f in ['solver.py','delayed_probe.py','delayed-probe-protocol.json']}
    source=data/(source_stem+'.npz');sm=json.loads(source.with_suffix('.json').read_text())
    if digest(source)!=sm['raw_sha256']:raise RuntimeError('Bad source checkpoint')
    if mp.exists():
        old=json.loads(mp.read_text())
        if old['sources']!=sources or old['source_raw_sha256']!=sm['raw_sha256'] or digest(rp)!=old['raw_sha256']:raise RuntimeError('Changed delay cache')
        print(stem+' verified cache',flush=True);return old
    tic=time.time();s=Solver(n,P['nu'],backend);xp=s.xp
    with np.load(source,allow_pickle=False) as z:
        initial=z['start_fields'][::2];ending=z['end_fields'][::2];probe=z['probe'];positions=z['start_positions'][::2]
    waiting_diag=[]
    if np.isclose(wait,.8):
        waited=ending;method='Exact stored unprobed endpoint at wait=0.8'
    else:
        h=s.unpack(initial);x=xp.asarray(positions);f=xp.asarray(np.tile(np.eye(3),(len(initial),P['centers'],1,1)))
        e0=.5*s.inner(h,h);budget=xp.zeros(len(h));cflmax=xp.zeros(len(h));t=0.
        while t<wait-1e-12:
            step=min(dt,wait-t);h,x,f,db,cfl=s.step(h,x,f,step);budget+=db;cflmax=xp.maximum(cflmax,cfl);t+=step
        waiting_diag=s.diagnostics(h,x,f,budget,e0,cflmax);waited=s.pack(h);method='Unforced continuation from stored initial fields'
    coefficients=np.stack([field+kicked*probe for field in waited for kicked in [0,1]])
    h=s.unpack(coefficients)
    # Same centers in every causal arm/history; native branches retain their own
    # original starting observation positions and remain a separate diagnostic.
    x=xp.asarray(np.repeat(positions,2,axis=0));f=xp.asarray(np.tile(np.eye(3),(len(h),P['centers'],1,1)))
    e0=.5*s.inner(h,h);budget=xp.zeros(len(h));cflmax=xp.zeros(len(h));t=0.
    times=[];xs=[];fs=[];diags=[];rt=[];rf=[]
    for target in np.arange(0,DESIGN['end']+1e-9,P['save_dt']):
        while t<target-1e-12:
            step=min(dt,target-t);h,x,f,db,cfl=s.step(h,x,f,step);budget+=db;cflmax=xp.maximum(cflmax,cfl);t+=step
        times.append(float(target));xs.append(s.cpu(x));fs.append(s.cpu(f));diags.append(s.diagnostics(h,x,f,budget,e0,cflmax))
        if any(abs(target-horizon)<1e-8 for horizon in P['horizons']):rt.append(float(target));rf.append(s.pack(h[1::2]-h[::2]))
    np.savez_compressed(rp,modes=s.modes_cpu,waited_fields=waited,start_fields=coefficients,time=times,positions=xs,tangent=fs,
                        diagnostics=diags,waiting_diagnostics=waiting_diag,response_time=rt,response_fields=rf,end_fields=s.pack(h))
    maxima={name:float(np.max(abs(np.array(diags)[...,i]))) if i==1 else float(np.max(np.array(diags)[...,i])) for i,name in enumerate(DIAG_NAMES)}
    if len(waiting_diag):
        for i,name in enumerate(DIAG_NAMES):maxima[name]=max(maxima[name],float(np.max(abs(waiting_diag[...,i]))))
    meta=dict(config=dict(seed=seed,n=n,dt=dt,wait=wait,backend=backend,amplitude=P['probe_velocity_rms']),branches=sm['branches'],sources=sources,
              source_raw_sha256=sm['raw_sha256'],raw_sha256=digest(rp),maxima=maxima,diagnostic_columns=DIAG_NAMES,
              checkpoint_method=method,seconds=time.time()-tic)
    dump(mp,meta);print(f'{stem} complete {meta["seconds"]:.1f}s',flush=True);return meta


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');ap.add_argument('--out',type=Path,default=ROOT/'delayed-probe-data');ap.add_argument('--backend',choices=['cpu','gpu'],default='gpu');args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    jobs=[(s,DESIGN['n'],DESIGN['dt'],w) for s in DESIGN['seeds'] for w in DESIGN['waits']]
    jobs += [(s,DESIGN['fine_n'],DESIGN['dt'],.8) for s in DESIGN['fine_seeds']]
    jobs += [(s,DESIGN['n'],DESIGN['dt']/2,.8) for s in DESIGN['half_step_seeds']]
    completed=[]
    for s,n,dt,w in jobs:
        m=run_one(s,n,dt,w,args.data,args.out,args.backend);completed.append(m['config'])
        dump(args.out/'execution.json',dict(planned=len(jobs),completed=len(completed),configs=completed,protocol_sha256=digest(ROOT/'delayed-probe-protocol.json')))
