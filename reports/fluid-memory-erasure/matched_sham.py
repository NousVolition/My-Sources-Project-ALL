"""Secondary equal-size structured displacement; match using initial fields only."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from solver import Solver
from run import ROOT,P,dump,digest,filename,DIAG_NAMES
from analyze import branch_index

DESIGN=json.loads((ROOT/'matched-sham-protocol.json').read_text())


def norm(a,modes):
    return np.sqrt(np.sum(abs(a)**2*np.where(modes[:,2]==0,1.,2.),axis=(-2,-1)))


def match_translation(retained,scrambles,modes):
    high=np.max(abs(modes),axis=1)>P['coarse_band_max']
    r=retained[:, :, high];ks=modes[high]
    target=float(np.mean(norm(scrambles[:,:,:,high]-r[None],ks)))
    direction=np.array([1.,-.7,.9]);direction/=np.linalg.norm(direction)
    angles=ks@direction
    def change(length):return float(np.mean(norm(r*(np.exp(1j*length*angles)-1)[None,None],ks)))
    samples=np.linspace(0,4*np.pi,513)
    crossing=next((b for a,b in zip(samples[:-1],samples[1:]) if change(a)<=target<=change(b)),None)
    if crossing is None:raise RuntimeError('No prespecified matched translation exists')
    left=crossing-(samples[1]-samples[0])
    root=brentq(lambda length:change(length)-target,left,crossing,xtol=1e-13)
    phase=np.ones(len(modes),complex);phase[high]=np.exp(1j*root*angles)
    transformed=retained*phase[None,None]
    return transformed,dict(length=float(root),direction=direction.tolist(),target_change_norm=target,
                             achieved_change_norm=change(root),absolute_matching_error=abs(change(root)-target),
                             high_band_norm=float(np.mean(norm(r,ks))))


def run_one(seed,n,dt,data,output,backend='gpu'):
    source_stem=filename(seed,n,dt,P['nu'],P['probe_velocity_rms'],max(P['horizons']),backend)
    stem=source_stem+'-matched-sham';raw=output/(stem+'.npz');meta_path=raw.with_suffix('.json')
    sources={f:digest(ROOT/f) for f in ['solver.py','matched_sham.py','matched-sham-protocol.json']}
    source_raw=data/(source_stem+'.npz');source_meta=json.loads((data/(source_stem+'.json')).read_text())
    if digest(source_raw)!=source_meta['raw_sha256']:raise RuntimeError('Source data hash mismatch')
    if meta_path.exists():
        old=json.loads(meta_path.read_text())
        if old['sources']!=sources or old['source_raw_sha256']!=source_meta['raw_sha256'] or digest(raw)!=old['raw_sha256']:
            raise RuntimeError('Matched-sham cache differs')
        print(stem+' verified cache',flush=True);return old
    tic=time.time()
    # Only initial branch coefficients and initial coordinates are read as inputs.
    with np.load(source_raw,allow_pickle=False) as z:
        modes=z['modes'];start=z['start_fields'];xstart=z['start_positions'];probe=z['probe']
    retained=np.stack([start[branch_index(source_meta,'retained',h,0)] for h in range(2)])
    scrambles=np.stack([np.stack([start[branch_index(source_meta,'scramble'+str(r),h,0)] for h in range(2)]) for r in range(P['phase_scramble_replicates'])])
    shifted,matching=match_translation(retained,scrambles,modes)
    s=Solver(n,P['nu'],backend);xp=s.xp
    coefficients=np.stack([shifted[h]+kick*probe for h in range(2) for kick in range(2)])
    h=s.unpack(coefficients)
    x=xp.asarray(np.stack([xstart[branch_index(source_meta,'retained',hist,0)] for hist in range(2) for _ in range(2)]))
    f=xp.asarray(np.tile(np.eye(3),(4,P['centers'],1,1)));e0=.5*s.inner(h,h);budget=xp.zeros(4);cflmax=xp.zeros(4)
    t=0.;times=[];positions=[];tangents=[];diagnostics=[];response_time=[];response_fields=[]
    for target in np.arange(0,DESIGN['end']+1e-9,P['save_dt']):
        while t<target-1e-12:
            step=min(dt,target-t);h,x,f,db,cfl=s.step(h,x,f,step);budget+=db;cflmax=xp.maximum(cflmax,cfl);t+=step
        times.append(float(target));positions.append(s.cpu(x));tangents.append(s.cpu(f));diagnostics.append(s.diagnostics(h,x,f,budget,e0,cflmax))
        if any(abs(target-horizon)<1e-8 for horizon in DESIGN['target_horizons']):
            response_time.append(float(target));response_fields.append(s.pack(h[1::2]-h[::2]))
    np.savez_compressed(raw,modes=modes,start_fields=coefficients,time=times,positions=positions,tangent=tangents,diagnostics=diagnostics,
                        end_fields=s.pack(h),response_time=response_time,response_fields=response_fields)
    maxima={name:float(np.max(abs(np.array(diagnostics)[...,i]))) if i==1 else float(np.max(np.array(diagnostics)[...,i])) for i,name in enumerate(DIAG_NAMES)}
    meta=dict(config=dict(seed=seed,n=n,dt=dt,backend=backend),status='complete',sources=sources,source_raw_sha256=source_meta['raw_sha256'],
              raw_sha256=digest(raw),matching=matching,maxima=maxima,diagnostic_columns=DIAG_NAMES,seconds=time.time()-tic,
              per_mode_power_error=float(abs(np.sum(abs(shifted)**2,axis=1)-np.sum(abs(retained)**2,axis=1)).max()),
              coarse_coefficient_error=float(abs((shifted-retained)[:,:,np.max(abs(modes),axis=1)<=2]).max()))
    dump(meta_path,meta);print(f'{stem} complete {meta["seconds"]:.1f}s; matching error {matching["absolute_matching_error"]:.3g}',flush=True)
    return meta


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');ap.add_argument('--out',type=Path,default=ROOT/'matched-sham-data')
    ap.add_argument('--backend',choices=['cpu','gpu'],default='gpu');args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    jobs=[(seed,n,P['dt']) for seed in DESIGN['seeds'] for n in DESIGN['grids']]+[(seed,DESIGN['half_step_n'],P['dt']/2) for seed in DESIGN['half_step_seeds']]
    completed=[]
    for seed,n,dt in jobs:
        result=run_one(seed,n,dt,args.data,args.out,args.backend);completed.append(result['config'])
        dump(args.out/'execution.json',dict(planned=len(jobs),completed=len(completed),configs=completed,design_sha256=digest(ROOT/'matched-sham-protocol.json')))
