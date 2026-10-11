"""Compare saved initial states with a separate CPU evolution implementation."""
import argparse,json,time
import numpy as np
from threadpoolctl import threadpool_limits
from independent_solver import IndependentSolver
from accelerated_experiment import ROOT,P,DrivenSolver,digest,dump
from analysis import read,filename,idx

def check(seed,regime='smooth',tight=False,pilot=False):
    n=32 if regime=='smooth' else 64;horizon=.4 if regime=='smooth' else .1
    source=filename(seed,regime,'coincident',n);a,m=read(source)
    indexes=list(range(len(m['branches']))) if regime=='smooth' else [j for j,b in enumerate(m['branches']) if b['arm'] in ['retained','scramble0'] and b['method'] in ['unprobed','fixed_rms','equal_work_equal_rms']]
    labels=[m['branches'][i] for i in indexes];coeff=a['start_fields'][indexes]
    positions=np.tile(a['start_positions'],(len(indexes),1,1))
    folder=ROOT/'independent-fast-data';folder.mkdir(exist_ok=True);stem=f'{regime}-s{seed}-n{n}-T{horizon:g}'+('-tight' if tight else '')
    raw=folder/(stem+'.npz');mp=raw.with_suffix('.json')
    hashes={f:digest(ROOT/f) for f in ['independent_solver.py','independent_check_fast.py']}
    if mp.exists():
        old=json.loads(mp.read_text())
        if old['sources']!=hashes or old['source_raw_sha256']!=m['raw_sha256'] or digest(raw)!=old['raw_sha256']:raise RuntimeError('Independent check cache mismatch')
        print(stem+' verified cache',flush=True);return old
    tic=time.time();config=P['regimes'][regime];solver=IndependentSolver(a['modes'],config['nu'],config['injection'])
    tol=(2e-12,2e-14) if tight else (2e-10,2e-12)
    pieces=[]
    # Bound CPU memory and keep the full algorithm/precision for every branch.
    for offset in range(0,len(indexes),8):
        r=solver.integrate(coeff[offset:offset+8],positions[offset:offset+8],horizon,*tol,max_step=.02)
        pieces.append(r);print(f'independent {stem} branches {min(offset+8,len(indexes))}/{len(indexes)}',flush=True)
    end=np.concatenate([r['end_fields'] for r in pieces]);x=np.concatenate([r['positions'] for r in pieces]);f=np.concatenate([r['tangent'] for r in pieces])
    ix=np.flatnonzero(np.isclose(a['time'],horizon))[0];reference_f=a['tangent'][ix,indexes]
    if horizon==P['horizon']:reference_end=a['end_fields'][indexes]
    else:
        # Reference restart at the same complete checkpoint; no different preparation.
        reference=DrivenSolver(n,config['nu'],config['injection'],'gpu');xp=reference.xp
        hh=reference.unpack(coeff);xx=xp.asarray(positions);ff=xp.asarray(np.tile(np.eye(3),(len(indexes),P['centers'],1,1)));tt=0.
        while tt<horizon-1e-12:
            dt=min(config['dt'],horizon-tt);hh,xx,ff,_,_=reference.step(hh,xx,ff,dt);tt+=dt
        reference_end=reference.pack(hh)
        if abs(reference.cpu(ff)-reference_f).max()>1e-11:raise RuntimeError('Reference short-prefix restart inconsistent')
        xp.get_default_memory_pool().free_all_blocks()
    fi=np.log(np.linalg.svd(f,compute_uv=False)[...,0])/horizon;fr=np.log(np.linalg.svd(reference_f,compute_uv=False)[...,0])/horizon
    weights=np.where(a['modes'][:,2]==0,1.,2.)
    norm=lambda v:np.sqrt(np.sum(abs(v)**2*weights,axis=(-2,-1)))
    response_errors=[];field_errors=[]
    for method in ['fixed_rms','equal_work_equal_rms']+(['equal_work_fixed_shape'] if regime=='smooth' else []):
        for arm in P['arms'] if regime=='smooth' else ['retained','scramble0']:
            for hist in range(2):
                kicked=next(i for i,b in enumerate(labels) if b['arm']==arm and b['history']==hist and b['method']==method)
                base=next(i for i,b in enumerate(labels) if b['arm']==arm and b['history']==hist and b['method']=='unprobed')
                response_errors.extend((fi[kicked]-fi[base])-(fr[kicked]-fr[base]))
                ri=end[kicked]-end[base];rr=reference_end[kicked]-reference_end[base]
                field_errors.append(float(norm(ri-rr)/max(norm(rr),1e-30)))
    comparison=dict(response_ftle_rmse=float(np.sqrt(np.mean(np.array(response_errors)**2))),response_field_relative=float(max(field_errors)),
                    absolute_ftle_rmse=float(np.sqrt(np.mean((fi-fr)**2))),position_error=float(abs(x-a['positions'][ix,indexes]).max()),
                    divergence=max(r['divergence'] for r in pieces),energy_budget=max(r['energy_budget'] for r in pieces),volume=max(r['volume'] for r in pieces))
    np.savez_compressed(raw,modes=a['modes'],start_fields=coeff,end_fields=end,reference_end_fields=reference_end,positions=x,tangent=f,reference_tangent=reference_f)
    result=dict(status='complete',seed=seed,regime=regime,n=n,horizon=horizon,tight=tight,pilot=pilot,sources=hashes,source_raw_sha256=m['raw_sha256'],raw_sha256=digest(raw),
                branches=labels,comparison=comparison,nfev=sum(r['nfev'] for r in pieces),accepted_steps=[r['accepted_steps'] for r in pieces],rtol=tol[0],atol=tol[1],seconds=time.time()-tic,
                solver='CPU full-complex Fourier/padded conservative flux/DOP853; independent implementation by same author',
                passed=comparison['response_ftle_rmse']<P['gates']['independent_response_ftle_rmse'] and comparison['response_field_relative']<P['gates']['independent_response_field_relative'] and comparison['energy_budget']<P['gates']['energy_budget_relative'])
    dump(mp,result);print(stem+' '+json.dumps(comparison),flush=True);return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);ap.add_argument('--regime',default='smooth');ap.add_argument('--tight',action='store_true');ap.add_argument('--pilot',action='store_true');a=ap.parse_args()
    with threadpool_limits(1):check(a.seed,a.regime,a.tight,a.pilot)
