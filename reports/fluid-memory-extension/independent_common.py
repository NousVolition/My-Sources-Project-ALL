import json,time,numpy as np
from independent_solver import IndependentSolver
from experiment import ROOT,P,digest,dump
from analysis import read,effective_regimes

def check(seed,regime,tight=False):
    n=32 if regime=='smooth' else max(64,min(effective_regimes()[regime]['grids']));dt=effective_regimes()[regime]['dt'];horizon=P['horizon']
    source=ROOT/'common-probe-data'/f'{regime}-coincident-s{seed}-n{n}-dt{dt:g}.npz';a,m=read(source)
    selected=[j for j,b in enumerate(m['branches']) if b['arm'] in ['retained','scramble0']]
    coefficients=np.stack([field for j in selected for field in [a['base_start_fields'][j],a['start_fields'][j]]])
    positions=np.stack([a['positions'][0,j] for j in selected for _ in range(2)])
    reference_fields=np.stack([field for j in selected for field in [a['base_end_fields'][j],a['end_fields'][j]]])
    reference_f=np.stack([field for j in selected for field in [a['base_tangent'][-1,j],a['tangent'][-1,j]]])
    folder=ROOT/'independent-common-data';folder.mkdir(exist_ok=True);stem=f'{regime}-s{seed}-n{n}'+('-tight' if tight else '')
    raw=folder/(stem+'.npz');mp=raw.with_suffix('.json');src={f:digest(ROOT/f) for f in ['independent_solver.py','independent_common.py']}
    if mp.exists():
        old=json.loads(mp.read_text())
        if old['sources']!=src or old['parent_raw_sha256']!=m['raw_sha256'] or old['raw_sha256']!=digest(raw):raise RuntimeError('Independent same-pulse cache mismatch')
        return old
    cfg=P['regimes'][regime];s=IndependentSolver(a['modes'],cfg['nu'],cfg['injection']);tol=(2e-12,2e-14) if tight else (2e-10,2e-12);tic=time.time()
    r=s.integrate(coefficients,positions,horizon,*tol)
    fi=np.log(np.linalg.svd(r['tangent'],compute_uv=False)[...,0])/horizon;fr=np.log(np.linalg.svd(reference_f,compute_uv=False)[...,0])/horizon
    di=(fi[1::2]-fi[::2])-(fr[1::2]-fr[::2]);weights=np.where(a['modes'][:,2]==0,1.,2.)
    norm=lambda h:np.sqrt(np.sum(abs(h)**2*weights,axis=(-2,-1)))
    ri=r['end_fields'][1::2]-r['end_fields'][::2];rr=reference_fields[1::2]-reference_fields[::2]
    comp=dict(response_ftle_rmse=float(np.sqrt(np.mean(di**2))),response_field_relative=float(np.max(norm(ri-rr)/norm(rr))),energy_budget=r['energy_budget'],volume=r['volume'],divergence=r['divergence'])
    np.savez_compressed(raw,modes=a['modes'],start_fields=coefficients,end_fields=r['end_fields'],reference_end_fields=reference_fields,tangent=r['tangent'],reference_tangent=reference_f,positions=r['positions'])
    result=dict(seed=seed,regime=regime,n=n,horizon=horizon,tight=tight,sources=src,parent_raw_sha256=m['raw_sha256'],raw_sha256=digest(raw),comparison=comp,seconds=time.time()-tic,
                passed=comp['response_ftle_rmse']<P['gates']['independent_response_ftle_rmse'] and comp['response_field_relative']<P['gates']['independent_response_field_relative'])
    dump(mp,result);print('independent common '+stem+' '+json.dumps(comp),flush=True);return result
