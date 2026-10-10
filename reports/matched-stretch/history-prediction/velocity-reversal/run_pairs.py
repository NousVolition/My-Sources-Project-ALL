"""Bounded physical velocity reversal using unchanged parent NS/Heun dynamics."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from simulate import Flow,L,initial,sample,sha,save_json
from features import index_at,future_target,neighbors,wrap


def advance(f,h,x,M,dt):
    u=f.real(h); v,J=sample(u,x)
    r1,a,s1=f.rhs(h); hp=h+dt*r1; xp=x+dt*v; dM=J@M; Mp=M+dt*dM
    vp,Jp=sample(f.real(hp),xp); r2,b,s2=f.rhs(hp)
    return f.project(h+dt/2*(r1+r2)),x+dt/2*(v+vp),M+dt/2*(dM+Jp@Mp),dt/2*(a[1]+b[1]),max(s1,s2)*dt/f.dx


def directional(x,v,J,k):
    idx=neighbors(x,k); r=wrap(x[idx]-x[:,None]); d=np.linalg.norm(r,axis=-1)
    e=r/np.maximum(d[...,None],1e-12)
    pair=np.sum(r*(v[idx]-v[:,None]),axis=-1)/np.maximum(d*d,1e-24)
    local=np.einsum('nki,nij,nkj->nk',e,J,e)
    C=np.einsum('nki,nkj->nij',r,r)/k
    _,axes=np.linalg.eigh(C); axis=axes[:,:,-1]
    axis_rate=np.einsum('ni,nij,nj->n',axis,J,axis)
    return np.column_stack([pair.mean(1),local.mean(1),axis_rate])


def identity(n): return np.broadcast_to(np.eye(3),(n,3,3)).copy()


def archive_path(seed,n,dt):
    name=f's{seed}_n{n}_nu0.02_dt{dt}.npz'
    candidates=[HERE.parent/'recorded-data'/name,HERE.parent/'data'/name,
                HERE.parent/'resolution-followup'/'recorded-data'/name]
    return next((x for x in candidates if x.exists()),None)


def run(seed,n,dt,p):
    out=HERE/'recorded-data'; out.mkdir(exist_ok=True)
    name=f's{seed}_n{n}_dt{dt}'
    dest=out/(name+'.npz'); metadata=out/(name+'.json'); checkpoint=out/(name+'_branch-state.npz')
    source={str(x.relative_to(HERE.parent.parent)).replace('\\','/'):sha(x) for x in [HERE/'run_pairs.py',HERE/'protocol.json',HERE.parent/'simulate.py',HERE.parent/'features.py',HERE.parent.parent/'numerics.py']}
    if dest.exists() and metadata.exists():
        meta=json.loads(metadata.read_text())
        if meta['sources']!=source or meta['sha256']!=sha(dest) or meta['checkpoint_sha256']!=sha(checkpoint):
            raise ValueError('Changed cached run '+name)
        return meta
    archive=archive_path(seed,n,dt)
    if archive is None: raise FileNotFoundError(name)
    am=json.loads(archive.with_suffix('.json').read_text())
    if sha(archive)!=am['sha256']: raise ValueError('Archive hash mismatch')
    with np.load(archive,allow_pickle=False) as z: old=dict(z)
    anchor=index_at(old['time'],p['branch_time'])
    f=Flow(n,p['nu'],workers=1); start=time.perf_counter()
    h=initial(f,seed); x=np.random.default_rng(seed+80000).uniform(-L/2,L/2,(p['markers'],3)); M=identity(len(x))
    if checkpoint.exists():
        # A checkpoint without verified metadata is preserved rather than silently reused.
        raise ValueError('Unverified existing branch checkpoint: '+str(checkpoint))
    for _ in range(round(p['branch_time']/dt)):
        h,x,M,_,cfl=advance(f,h,x,M,dt)
        if cfl>=.5 or not np.isfinite(h).all() or not np.isfinite(M).all(): raise FloatingPointError('Replay failure')
    v,J=sample(f.real(h),x)
    match={key:float(np.max(abs(a-old[key][anchor]))) for key,a in [('positions',x),('velocity',v),('gradient',J),('tangent',M)]}
    if max(match.values())>=1e-10: raise ValueError('Replay mismatch '+str(match))
    np.savez_compressed(checkpoint,field=h,positions=x,prior_tangent=M)
    pair={'initial_directional_rates':directional(x,v,J,p['neighbors'])}
    branch_diagnostics={}; comparisons={}
    for sign,label in [(1,'normal'),(-1,'reversed')]:
        bh=sign*h.copy(); bx=x.copy(); bM=identity(len(x)); energy0=f.inner(bh,bh)/2
        arrays={key:[] for key in ['time','positions','velocity','gradient','tangent']}; diags=[]; diss=0.; max_cfl=0.
        try:
            for step in range(round(p['duration']/dt)+1):
                bv,bJ=sample(f.real(bh),bx)
                if step%round(p['save_dt']/dt)==0:
                    for key,val in zip(arrays,[step*dt,bx.copy(),bv,bJ,bM.copy()]): arrays[key].append(val)
                    energy=f.inner(bh,bh)/2; div=f.real(1j*sum(k*c for k,c in zip(f.k,bh)))
                    diags.append({'elapsed':step*dt,'energy':energy,'energy_budget_error':(energy-energy0+diss)/energy0,
                                  'spectral_divergence_max':float(abs(div).max()),'interpolated_trace_rms':float(np.sqrt(np.mean(np.trace(bJ,axis1=1,axis2=2)**2))),
                                  'tangent_det_max_error':float(abs(np.linalg.det(bM)-1).max())})
                if step==round(p['duration']/dt): break
                bh,bx,bM,loss,cfl=advance(f,bh,bx,bM,dt); diss+=loss; max_cfl=max(max_cfl,cfl)
                if cfl>=.5 or not np.isfinite(bh).all() or not np.isfinite(bM).all(): raise FloatingPointError('Branch failure')
            if max(abs(q['energy_budget_error']) for q in diags)>=.005 or max(q['spectral_divergence_max'] for q in diags)>=1e-10:
                raise FloatingPointError('Budget/divergence failure')
        except Exception as exc:
            np.savez_compressed(out/(name+'_'+label+'_FAILED.npz'),field=bh,positions=bx,tangent=bM,step=step)
            save_json(out/(name+'_'+label+'_FAILED.json'),{'error':repr(exc),'sources':source,'diagnostics':diags})
            raise
        arrays={k:np.asarray(v) for k,v in arrays.items()}
        pair.update({label+'_'+k:val for k,val in arrays.items()})
        branch_diagnostics[label]={'max_cfl':max_cfl,'rows':diags}
    for key in ['positions','velocity','gradient']:
        inds=[index_at(old['time'],p['branch_time']+t) for t in pair['normal_time']]
        comparisons['normal_'+key+'_archive_max_error']=float(np.max(abs(pair['normal_'+key]-old[key][inds])))
    comparisons['initial_positions_equal']=bool(np.array_equal(pair['normal_positions'][0],pair['reversed_positions'][0]))
    comparisons['initial_velocity_negation_max_error']=float(abs(pair['normal_velocity'][0]+pair['reversed_velocity'][0]).max())
    comparisons['initial_gradient_negation_max_error']=float(abs(pair['normal_gradient'][0]+pair['reversed_gradient'][0]).max())
    comparisons['directional_rate_negation_max_error']=float(abs(pair['initial_directional_rates']+directional(x,-v,-J,p['neighbors'])).max())
    for horizon in p['horizons']:
        b=index_at(pair['normal_time'],horizon); e=index_at(old['time'],p['branch_time']+horizon)
        target=future_target(old['tangent'][anchor],old['tangent'][e],horizon)
        actual=future_target(pair['normal_tangent'][0],pair['normal_tangent'][b],horizon)
        comparisons[f'normal_target_archive_error_{horizon}']=float(abs(target-actual).max())
    if not comparisons['initial_positions_equal'] or max(v for v in comparisons.values() if not isinstance(v,bool))>=1e-9:
        raise ValueError('Branch control mismatch '+str(comparisons))
    # Reference positions for an actual retracing diagnostic, not a target for fitting.
    pair['recorded_past_positions']=old['positions'][[index_at(old['time'],p['branch_time']-t) for t in pair['normal_time']]]
    np.savez_compressed(dest,**pair)
    meta={'name':name,'seed':seed,'n':n,'dt':dt,'sources':source,'sha256':sha(dest),'checkpoint_sha256':sha(checkpoint),
          'archive_name':archive.name,'archive_sha256':sha(archive),'archive_metadata_sha256':sha(archive.with_suffix('.json')),
          'replay_match':match,'checks':comparisons,'diagnostics':branch_diagnostics,'seconds':time.perf_counter()-start,'complete':True}
    save_json(metadata,meta); print(name+' both branches verified',flush=True)
    return meta


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--pilot',action='store_true'); args=ap.parse_args()
    p=json.loads((HERE/'protocol.json').read_text())
    jobs=[(s,n,p['dt']) for n in p['grids'] for s in p['seeds']]
    c=p['half_step']; jobs.append((c['seed'],c['n'],c['dt']))
    if args.pilot: jobs=jobs[:1]
    records=[]
    for s,n,dt in jobs:
        m=run(s,n,dt,p); records.append({k:m[k] for k in ['name','sha256','checkpoint_sha256','seconds','complete']})
        save_json(HERE/('pilot.json' if args.pilot else 'execution.json'),{'jobs':records,'complete':len(records)==len(jobs),'protocol_sha256':sha(HERE/'protocol.json')})


if __name__=='__main__':main()
