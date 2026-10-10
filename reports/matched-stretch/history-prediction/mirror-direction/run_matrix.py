"""Mirror passive geometry in both physical velocity directions; reuse checkpoints."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json
import argparse
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'velocity-reversal'))
from run_pairs import advance,identity,Flow,sample,sha,save_json,neighbors,wrap
R=np.diag([-1.,1.,1.])


def mirror_field(f,h):
    field=f.real(h)[:,(-np.arange(f.n))%f.n,:,:].copy()
    field[0]*=-1
    return f.hat(field)


def cloud_target(r0,r1,duration):
    gram=np.einsum('nki,nkj->nij',r0,r0)
    cross=np.einsum('nki,nkj->nij',r0,r1)
    reg=1e-10*np.maximum(np.trace(gram,axis1=1,axis2=2),1e-12)[:,None,None]
    B=np.linalg.solve(gram+reg*np.eye(3),cross)
    return np.log(np.maximum(np.linalg.svd(B,compute_uv=False)[:,0],1e-12))/duration


def prepare(x,k):
    idx=neighbors(x,k); offsets=wrap(x[idx]-x[:,None])
    local=x[:,None]+offsets@R
    return idx,offsets,local


def sources():
    root=HERE.parent.parent
    paths=[HERE/'run_matrix.py',HERE/'protocol.json',HERE.parent/'velocity-reversal'/'run_pairs.py',HERE.parent/'simulate.py',HERE.parent/'features.py',root/'numerics.py']
    return {str(p.relative_to(root)).replace('\\','/'):sha(p) for p in paths}


def source_matches(path,digest):
    import hashlib
    b=Path(path).read_bytes(); lf=b.replace(b'\r\n',b'\n')
    return digest in {hashlib.sha256(s).hexdigest() for s in [b,lf,lf.replace(b'\n',b'\r\n')]}


def run(seed,n,dt,p):
    name=f's{seed}_n{n}_dt{dt}'; out=HERE/'recorded-data'; out.mkdir(exist_ok=True)
    dest=out/(name+'.npz'); meta_path=dest.with_suffix('.json'); src=sources()
    if dest.exists() and meta_path.exists():
        old=json.loads(meta_path.read_text())
        if old['sources']!=src or old['sha256']!=sha(dest): raise ValueError('Changed cached result '+name)
        return old
    previous=HERE.parent/'velocity-reversal'/'recorded-data'
    archived=previous/(name+'.npz'); checkpoint=previous/(name+'_branch-state.npz')
    old_meta=json.loads(archived.with_suffix('.json').read_text())
    if sha(archived)!=old_meta['sha256'] or sha(checkpoint)!=old_meta['checkpoint_sha256']: raise ValueError('Changed checkpoint/archive')
    for rel,digest in old_meta['sources'].items():
        if not source_matches(HERE.parent.parent/rel,digest): raise ValueError('Original source changed')
    with np.load(checkpoint,allow_pickle=False) as z: h0=z['field'].copy(); x0=z['positions'].copy()
    with np.load(archived,allow_pickle=False) as z: previous_arrays=dict(z)
    f=Flow(n,p['nu'],workers=1); idx,r0,local=prepare(x0,p['neighbors']); m=len(x0)
    arrays={'neighbor_indices':idx,'initial_offsets':r0,'initial_local_mirror_offsets':r0@R}
    checks={'mirror_distance_error':float(abs(np.linalg.norm(r0,axis=-1)-np.linalg.norm(r0@R,axis=-1)).max()),
            'mirror_angle_gram_error':float(abs(np.einsum('nki,nli->nkl',r0,r0)-np.einsum('nki,nli->nkl',r0@R,r0@R)).max())}
    diagnostics={}; start=time.perf_counter()
    for sign,label in [(1,'normal'),(-1,'reversed')]:
        h=sign*h0.copy(); x=np.concatenate([x0,x0@R,local.reshape(-1,3)]); M=identity(len(x))
        values={k:[] for k in ['time','global_positions','global_tangent','local_neighbor_positions']}
        diags=[]; cflmax=0.; diss=0.; energy0=f.inner(h,h)/2
        try:
            for step in range(round(p['duration']/dt)+1):
                if step%round(p['save_dt']/dt)==0:
                    save=step//round(p['save_dt']/dt)
                    error=float(abs(x[:m]-previous_arrays[label+'_positions'][save]).max())
                    terr=float(abs(M[:m]-previous_arrays[label+'_tangent'][save]).max())
                    checks[label+'_original_position_error']=max(checks.get(label+'_original_position_error',0),error)
                    checks[label+'_original_tangent_error']=max(checks.get(label+'_original_tangent_error',0),terr)
                    if max(error,terr)>=1e-9: raise ValueError('Original trajectory control failed')
                    for key,val in zip(values,[step*dt,x[m:2*m].copy(),M[m:2*m].copy(),x[2*m:].reshape(m,p['neighbors'],3).copy()]): values[key].append(val)
                    energy=f.inner(h,h)/2; div=f.real(1j*sum(k*c for k,c in zip(f.k,h)))
                    diags.append({'elapsed':step*dt,'energy':energy,'energy_budget_error':(energy-energy0+diss)/energy0,
                                  'spectral_divergence_max':float(abs(div).max()),'probe_tangent_volume_error_max':float(abs(np.linalg.det(M)-1).max())})
                if step==round(p['duration']/dt): break
                h,x,M,loss,cfl=advance(f,h,x,M,dt); diss+=loss; cflmax=max(cflmax,cfl)
                if cfl>=.5 or not np.isfinite(h).all() or not np.isfinite(M).all(): raise FloatingPointError('Numerical branch failure')
            if max(abs(d['energy_budget_error']) for d in diags)>=.005 or max(d['spectral_divergence_max'] for d in diags)>=1e-10: raise FloatingPointError('Budget/divergence failure')
        except Exception as exc:
            np.savez_compressed(out/(name+'_'+label+'_FAILED.npz'),field=h,positions=x,tangent=M,step=step)
            save_json(out/(name+'_'+label+'_FAILED.json'),{'error':repr(exc),'sources':src,'diagnostics':diags})
            raise
        arrays.update({label+'_'+key:np.asarray(val) for key,val in values.items()})
        diagnostics[label]={'max_cfl':cflmax,'rows':diags}
    # One complete-state reflection control is enough to distinguish symmetry from intervention.
    if seed==p['seeds'][0] and n==p['grids'][0] and dt==p['dt']:
        for sign,label in [(1,'normal'),(-1,'reversed')]:
            h=sign*mirror_field(f,h0); x=x0@R; M=identity(m)
            values=[]; tangents=[]
            for step in range(round(p['duration']/dt)+1):
                if step%round(p['save_dt']/dt)==0:
                    save=step//round(p['save_dt']/dt)
                    pe=float(abs(x-previous_arrays[label+'_positions'][save]@R).max())
                    te=float(abs(M-R@previous_arrays[label+'_tangent'][save]@R).max())
                    checks[label+'_full_mirror_position_error']=max(checks.get(label+'_full_mirror_position_error',0),pe)
                    checks[label+'_full_mirror_tangent_error']=max(checks.get(label+'_full_mirror_tangent_error',0),te)
                    values.append(x.copy()); tangents.append(M.copy())
                    if max(pe,te)>=1e-9: raise ValueError('Whole-state symmetry control failed')
                if step==round(p['duration']/dt): break
                h,x,M,_,cfl=advance(f,h,x,M,dt)
                if cfl>=.5 or not np.isfinite(h).all() or not np.isfinite(M).all(): raise FloatingPointError('Full mirror failed')
            arrays[label+'_full_mirror_positions']=np.asarray(values); arrays[label+'_full_mirror_tangent']=np.asarray(tangents)
    np.savez_compressed(dest,**arrays)
    meta={'name':name,'seed':seed,'n':n,'dt':dt,'sources':src,'sha256':sha(dest),
          'reused_archive_sha256':sha(archived),'reused_checkpoint_sha256':sha(checkpoint),'checks':checks,
          'diagnostics':diagnostics,'seconds':time.perf_counter()-start,'complete':True}
    save_json(meta_path,meta); print(name+' mirrored probes verified in both directions',flush=True)
    return meta


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--pilot',action='store_true'); args=ap.parse_args()
    p=json.loads((HERE/'protocol.json').read_text())
    jobs=[(s,n,p['dt']) for n in p['grids'] for s in p['seeds']]
    c=p['half_step']; jobs.append((c['seed'],c['n'],c['dt']))
    if args.pilot: jobs=jobs[:1]
    records=[]
    for s,n,dt in jobs:
        meta=run(s,n,dt,p); records.append({k:meta[k] for k in ['name','sha256','seconds','complete']})
        save_json(HERE/('pilot.json' if args.pilot else 'execution.json'),{'jobs':records,'complete':len(records)==len(jobs),'sources':sources()})


if __name__=='__main__':main()
