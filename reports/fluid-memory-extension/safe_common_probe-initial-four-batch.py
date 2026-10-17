import json,time
import numpy as np
from final_experiment import ROOT,P,DrivenSolver,basis,digest,dump,DIAG
from analysis import read,filename,idx

DESIGN=json.loads((ROOT/'common-probe-protocol.json').read_text())

def construct(s,states,center):
    xp=s.xp
    span=xp.concatenate([basis(s,np.array(center)+off,DESIGN['width'],DESIGN['cutoff']) for off in np.array(DESIGN['offsets'])])
    gram=s.cpu(xp.stack([s.inner(span[i],span[j]) for i in range(len(span)) for j in range(len(span))])).reshape(len(span),len(span))
    eig,v=np.linalg.eigh(gram);keep=eig>eig[-1]*1e-10
    transform=(v[:,keep]/np.sqrt(eig[keep])).T
    orth=xp.einsum('ij,jcxyz->icxyz',xp.asarray(transform),span)
    constraints=s.cpu(xp.stack([s.inner(states,q) for q in orth],axis=1))
    preferred=s.cpu(xp.stack([s.inner(span[2],q) for q in orth]))
    u,singular,vh=np.linalg.svd(constraints,full_matrices=True)
    rank=int(np.sum(singular>min(max(singular[0],1e-300)*1e-9,1e-10)));null=vh[rank:]
    projected=null.T@(null@preferred);alignment=float(np.linalg.norm(projected));fallback=False;fallback_basis_index=None
    if alignment<1e-6:
        # Select a fixed physical basis direction, never an arbitrary SVD
        # vector whose orientation can depend on the order of state labels.
        fallback=True
        for candidate in range(len(span)):
            anchor=s.cpu(xp.stack([s.inner(span[candidate],q) for q in orth]))
            proposed=null.T@(null@anchor)
            if np.linalg.norm(proposed)>=1e-6:
                projected=proposed;fallback_basis_index=candidate;break
        else:raise RuntimeError('No label-independent feasible common impulse')
    pulse=xp.einsum('i,icxyz->cxyz',xp.asarray(projected),orth)
    pulse*=DESIGN['strength']/xp.sqrt(s.inner(pulse,pulse))
    linear=s.cpu(s.inner(states,pulse));norm=float(s.cpu(xp.sqrt(s.inner(pulse,pulse))))
    increments=linear+.5*norm**2
    axes=xp.arange(s.n)*s.dx
    coordinates=xp.meshgrid(axes,axes,axes,indexing='ij',sparse=True)
    radius2=sum(((coordinates[j]-center[j]+np.pi)%(2*np.pi)-np.pi)**2 for j in range(3))
    energy=xp.sum(s.real(pulse[None])[0]**2,axis=0);reference_energy=xp.sum(s.real(span[2][None])[0]**2,axis=0)
    localization=dict(energy_fraction_radius2=float(s.cpu(xp.sum(energy*(radius2<4))/xp.sum(energy))),
                      rms_periodic_radius=float(s.cpu(xp.sqrt(xp.sum(energy*radius2)/xp.sum(energy)))),
                      preferred_energy_fraction_radius2=float(s.cpu(xp.sum(reference_energy*(radius2<4))/xp.sum(reference_energy))))
    audit=dict(span_dimension=len(span),retained_span_dimension=int(keep.sum()),constraint_rank=rank,null_dimension=len(null),preferred_projection_norm=alignment,fallback=fallback,fallback_basis_index=fallback_basis_index,
               linear_work=linear.tolist(),increments=increments.tolist(),rms=norm,work_error=float(abs(increments-DESIGN['energy_increment']).max()),rms_error=abs(norm-DESIGN['strength']),
               preferred_cosine=float(s.cpu(s.inner(pulse,span[2])))/norm,localization=localization,
               gram_eigenvalues=eig.tolist(),constraint_singular_values=singular.tolist())
    if audit['work_error']>P['gates']['work_error'] or audit['rms_error']>P['gates']['rms_error']:raise RuntimeError('Common-probe work/strength construction failed')
    return pulse,audit

def run_one(seed,regime,layout,n,dt=None):
    cfg=P['regimes'][regime];dt=cfg['dt'] if dt is None else dt;parent=filename(seed,regime,layout,n,dt)
    a,m=read(parent);folder=ROOT/'common-probe-data';folder.mkdir(exist_ok=True);stem=f'{regime}-{layout}-s{seed}-n{n}-dt{dt:g}'
    raw=folder/(stem+'.npz');mp=raw.with_suffix('.json');src={f:digest(ROOT/f) for f in ['safe_common_probe.py','common-probe-protocol.json','common-identity-amendment.json','corrected_solver.py','corrected_observation.py','final_experiment.py']}
    if mp.exists():
        old=json.loads(mp.read_text())
        if old['sources']!=src or old['parent_raw_sha256']!=m['raw_sha256'] or old['raw_sha256']!=digest(raw):raise RuntimeError('Common-probe cache mismatch')
        print(stem+' common probe verified cache',flush=True);return old
    tic=time.time();s=DrivenSolver(n,cfg['nu'],cfg['injection'],'gpu');xp=s.xp
    base_indexes=[idx(m,arm,h,'unprobed') for arm in P['arms'] for h in range(2)]
    states=s.unpack(a['start_fields'][base_indexes]);center=np.array([np.pi]*3)+P['layouts'][layout]['probe_shift'];pulse,audit=construct(s,states,center)
    start=s.pack(states+pulse);times=a['time'];pieces=[];ends=[]
    batch=4 if n>=128 else 8
    xp.get_default_memory_pool().free_all_blocks()
    for off in range(0,len(states),batch):
        h=(states[off:off+batch]+pulse).copy();b=len(h);x=xp.asarray(np.tile(a['start_positions'],(b,1,1)));f=xp.asarray(np.tile(np.eye(3),(b,P['centers'],1,1)))
        e0=.5*s.inner(h,h);budget=xp.zeros(b);cfl=xp.zeros(b);t=0.;xx=[];ff=[];dd=[]
        for target in times:
            while t<target-1e-12:
                step=min(dt,target-t);h,x,f,db,cf=s.step(h,x,f,step);budget+=db;cfl=xp.maximum(cfl,cf);t+=step
            xx.append(s.cpu(x));ff.append(s.cpu(f));dd.append(s.diagnostics(h,x,f,budget,e0,cfl))
        pieces.append((xx,ff,dd));ends.append(s.pack(h))
    positions=np.concatenate([r[0] for r in pieces],axis=1);tangent=np.concatenate([r[1] for r in pieces],axis=1);diagnostics=np.concatenate([r[2] for r in pieces],axis=1)
    np.savez_compressed(raw,modes=a['modes'],pulse=s.pack(pulse[None])[0],start_fields=start,end_fields=np.concatenate(ends),
                        base_start_fields=a['start_fields'][base_indexes],base_end_fields=a['end_fields'][base_indexes],time=times,
                        positions=positions,tangent=tangent,base_positions=a['positions'][:,base_indexes],base_tangent=a['tangent'][:,base_indexes],diagnostics=diagnostics)
    maxima={name:float(np.max(abs(diagnostics[...,j]))) if j==1 else float(np.max(diagnostics[...,j])) for j,name in enumerate(DIAG)}
    result=dict(config=dict(seed=seed,regime=regime,layout=layout,n=n,dt=dt),sources=src,parent_raw_sha256=m['raw_sha256'],raw_sha256=digest(raw),
                branches=[m['branches'][i] for i in base_indexes],audit=audit,maxima=maxima,seconds=time.time()-tic)
    dump(mp,result);print('common probe '+stem+f' {result["seconds"]:.1f}s work_error={audit["work_error"]:.2g}',flush=True);xp.get_default_memory_pool().free_all_blocks();return result
