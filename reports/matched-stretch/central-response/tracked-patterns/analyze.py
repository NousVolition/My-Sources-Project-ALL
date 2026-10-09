"""Measure passive pair relationships from immutable saved fields; no evolution."""
from pathlib import Path
import os,sys,json,hashlib
import numpy as np
from scipy import fft,ndimage
ROOT=Path(__file__).resolve().parent
STUDY=ROOT.parents[2]/'adversarial-vortex-study'
L=6.
def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()
def wrap(x):return (x+L/2)%L-L/2
def trilinear(u,p):
    n=u.shape[1];q=((p+L/2)%L)*n/L;a=np.floor(q).astype(int);b=q-a;out=np.zeros((len(p),3))
    for i in (0,1):
        for j in (0,1):
            for k in (0,1):
                shift=np.array([i,j,k]);ind=(a+shift)%n;w=np.prod(np.where(shift,b,1-b),axis=1)
                out+=u[:,ind[:,0],ind[:,1],ind[:,2]].T*w[:,None]
    return out
def edges_for(p):
    groups=np.r_[np.zeros(64,int),np.repeat(np.arange(1,5),32)]
    edges={}
    for g in range(5):
        ids=np.flatnonzero(groups==g)
        for i,j in zip(ids,np.roll(ids,-1)):edges[tuple(sorted((int(i),int(j))))]='central' if g==0 else 'outer'
    # One closest differently seeded label per point, with lowest ID breaking ties.
    for i in range(len(p)):
        d=np.linalg.norm(wrap(p-p[i]),axis=1);d[groups==groups[i]]=np.inf
        j=int(np.argmin(d));edges[tuple(sorted((i,j))) ]='across'
    return groups,[{'i':i,'j':j,'kind':kind} for (i,j),kind in sorted(edges.items())]
def relations(p,v,edges,speed_floor):
    i=np.array([e['i'] for e in edges]);j=np.array([e['j'] for e in edges])
    r=wrap(p[j]-p[i]);distance=np.linalg.norm(r,axis=1);dv=v[j]-v[i];common=(v[i]+v[j])/2
    assert np.min(distance)>1e-10
    radial=np.sum(r*dv,axis=1)/distance
    speedi=np.linalg.norm(v[i],axis=1);speedj=np.linalg.norm(v[j],axis=1)
    valid=(speedi>speed_floor)&(speedj>speed_floor)
    cosine=np.full(len(edges),np.nan);cosine[valid]=np.clip(np.sum(v[i][valid]*v[j][valid],axis=1)/(speedi[valid]*speedj[valid]),-1,1)
    transverse=dv-r*radial[:,None]/distance[:,None]
    return dict(distance=distance,radial_speed=radial,separation_rate=radial/distance,
        relative_speed=np.linalg.norm(dv,axis=1),shared_speed=np.linalg.norm(common,axis=1),
        direction_cosine=cosine,turn_rate=np.linalg.norm(transverse,axis=1)/distance)
def checks():
    p=np.array([[2.9,0,0],[-2.9,0,0]]);edge=[dict(i=0,j=1)]
    t=relations(p,np.array([[1,2,0],[1,2,0]]),edge,1e-12)
    assert abs(t['distance'][0]-.2)<1e-12 and t['radial_speed'][0]==0
    p=np.array([[-.2,0,0],[.2,0,0]])
    t=relations(p,np.array([[0,-.4,0],[0,.4,0]]),edge,1e-12)
    assert abs(t['turn_rate'][0]-2)<1e-12 and t['radial_speed'][0]==0 and t['direction_cosine'][0]==-1
    v=np.array([[-.6,0,0],[.6,0,0]])
    t=relations(p,v,edge,1e-12);s=relations(p,v+np.array([7,2,3]),edge,1e-12)
    assert abs(t['separation_rate'][0]-3)<1e-12 and np.allclose(t['relative_speed'],s['relative_speed'])
    return ['periodic-neighbor distance','uniform translation','solid-body rotation','linear extension','relative-velocity invariance under uniform translation']
def extract():
    hashes={k:sha(STUDY/k) for k in ('solver.py','initial_design.py','protocol.py','run_study.py')}
    jobs=[f'{case}-{method}-n{n}-{step}' for case in ('aligned','compressive') for method,n,step in [('fourier',112,'base'),('fourier',112,'half'),('fourier',160,'base'),('fd4',160,'base')]]
    data={'created_from':'existing complete unforced surroundings runs','new_simulations':0,'domain_side':L,'viscosity':.001,'external_force':0,'runs':{},'source_hashes':hashes,'verified_fields':[], 'analytic_checks':checks()}
    seed=None
    for rid in jobs:
        folder=STUDY/'runs'/rid;record=json.loads((folder/'result.json').read_text());assert record['status']=='complete' and record['source_hashes']==hashes
        job=record['job'];assert job['nu']==.001
        frames=[]
        for path in sorted(folder.glob('field-t*.npz')):
            with np.load(path,allow_pickle=False) as z:
                h=z['h'];p=z['markers'];t=float(z['t']);assert json.loads(str(z['source_hashes']))==hashes and json.loads(str(z['job_json']))==job
            assert np.isfinite(h).all() and np.isfinite(p).all() and p.shape==(192,3)
            if t==0:
                if seed is None:
                    seed=p.copy();groups,edges=edges_for(seed);data.update(groups=groups.tolist(),edges=edges,initial_positions=seed.tolist())
                else:assert np.array_equal(seed,p)
            u=fft.irfftn(h,s=(job['n'],)*3,axes=(-3,-2,-1),workers=1)
            v=trilinear(u,p)
            independent=np.stack([ndimage.map_coordinates(ui,((p+3)%6*job['n']/6).T,order=1,mode='grid-wrap') for ui in u],axis=1)
            error=float(np.max(abs(v-independent)));assert error<1e-11
            if t==0:floor=max(float(np.linalg.norm(v,axis=1).max())*1e-10,1e-12)
            rel=relations(p,v,edges,floor)
            metrics={key:[None if not np.isfinite(x) else float(x) for x in value] for key,value in rel.items()}
            summaries={}
            for kind in ('central','outer','across'):
                mask=np.array([e['kind']==kind for e in edges]);summaries[kind]={k:float(np.nanmedian(values[mask])) for k,values in rel.items()}
                summaries[kind]['valid_direction_pairs']=int(np.isfinite(rel['direction_cosine'][mask]).sum())
            frames.append(dict(t=round(t,8),positions=p.tolist(),velocity=v.tolist(),speed_floor=floor,metrics=metrics,summary=summaries))
            data['verified_fields'].append(dict(run=rid,file=path.name,sha256=sha(path),independent_interpolation_max_error=error))
            print(rid+' time '+str(round(t,2)),flush=True)
        assert [f['t'] for f in frames]==[0,.1,.2,.3,.4]
        data['runs'][rid]=dict(job=job,frames=frames,source_result_sha256=sha(folder/'result.json'))
        save(ROOT/'measurements.partial.json',data)
    save(ROOT/'measurements.json',data)
    print(json.dumps({'runs':len(data['runs']),'fields':len(data['verified_fields']),'labels':192,'pairs':len(edges),'pair_types':{k:sum(e['kind']==k for e in edges) for k in ('central','outer','across')}}))
if __name__=='__main__':extract()
