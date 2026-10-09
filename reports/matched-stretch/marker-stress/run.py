"""Stress passive labels using immutable stored NS velocities; never evolve fluid."""
from pathlib import Path
import os,sys,json,hashlib,time
import numpy as np
from scipy import fft,ndimage
ROOT=Path(__file__).resolve().parent;STUDY=ROOT.parent;L=6.
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def wrap(x):return (x+3)%6-3
def interp(u,p):
    n=u.shape[1];q=((p+3)%6)*n/6;a=np.floor(q).astype(int);b=q-a;out=np.zeros((len(p),3))
    for i in (0,1):
        for j in (0,1):
            for k in (0,1):
                s=np.array([i,j,k]);idx=(a+s)%n;w=np.prod(np.where(s,b,1-b),axis=1)
                out+=u[:,idx[:,0],idx[:,1],idx[:,2]].T*w[:,None]
    return out
def rk4(x,dt,v,t):
    a=v(x,t);b=v(x+dt*a/2,t+dt/2);c=v(x+dt*b/2,t+dt/2);d=v(x+dt*c,t+dt)
    return x+dt*(a+2*b+2*c+d)/6
def known_checks():
    initial=np.array([[1.,0.,0.]]);errors=[]
    for dt in (.1,.05):
        p=initial.copy()
        for k in range(round(1/dt)):p=rk4(p,dt,lambda x,t:np.stack([-x[:,1],x[:,0],0*x[:,2]],axis=1),k*dt)
        errors.append(float(np.linalg.norm(p-[np.cos(1),np.sin(1),0])))
    assert errors[0]/errors[1]>15
    field=np.broadcast_to(np.array([1.,-2.,.5])[:,None,None,None],(3,8,8,8));p=np.array([[2.99,0,0],[-3.01,0,0]])
    assert np.allclose(interp(field,p),[1,-2,.5]) and np.allclose(wrap(p[1]-p[0]),0)
    q=rk4(p,.1,lambda x,t:interp(field,x),0);assert np.allclose(q-p,[.1,-.2,.05])
    return {'rotation_error_ratio':errors[0]/errors[1],'uniform_translation_and_periodic_coordinates':'passed'}
def seeds(center):
    offset=np.array(np.meshgrid(*([np.arange(-2,3)*L/64]*3),indexing='ij')).reshape(3,-1).T
    base=np.vstack([offset,offset+center]);rng=np.random.default_rng(104);d=rng.normal(size=base.shape);d*=L/64*.01/np.linalg.norm(d,axis=1)[:,None]
    edges=[]
    for g in (0,1):
        for idx in np.ndindex(5,5,5):
            for axis in range(3):
                if idx[axis]<4:
                    nxt=list(idx);nxt[axis]+=1;edges.append([125*g+np.ravel_multi_index(idx,(5,5,5)),125*g+np.ravel_multi_index(tuple(nxt),(5,5,5))])
    return np.vstack([base,base+d,base[:5]]),np.array(edges)
def metrics(x,v,edges,initial_neighbors):
    groups=[]
    for g in (0,1):
        sl=slice(g*125,(g+1)*125);p=x[sl];vel=v[sl];ij=edges[g*300:(g+1)*300];r=wrap(x[ij[:,1]]-x[ij[:,0]]);dv=v[ij[:,1]]-v[ij[:,0]];distance=np.linalg.norm(r,axis=1)
        radial=np.sum(r*dv,axis=1)/np.maximum(distance,1e-30);mean=vel.mean(axis=0);c=vel-mean
        dd=np.linalg.norm(wrap(p[:,None]-p[None,:]),axis=2);np.fill_diagonal(dd,np.inf);neighbors=np.argsort(dd,axis=1,kind='stable')[:,:6]
        retention=np.mean([len(set(a)&set(b))/6 for a,b in zip(neighbors,initial_neighbors[g])]);amp=np.linalg.norm(wrap(x[250+g*125:250+(g+1)*125]-p),axis=1)/(L/64*.01)
        groups.append({'mean_velocity':mean.tolist(),'mean_speed':float(np.linalg.norm(mean)),'velocity_deviation_rms':float(np.sqrt(np.mean(np.sum(c*c,axis=1)))),'speed_quantiles':np.quantile(np.linalg.norm(vel,axis=1),[0,.1,.5,.9,1]).tolist(),'neighbor_distance_quantiles':np.quantile(distance,[0,.1,.5,.9,1]).tolist(),'approaching_fraction':float(np.mean(radial<0)),'median_radial_speed':float(np.median(radial)),'median_turn_rate':float(np.median(np.linalg.norm(dv-r*radial[:,None]/np.maximum(distance[:,None],1e-30),axis=1)/np.maximum(distance,1e-30))),'initial_neighbor_retention':float(retention),'seed_amplification_rms':float(np.sqrt(np.mean(amp*amp))),'seed_amplification_max':float(amp.max())})
    return groups
def main():
    if (ROOT/'state.json').exists():
        old=json.loads((ROOT/'state.json').read_text());raise RuntimeError('Existing stress study state: inspect before rerunning '+str(old))
    lock=ROOT/'run.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        protocol=json.loads((ROOT/'protocol.json').read_text());sourcehashes={n:sha(STUDY/n) for n in ('numerics.py','run_suite.py')};ownhash=sha(Path(__file__))
        out={'protocol':protocol,'source_hashes':sourcehashes,'run_script_sha256':ownhash,'analytic_checks':known_checks(),'fields':{},'runs':{},'new_fluid_simulations':0}
        record0=json.loads((STUDY/'runs/baseline-n64-base/result.json').read_text());seed,edges=seeds(record0['series'][0]['peak_location'])
        out.update(initial_positions=seed[:250].tolist(),perturbed_initial_positions=seed[250:500].tolist(),edges=edges.tolist(),groups=['central','periodic_join'],markers=250,companions=250,clones=5)
        neighbors=[]
        for g in (0,1):
            p=seed[g*125:(g+1)*125];dd=np.linalg.norm(wrap(p[:,None]-p[None,:]),axis=2);np.fill_diagonal(dd,np.inf);neighbors.append(np.argsort(dd,axis=1,kind='stable')[:,:6])
        configs=[('tracking_coarse','baseline-n64-base',1,16),('tracking_medium','baseline-n64-base',1,32),('reference','baseline-n64-base',1,64),('saved_time_coarse','baseline-n64-base',2,128),('fluid_half','baseline-n64-half',1,64),('grid128_medium','baseline-n128-base',1,32),('grid128','baseline-n128-base',1,64)]
        for name,rid,stride,substeps in configs:
            save(ROOT/'state.json',{'status':'running','pid':os.getpid(),'active':name,'completed':list(out['runs'])});folder=STUDY/'runs'/rid;record=json.loads((folder/'result.json').read_text());assert record['series'][-1]['t']>=.30 and record['job']['nu']==.001
            assert record['source_hashes']['numerics.py']==sourcehashes['numerics.py'];n=record['job']['n'];source_rows=[r for r in record['series'] if r['t']<=.30+1e-12];assert len(source_rows)==31
            (ROOT/'source-records').mkdir(exist_ok=True)
            record_path=ROOT/'source-records'/f'{rid}.json'
            if not record_path.exists():save(record_path,{'job':record['job'],'dt':record['dt'],'source_hashes':record['source_hashes'],'series':source_rows,'scope':'frozen available observations through0.30; not a full-run completion claim'})
            cache={}
            def field(i):
                if i not in cache:
                    path=folder/f'field-{i:03d}.npy';h=np.load(path,allow_pickle=False);assert np.isfinite(h).all();digest=sha(path);key=rid+'/'+path.name
                    if key in out['fields']:assert out['fields'][key]==digest
                    out['fields'][key]=digest;cache[i]=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
                return cache[i]
            x=seed.copy();times=[0.];pos=[x.copy()];u0=field(0);vel=[interp(u0,x)];rows=[metrics(x,vel[-1],edges,neighbors)]
            q=((seed+3)%6)*n/6;independent=np.stack([ndimage.map_coordinates(ui,q.T,order=1,mode='grid-wrap') for ui in u0],axis=1);interpolation_error=float(np.max(abs(vel[-1]-independent)));assert interpolation_error<1e-11
            for i in range(0,30,stride):
                u0=field(i);u1=field(i+stride);duration=.01*stride;dt=duration/substeps
                def v(p,t):
                    a=np.clip(t/duration,0,1);return (1-a)*interp(u0,p)+a*interp(u1,p)
                for k in range(substeps):x=rk4(x,dt,v,k*dt)
                assert np.isfinite(x).all() and np.array_equal(x[500:],x[:5])
                times.append(round((i+stride)*.01,8));pos.append(x.copy());vel.append(interp(u1,x));rows.append(metrics(x,vel[-1],edges,neighbors));cache.pop(i,None)
            np.savez_compressed(ROOT/(name+'.npz'),time=times,positions=pos,velocities=vel)
            out['runs'][name]={'source':rid,'n':n,'tracking_dt':dt,'snapshot_dt':duration,'times':times,'metrics':rows,'mean_fluid_velocity':source_rows[0]['mean_velocity'],'initial_fluid_energy':source_rows[0]['energy'],'interpolation_check_max_error':interpolation_error,'exact_clones_agree':True,'trajectory_sha256':sha(ROOT/(name+'.npz'))}
            save(ROOT/'results.json',out);print(name+' complete through 0.30',flush=True)
        assert all(sha(STUDY/n)==s for n,s in sourcehashes.items()) and sha(Path(__file__))==ownhash
        save(ROOT/'state.json',{'status':'complete','pid':os.getpid(),'completed':list(out['runs']),'saved_fields_verified':len(out['fields'])});print('Stress tracking complete.',flush=True)
    except Exception as e:
        save(ROOT/'failure.json',{'error':repr(e),'pid':os.getpid(),'at':time.time()});raise
    finally:lock.unlink(missing_ok=True)
if __name__=='__main__':main()
