"""Joint marker/state pilot. Reads immutable saved fluid; never evolves fluid."""
from pathlib import Path
import os, sys, json, hashlib, time
import numpy as np
from scipy import fft, linalg
ROOT=Path(__file__).resolve().parent
FLUID=ROOT.parent/'matched-stretch-study'
sys.path.insert(0,str(FLUID/'marker-stress'))
from run import interp, wrap

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def save(p,obj):
    temp=p.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8');os.replace(temp,p)

def koch(depth):
    points=np.exp(2j*np.pi*np.arange(3)/3)
    for _ in range(depth):
        out=[]
        for a,b in zip(points,np.roll(points,-1)):
            d=(b-a)/3
            out.extend([a,a+d,a+d+d*np.exp(-1j*np.pi/3),a+2*d])
        points=np.array(out)
    return np.c_[points.real,points.imag,np.zeros(len(points))]

def menger(depth):
    offsets=np.array([v for v in np.ndindex(3,3,3) if sum(x==1 for x in v)<=1]);cells=np.zeros((1,3),int)
    for _ in range(depth):cells=(3*cells[:,None]+offsets[None]).reshape(-1,3)
    return (cells+.5)/3**depth-.5

def normalize(x):
    x=x-x.mean(0);return .15*x/np.sqrt(np.mean(np.sum(x*x,axis=1)))

def geometries():
    rng=np.random.default_rng(617);m=menger(2);k=koch(3)
    ms=m[rng.choice(len(m),192,replace=False)]
    mix=np.vstack([normalize(ms[:96]),normalize(k[::2])])
    angle=2*np.pi*np.arange(192)/192;circle=np.c_[np.cos(angle),np.sin(angle),np.zeros(192)]
    ball=rng.normal(size=(192,3));ball/=np.linalg.norm(ball,axis=1)[:,None];ball*=rng.random(192)[:,None]**(1/3)
    return {n:normalize(x) for n,x in [('menger',ms),('koch',k),('mixed',mix),('circle',circle),('ball',ball)]}

def graph(x):
    d=np.linalg.norm(wrap(x[:,None]-x[None,:]),axis=2);np.fill_diagonal(d,np.inf)
    cutoff=np.sort(d,axis=1)[:,7]
    a=d<=cutoff[:,None]+1e-12;a=a|a.T
    return np.array(np.where(np.triu(a,1))).T

def make_cases():
    geo=geometries();rng=np.random.default_rng(719);noise=rng.uniform(-.1,.1,192);noise-=noise.mean()
    initials={'balanced':noise,'positive':.2+noise,'negative':-(.2+noise)}
    types=np.r_[np.full(96,20.),np.full(96,-20.)];np.random.default_rng(701).shuffle(types)
    cases=[];positions=[];states=[];aa=[];bb=[];mu=[];edges=[];base_degrees=[]
    for name,x in geo.items():
        ed=graph(x);r=wrap(x[ed[:,1]]-x[ed[:,0]]);w=np.exp(-np.sum(r*r,axis=1)/.02);degree=2*w.sum()/192
        for profile in ['fading','mutual_amplification','opposing_responses']:
            a=0 if profile=='opposing_responses' else -20
            b=types if profile=='opposing_responses' else np.full(192,10 if profile=='fading' else 30)
            for state,z in initials.items():
                for mobility in [0,20]:
                    i=len(cases);cases.append({'id':f'{name}_{profile}_{state}_m{mobility}','geometry':name,'profile':profile,'initial_state':state,'mu':mobility,'edges':len(ed)})
                    positions.append(x);states.append(z);aa.extend([a]*192);bb.extend(b);mu.extend([mobility/degree]*192);edges.append(ed+i*192);base_degrees.append(degree)
    return cases,np.array(positions).reshape(-1,3),np.array(states).reshape(-1),np.array(aa),np.array(bb),np.array(mu),np.vstack(edges),geo

def relationship_rhs(x,z,a,b,mobility,edges):
    i,j=edges.T;r=wrap(x[j]-x[i]);w=np.exp(-np.minimum(np.sum(r*r,axis=1)/.02,50))
    degree=np.bincount(i,weights=w,minlength=len(z))+np.bincount(j,weights=w,minlength=len(z))
    neighbor=np.bincount(i,weights=w*z[j],minlength=len(z))+np.bincount(j,weights=w*z[i],minlength=len(z))
    dz=a*z+b*neighbor/degree
    pair=(w*np.tanh((z[i]+z[j])/2))[:,None]*r
    drift=np.stack([np.bincount(i,weights=pair[:,d],minlength=len(z))-np.bincount(j,weights=pair[:,d],minlength=len(z)) for d in range(3)],axis=1)
    return mobility[:,None]*drift,dz

def analytic_checks():
    counts={'menger':[len(menger(d)) for d in [1,2,3]],'koch':[len(koch(d)) for d in [0,1,2,3]]}
    assert counts=={'menger':[20,400,8000],'koch':[3,12,48,192]}
    errors={}
    for name,A in [('fading',np.array([[-20.,10.],[10.,-20.]])),('mutual',np.array([[-20.,30.],[30.,-20.]])),('opposed',np.array([[0.,20.],[-20.,0.]]))]:
        z=np.array([.3,-.1]);initial=z.copy();dt=.0003125
        for _ in range(320):
            f=A@z;g=A@(z+dt*f/2);h=A@(z+dt*g/2);k=A@(z+dt*h);z+=dt*(f+2*g+2*h+k)/6
        errors[name]=float(np.max(abs(z-linalg.expm(.1*A)@initial)));assert errors[name]<1e-9
    x=geometries()['mixed'];z=np.random.default_rng(41).normal(size=192);e=graph(x);a=np.full(192,-20.);b=np.full(192,30.);mu=np.ones(192)
    drift,dz=relationship_rhs(x,z,a,b,mu,e);balance=float(abs(drift.sum(0)).max());assert balance<1e-12
    perm=np.random.default_rng(42).permutation(192);inverse=np.argsort(perm);new_e=inverse[e]
    new_drift,new_dz=relationship_rhs(x[perm],z[perm],a[perm],b[perm],mu[perm],new_e)
    pe=float(max(abs(new_drift-drift[perm]).max(),abs(new_dz-dz[perm]).max()));assert pe<1e-12
    zero,zero_z=relationship_rhs(x,np.zeros(192),a,b,mu,e);assert not np.any(zero) and not np.any(zero_z)
    u=np.broadcast_to(np.array([1.,-2.,.5])[:,None,None,None],(3,8,8,8));assert np.max(abs(interp(u,x)-[1,-2,.5]))<1e-12
    return {'construction_counts':counts,'pair_exact_solution_errors':errors,'pair_drift_sum_error':balance,'permutation_error':pe,'zero_state_zero_extra_motion':True,'constant_velocity_interpolation':True}

def main():
    if (ROOT/'state.json').exists():raise RuntimeError('Existing state must be inspected before any restart')
    fd=os.open(ROOT/'run.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    cases,x0,z0,a,b,mobility,edges,geo=make_cases();checks=analytic_checks();protocol=json.loads((ROOT/'protocol.json').read_text())
    source_hashes={str(p):sha(p) for p in [Path(__file__),ROOT/'protocol.json',FLUID/'numerics.py',FLUID/'run_suite.py',FLUID/'marker-stress/run.py']}
    np.savez_compressed(ROOT/'construction.npz',**geo,edges=edges,a=a,b=b,mobility=mobility,initial_states=z0.reshape(90,192))
    output={'protocol':protocol,'cases':cases,'checks':checks,'source_hashes':source_hashes,'fields':{},'runs':{},'new_fluid_simulations':0}
    configs=[('base64','baseline-n64-base',1,32),('step64','baseline-n64-base',1,64),('grid128','baseline-n128-base',1,64),('snapshot_coarse','baseline-n64-base',2,128)]
    try:
        for name,rid,stride,steps in configs:
            save(ROOT/'state.json',{'status':'running','pid':os.getpid(),'configuration':name,'completed':list(output['runs'])})
            record=json.loads((FLUID/'runs'/rid/'result.json').read_text());n=record['job']['n'];assert record['status']=='complete';cache={};x=x0.copy();z=z0.copy();times=[0.];xs=[x.reshape(90,192,3).copy()];zs=[z.reshape(90,192).copy()];tic=time.time();max_balance=0
            def field(index):
                if index not in cache:
                    path=FLUID/'runs'/rid/f'field-{index:03d}.npy';digest=sha(path);key=rid+'/'+path.name
                    if key in output['fields']:assert output['fields'][key]==digest
                    output['fields'][key]=digest;h=np.load(path,allow_pickle=False);assert np.isfinite(h).all();cache[index]=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
                return cache[index]
            for start in range(0,10,stride):
                u0=field(start);u1=field(start+stride);duration=.01*stride;dt=duration/steps
                def rhs(x,z,t):
                    blend=np.clip(t/duration,0,1);drift,dz=relationship_rhs(x,z,a,b,mobility,edges)
                    return (1-blend)*interp(u0,x)+blend*interp(u1,x)+drift,dz
                for k in range(steps):
                    local_t=k*dt;kx,kz=rhs(x,z,local_t);lx,lz=rhs(x+dt*kx/2,z+dt*kz/2,local_t+dt/2);mx,mz=rhs(x+dt*lx/2,z+dt*lz/2,local_t+dt/2);nx,nz=rhs(x+dt*mx,z+dt*mz,local_t+dt)
                    x+=dt*(kx+2*lx+2*mx+nx)/6;z+=dt*(kz+2*lz+2*mz+nz)/6
                assert np.isfinite(x).all() and np.isfinite(z).all()
                drift,_=relationship_rhs(x,z,a,b,mobility,edges);max_balance=max(max_balance,float(abs(drift.reshape(90,192,3).sum(1)).max()));assert max_balance<1e-10
                times.append(round(.01*(start+stride),8));xs.append(x.reshape(90,192,3).copy());zs.append(z.reshape(90,192).copy());cache.pop(start,None)
                save(ROOT/'state.json',{'status':'running','pid':os.getpid(),'configuration':name,'time':times[-1],'completed':list(output['runs'])});print(name,times[-1],flush=True)
            path=ROOT/(name+'.npz');np.savez_compressed(path,time=times,positions=xs,states=zs)
            output['runs'][name]={'source':rid,'n':n,'snapshot_dt':duration,'dt':dt,'sha256':sha(path),'seconds':time.time()-tic,'max_sum_added_drifts':max_balance,'cases':90}
            save(ROOT/'results.json',output)
        for path,digest in source_hashes.items():assert sha(Path(path))==digest
        save(ROOT/'state.json',{'status':'complete','configurations':4,'case_trajectories':360,'new_fluid_simulations':0})
    except Exception as e:
        save(ROOT/'failure.json',{'error':repr(e),'pid':os.getpid(),'completed':list(output['runs'])});raise
    finally:
        (ROOT/'run.lock').unlink(missing_ok=True)

if __name__=='__main__':main()
