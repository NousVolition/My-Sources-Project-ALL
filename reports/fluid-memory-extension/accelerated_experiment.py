"""Follow-up experiments. The reference evolution code is preserved byte for byte."""
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from reference_solver import random_field, interventions
from fast_observation import FastSolver as Solver

ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'protocol.json').read_text())
DIAG=['energy','energy_budget_relative','divergence_rms','tail_enstrophy','tangent_volume_error','gradient_trace','max_cfl']
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')
def sources():return {f:digest(ROOT/f) for f in ['protocol.json','reference_solver.py','fast_observation.py','accelerated_experiment.py']}

class DrivenSolver(Solver):
    def __init__(self,n,nu,injection=0,backend='gpu'):
        super().__init__(n,nu,backend);self.injection=injection
        self.forced_band=self.xp.asarray((np.max(abs(self.k_cpu),axis=0)<=2)&(np.sum(self.k_cpu**2,axis=0)>0))
    def rhs(self,h,force):
        if self.injection:
            low=h*self.forced_band
            power=self.inner(low,low)
            if float(self.cpu(self.xp.min(power)))<1e-20:raise RuntimeError('Low-mode forcing singularity')
            force=force+low*(self.injection/power)[:,None,None,None,None]
        return super().rhs(h,force)

def basis(s,center,width=.65,cutoff=4,high=False):
    xp=s.xp
    psi=xp.exp(-.5*width**2*s.k2)*xp.exp(-1j*sum(s.k[j]*center[j] for j in range(3)))*s.n**3
    band=xp.max(abs(s.k),axis=0)<=cutoff
    if high:band &= xp.max(abs(s.k),axis=0)>2
    psi*=band
    kx,ky,kz=s.k;zero=xp.zeros_like(psi)
    q=xp.stack([xp.stack([zero,1j*kz*psi,-1j*ky*psi]),xp.stack([-1j*kz*psi,zero,1j*kx*psi]),xp.stack([1j*ky*psi,-1j*kx*psi,zero])])
    q=s.project(q)
    gram=s.cpu(xp.stack([s.inner(q[i],q[j]) for i in range(3) for j in range(3)])).reshape(3,3)
    eig,v=np.linalg.eigh(gram);transform=(v/np.sqrt(eig))@v.T
    return xp.einsum('ij,jcxyz->icxyz',xp.asarray(transform),q)

def pulse_coefficients(c):
    """c=<u,q_j>; q_j are orthonormal, real, localized divergence-free fields."""
    c=np.asarray(c,float);strength=P['probe_rms'];work=.5*strength**2
    out=np.zeros((len(c),3,3));out[:,0,2]=strength
    dot=c[:,2];root=np.sqrt(dot*dot+2*work)
    amplitude=np.where(dot>=0,2*work/(root+dot),root-dot)
    out[:,1,2]=amplitude
    preferred=np.array([0.,0.,1.]);fallback=[]
    for i,v in enumerate(c):
        vv=float(v@v)
        direction=preferred.copy() if vv<1e-28 else preferred-v*(v@preferred)/vv
        use=np.linalg.norm(direction)<1e-6
        if use:
            axis=np.eye(3)[np.argmin(abs(v))];direction=axis-v*(v@axis)/vv
        out[i,2]=strength*direction/np.linalg.norm(direction);fallback.append(bool(use))
    increments=np.einsum('bmj,bj->bm',out,c)+.5*np.sum(out*out,axis=-1)
    return out,dict(increments=increments.tolist(),rms=np.linalg.norm(out,axis=-1).tolist(),linear_work=np.einsum('bmj,bj->bm',out,c).tolist(),
                   equal_work_error=float(abs(increments[:,1:]-work).max()),equal_rms_error=float(abs(np.linalg.norm(out[:,2],axis=1)-strength).max()),
                   orthogonal_linear_work_error=float(abs(np.einsum('bj,bj->b',out[:,2],c)).max()),fallback=fallback,
                   preferred_cosine=(out[:,2,2]/strength).tolist())

def cloud(seed,center):return np.array(center)+np.random.default_rng(seed+170000).uniform(-.7,.7,(P['centers'],3))
def restrict(s,coeff,modes):
    result=np.zeros((len(coeff),3,len(s.modes_cpu)),complex);lookup={tuple(k):i for i,k in enumerate(modes)}
    for j,k in enumerate(s.modes_cpu):
        if tuple(k) in lookup:result[:,:,j]=coeff[:,:,lookup[tuple(k)]]
    return s.unpack(result)

def stats(s,h,cfl=0):
    xp=s.xp;energy=.5*s.inner(h,h);w2=s.inner(s.curl(h),s.curl(h));diss=s.nu*w2
    force=xp.zeros_like(h);rhs,_,_=s.rhs(h,force)
    low=h*s.forced_band
    forcing=low*(s.injection/xp.maximum(s.inner(low,low),1e-300))[:,None,None,None,None] if s.injection else 0
    nonlinear=rhs+s.nu*s.k2*h-forcing
    flux={str(k):float(s.cpu(-s.inner(h*(xp.sqrt(s.k2)<=k),nonlinear))[0]) for k in [2,4,8,12,16,24]}
    kmax=int(np.max(abs(s.modes_cpu)))
    eps=float(s.cpu(diss)[0]);u=float(s.cpu(xp.sqrt(2*energy))[0]);om=float(s.cpu(w2)[0])
    shell=np.rint(np.sqrt(np.sum(s.modes_cpu**2,axis=1))).astype(int)
    coeff=s.pack(h)[0];weights=np.where(s.modes_cpu[:,2]==0,1.,2.)
    spectrum=np.bincount(shell,weights=.5*np.sum(abs(coeff)**2,axis=0)*weights,minlength=int(shell.max())+1)
    return dict(energy=float(s.cpu(energy)[0]),dissipation=eps,velocity_rms=u,Re_lambda=u/np.sqrt(3)*np.sqrt(5*u*u/max(om,1e-300))/s.nu,
                kmax_eta=kmax*(s.nu**3/max(eps,1e-300))**.25,
                tail_enstrophy=float(s.cpu(s.inner(s.curl(h)*s.high,s.curl(h))/xp.maximum(w2,1e-300))[0]),
                cfl=float(cfl),energy_spectrum=spectrum.tolist(),flux=flux)

def spinup(seed,regime,backend='gpu'):
    cfg=P['regimes'][regime];n=max(cfg['grids']);folder=ROOT/'fast-spinup';folder.mkdir(exist_ok=True)
    name=f'{regime}-s{seed}-n{n}-{backend}';raw=folder/(name+'.npz');mp=raw.with_suffix('.json')
    original=ROOT/'spinup'/(name+'.npz')
    if original.with_suffix('.json').exists():
        import experiment as original_version
        om=json.loads(original.with_suffix('.json').read_text())
        if om['sources']!=original_version.sources() or digest(original)!=om['raw_sha256']:raise RuntimeError('Original frozen spinup input mismatch')
        return original,om
    if mp.exists():
        m=json.loads(mp.read_text())
        if m['sources']!=sources() or digest(raw)!=m['raw_sha256']:raise RuntimeError('Spinup cache mismatch')
        return raw,m
    tic=time.time();s=DrivenSolver(n,cfg['nu'],cfg['injection'],backend);xp=s.xp
    h=random_field(s,seed,cfg['initial_rms'])[None];x=xp.zeros((1,0,3));f=xp.zeros((1,0,3,3))
    e0=.5*s.inner(h,h);budget=xp.zeros(1);t=0.;rows=[];maxcfl=0.
    for target in np.arange(0,cfg['spinup']+1e-9,.5):
        while t<target-1e-12:
            step=min(cfg['dt'],target-t);h,x,f,db,cf=s.step(h,x,f,step);budget+=db;t+=step;maxcfl=max(maxcfl,float(s.cpu(cf)[0]))
        row=stats(s,h,maxcfl);row.update(time=float(target),energy_budget_error=float(s.cpu((.5*s.inner(h,h)-e0-budget)/e0)[0]));rows.append(row)
    np.savez_compressed(raw,modes=s.modes_cpu,fields=s.pack(h))
    m=dict(seed=seed,regime=regime,n=n,status='complete',sources=sources(),raw_sha256=digest(raw),stats=rows,seconds=time.time()-tic)
    dump(mp,m);print(f'spinup {name}: {m["seconds"]:.1f}s Re_lambda={rows[-1]["Re_lambda"]:.2f} kmax_eta={rows[-1]["kmax_eta"]:.2f}',flush=True)
    if backend=='gpu':xp.get_default_memory_pool().free_all_blocks()
    return raw,m

def matched_translation(s,retained,scrambles):
    xp=s.xp;high=xp.max(abs(s.k),axis=0)>2;r=retained*high
    delta=scrambles-retained[None]
    target=float(s.cpu(xp.mean(xp.sqrt(s.inner(delta,delta)))))
    direction=np.array([1.,-.7,.9]);direction/=np.linalg.norm(direction)
    angles=sum(s.k[j]*direction[j] for j in range(3))
    def change(length):
        delta=r*(xp.exp(1j*length*angles)-1)
        return float(s.cpu(xp.mean(xp.sqrt(s.inner(delta,delta)))))
    samples=np.linspace(0,4*np.pi,513);values=[change(a) for a in samples]
    at=next((i for i in range(1,len(samples)) if values[i-1]<=target<=values[i]),None)
    if at is None:raise RuntimeError('No matched translation')
    length=brentq(lambda a:change(a)-target,samples[at-1],samples[at],xtol=1e-12)
    phase=xp.where(high,xp.exp(1j*length*angles),1)
    return retained*phase,dict(length=length,target_norm=target,achieved_norm=change(length),matching_error=abs(change(length)-target))

def run_one(seed,regime,layout,n,dt=None,backend='gpu',translation=None):
    cfg=P['regimes'][regime];dt=cfg['dt'] if dt is None else dt
    shift=np.zeros(3) if translation is None else np.array(translation)
    stem=f'{regime}-{layout}-s{seed}-n{n}-dt{dt:g}-{backend}'+('-translated' if translation is not None else '')
    folder=ROOT/'fast-data';folder.mkdir(exist_ok=True);raw=folder/(stem+'.npz');mp=raw.with_suffix('.json')
    if mp.exists():
        m=json.loads(mp.read_text())
        if m['sources']!=sources() or digest(raw)!=m['raw_sha256']:raise RuntimeError('Run cache mismatch')
        print(stem+' verified cache',flush=True);return m
    tic=time.time();s=DrivenSolver(n,cfg['nu'],cfg['injection'],backend);xp=s.xp
    spin_meta=None
    if cfg['spinup']:
        sr,spin_meta=spinup(seed,regime,backend)
        with np.load(sr) as z:initial=restrict(s,z['fields'],z['modes'])[0]
    else:initial=random_field(s,seed,cfg['initial_rms'])
    initial*=xp.exp(-1j*sum(s.k[j]*shift[j] for j in range(3)))
    writer_center=np.array([np.pi]*3)+np.array(P['layouts'][layout]['writer_shift'])+shift
    probe_center=np.array([np.pi]*3)+np.array(P['layouts'][layout]['probe_shift'])+shift
    q=basis(s,probe_center);writer=basis(s,writer_center,.45,4,True)[2]*P['writer_rms']
    x0=cloud(seed,probe_center);h=xp.stack([initial,initial]);x=xp.asarray(np.tile(x0,(2,1,1)));f=xp.asarray(np.tile(np.eye(3),(2,P['centers'],1,1)))
    force=xp.stack([writer,-writer]);e0=.5*s.inner(h,h);budget=xp.zeros(2);cfl=xp.zeros(2);pre_diag=[];t=0.
    for target in np.arange(0,P['prepare_end']+1e-9,P['save_dt']):
        while t<target-1e-12:
            step=min(dt,target-t)
            if t<P['write_stop']-1e-12:step=min(step,P['write_stop']-t);ff=force
            else:ff=xp.zeros_like(force)
            h,x,f,db,cf=s.step(h,x,f,step,ff);budget+=db;cfl=xp.maximum(cfl,cf);t+=step
        pre_diag.append(s.diagnostics(h,x,f,budget,e0,cfl))
    prepared=s.pack(h);variants,checks=interventions(s,h,seed,[.13,-.09,.07],P['scrambles'])
    if translation is not None:
        low=xp.max(abs(s.k),axis=0)<=2
        phase=xp.exp(-1j*sum(s.k[j]*shift[j] for j in range(3)))
        variants['reset']=variants['reset']*low+variants['reset']*(~low)*phase
    scr=xp.stack([variants['scramble'+str(i)] for i in range(P['scrambles'])]);shifted,match=matched_translation(s,variants['retained'],scr)
    variants['matched_translation']=shifted;variants.pop('sham')
    bases=xp.stack([variants[a][i] for a in P['arms'] for i in range(2)])
    c=s.cpu(xp.stack([s.inner(bases,q[j]) for j in range(3)],axis=1));coeff,audit=pulse_coefficients(c)
    states=[];labels=[]
    for i,(arm,hist) in enumerate((a,j) for a in P['arms'] for j in range(2)):
        states.append(bases[i]);labels.append(dict(arm=arm,history=hist,method='unprobed',base=i))
        for j,method in enumerate(P['probe_methods']):
            pulse=xp.einsum('i,icxyz->cxyz',xp.asarray(coeff[i,j]),q)
            states.append(bases[i]+pulse);labels.append(dict(arm=arm,history=hist,method=method,base=i))
    initial_branches=xp.stack(states);del states
    # Stream eight branches at a time at high resolution to bound working memory.
    batch=8 if n>=64 else len(labels);rows=[];ends=[]
    times=np.arange(0,P['horizon']+1e-9,P['save_dt'])
    for offset in range(0,len(labels),batch):
        h=initial_branches[offset:offset+batch].copy();b=len(h);x=xp.asarray(np.tile(x0,(b,1,1)));f=xp.asarray(np.tile(np.eye(3),(b,P['centers'],1,1)))
        e0=.5*s.inner(h,h);budget=xp.zeros(b);cfl=xp.zeros(b);t=0.;xx=[];ff=[];dd=[]
        for target in times:
            while t<target-1e-12:
                step=min(dt,target-t);h,x,f,db,cf=s.step(h,x,f,step);budget+=db;cfl=xp.maximum(cfl,cf);t+=step
            xx.append(s.cpu(x));ff.append(s.cpu(f));dd.append(s.diagnostics(h,x,f,budget,e0,cfl))
        rows.append((np.array(xx),np.array(ff),np.array(dd)));ends.append(s.pack(h))
    positions=np.concatenate([r[0] for r in rows],axis=1);tangent=np.concatenate([r[1] for r in rows],axis=1);diagnostics=np.concatenate([r[2] for r in rows],axis=1)
    end_fields=np.concatenate(ends);start_fields=s.pack(initial_branches)
    all_diag=np.concatenate([np.array(pre_diag).reshape(-1,7),diagnostics.reshape(-1,7)])
    maxima={name:float(np.max(abs(all_diag[:,j]))) if j==1 else float(np.max(all_diag[:,j])) for j,name in enumerate(DIAG)}
    field_stats=stats(s,s.unpack(prepared[:1]))
    np.savez_compressed(raw,modes=s.modes_cpu,prepared_fields=prepared,start_fields=start_fields,end_fields=end_fields,probe_basis=s.pack(q),pulse_coefficients=coeff,
                        time=times,positions=positions,tangent=tangent,diagnostics=diagnostics,prep_diagnostics=pre_diag,start_positions=x0)
    m=dict(status='complete',config=dict(seed=seed,regime=regime,layout=layout,n=n,dt=dt,backend=backend,translation=shift.tolist()),sources=sources(),raw_sha256=digest(raw),
           branches=labels,diagnostic_columns=DIAG,maxima=maxima,work_audit=audit,matching=match,intervention_checks=checks,field_stats=field_stats,
           spinup_raw_sha256=None if spin_meta is None else spin_meta['raw_sha256'],seconds=time.time()-tic,python=platform.python_version(),numpy=np.__version__)
    if backend=='gpu':m.update(cupy=xp.__version__,device=xp.cuda.runtime.getDeviceProperties(0)['name'].decode())
    dump(mp,m);print(f'{stem}: {m["seconds"]:.1f}s, CFL {maxima["max_cfl"]:.3f}, tail {maxima["tail_enstrophy"]:.4g}, work error {audit["equal_work_error"]:.2g}',flush=True)
    if backend=='gpu':xp.get_default_memory_pool().free_all_blocks()
    return m

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);ap.add_argument('--regime',default='smooth');ap.add_argument('--layout',default='coincident');ap.add_argument('--n',type=int,default=32);ap.add_argument('--dt',type=float);ap.add_argument('--backend',default='gpu');a=ap.parse_args()
    run_one(a.seed,a.regime,a.layout,a.n,a.dt,a.backend)
