"""Matched localized impulses in resolved periodic 3D Navier-Stokes fields."""
import argparse,json,time,hashlib
import numpy as np
from scipy.ndimage import map_coordinates
from physics import PilotFlow,initial,impulse,energy,L
from common import ROOT,save,digest,read_protocol

def markers(seed,clusters,edge):
    rng=np.random.default_rng(seed+830000);centers=rng.uniform(-L/2,L/2,(clusters,3))
    offsets=np.r_[np.zeros((1,3)),np.eye(3)*edge]
    return centers[:,None,:]+offsets[None,:,:]
def marker_sample(u,x,dx):
    coords=(((x.reshape(-1,3)+L/2)%L)/dx).T
    return np.stack([map_coordinates(a,coords,order=1,mode='grid-wrap',prefilter=False) for a in u],axis=-1).reshape(x.shape)
def geometry(x,edge):
    edges=x[:,1:]-x[:,:1];length=np.linalg.norm(edges,axis=-1)/edge
    logs=np.log(np.maximum(length,1e-12));volume=abs(np.linalg.det(edges))/edge**3
    lv=np.log(np.maximum(volume,1e-12));anis=logs.max(1)-logs.min(1)
    return np.array([logs.mean(),logs.std(),lv.mean(),lv.std(),anis.mean(),anis.std()])
def geometry_target(x,edge):
    return float(np.log(np.maximum(np.linalg.norm(x[:,1:]-x[:,:1],axis=-1)/edge,1e-12)).mean())
def stage(f,h,x):
    u=f.real(h);wh=f.curl(h);w=f.real(wh)
    adv=f.project(f.hat(np.cross(u,w,axisa=0,axisb=0,axisc=0)))
    rhs=adv-f.nu*f.k2*h+f.force
    rates=np.array([f.inner(h,f.force),f.nu*f.inner(wh,wh)])/L**3
    speed=float(np.sqrt((u*u).sum(0)).max())
    return rhs,marker_sample(u,x,f.dx),rates,speed
def advance(f,hs,xs,dt):
    kh=[];kx=[];rates=[];pd=[];speed=0.
    for j,a in enumerate((0.,.5,.5,1.)):
        hh=[h+a*dt*kh[-1][i] if j else h for i,h in enumerate(hs)]
        xx=[x+a*dt*kx[-1][i] if j else x for i,x in enumerate(xs)]
        stages=[stage(f,h,x) for h,x in zip(hh,xx)]
        kh.append([q[0] for q in stages]);kx.append([q[1] for q in stages]);rates.append(np.array([q[2] for q in stages]));speed=max(speed,max(q[3] for q in stages))
        row=[]
        for i in range(1,len(hs)):
            d=hh[i]-hh[0];diss=f.nu*f.inner(f.curl(d),f.curl(d))/L**3
            net=f.inner(d,kh[-1][i]-kh[-1][0])/L**3;row.append([net+diss,diss])
        pd.append(np.array(row).reshape(-1,2))
    weights=(1,2,2,1)
    hn=[f.project(h+dt*sum(w*kh[j][i] for j,w in enumerate(weights))/6) for i,h in enumerate(hs)]
    xn=[x+dt*sum(w*kx[j][i] for j,w in enumerate(weights))/6 for i,x in enumerate(xs)]
    return hn,xn,dt*sum(w*r for w,r in zip(weights,rates))/6,dt*sum(w*r for w,r in zip(weights,pd))/6,speed

def observe(f,h,x,edge):
    u=f.real(h);wh=f.curl(h);w=f.real(wh);adv=f.project(f.hat(np.cross(u,w,axisa=0,axisb=0,axisc=0)))
    e=energy(f,h);ens=.5*f.inner(wh,wh)/L**3;eps=2*f.nu*ens;urms=np.sqrt(2*e)
    rhs=adv-f.nu*f.k2*h+f.force;norm=lambda a:np.sqrt(max(f.inner(a,a)/L**3,0))
    forcecorr=f.inner(h,f.force)/max(np.sqrt(f.inner(h,h)*f.inner(f.force,f.force)),1e-30)
    shells=np.sqrt(f.k2)/(2*np.pi/L);total=max(f.inner(h,h),1e-30)
    bands=np.array([f.inner(h*(shells<2),h),f.inner(h*((shells>=2)&(shells<4)),h),f.inner(h*(shells>=4),h)])/total
    helicity=f.inner(h,wh)/L**3;speed=np.sqrt((u*u).sum(0));stretch=f.inner(wh,f.curl(adv))/L**3
    center=(f.n//2,)*3;localu=u[(slice(None),)+center];localw=w[(slice(None),)+center]
    modes=np.fft.fftfreq(f.n)*f.n;pos=np.fft.rfftfreq(f.n)*f.n
    phase=np.cos(np.pi*(modes[:,None,None]+modes[None,:,None]+pos[None,None,:]))
    grad=np.array([[float(np.sum((1j*k*h[i]*phase).real*f.weights)/f.n**3) for k in f.k] for i in range(3)])
    strain=(grad+grad.T)/2
    global8=np.array([e,ens,helicity,forcecorr,bands[0],speed.max(),np.mean(speed**4),norm(rhs)])
    local=np.array([np.linalg.norm(localu),np.linalg.norm(localw),np.sqrt((strain*strain).sum()),float(localw@strain@localw)])
    strong=np.r_[global8,bands[1:],stretch,local,geometry(x,edge)]
    integral=0.;spectrum=[]
    # Isotropic-spectrum integral-length convention is a diagnostic approximation here.
    for shell in range(1,int(np.ceil(shells.max()))+1):
        mask=(shells>=shell-.5)&(shells<shell+.5)
        shell_energy=f.inner(h*mask,h);integral+=shell_energy/(shell*2*np.pi/L);spectrum.append(shell_energy/(2*L**3))
    ell=3*np.pi/4*integral/total
    ucomp=urms/np.sqrt(3);lam=np.sqrt(15*f.nu*ucomp**2/max(eps,1e-30));relam=ucomp*lam/f.nu
    eta=(f.nu**3/max(eps,1e-30))**.25;cut=(int(np.ceil(f.n/3))-1)*2*np.pi/L
    low=shells<4;flux=-f.inner(h*low,adv)/L**3
    div=float(abs(f.real(1j*sum(k*v for k,v in zip(f.k,h)))).max())
    return strong,{'energy':e,'enstrophy':ens,'dissipation':eps,'urms':urms,'integral_length':ell,'integral_Re':urms*ell/f.nu,'Taylor_Re':relam,'kcut_eta':cut*eta,'high_fraction':f.inner(h*f.high,h)/total,'fraction_k_ge_4':float(bands[-1]),'vortex_stretching':stretch,'flux_k4':flux,'nonlinear_over_viscous_rms':norm(adv)/max(norm(f.nu*f.k2*h),1e-30),'divergence':div,'spectrum':spectrum}

def run(seed,regime,label='base',screen=False,n=None,dt=None):
    p=read_protocol(screen);c={**p,**p['regimes'][regime]};n=n or c['n'];dt=dt or c['dt']
    config={k:c[k] for k in ['nu','burn','end','save_dt','forcing_rms','initial_rms_range','kick_rms_over_baseline_rms','impulse_width','marker_clusters','marker_edge']};config.update(n=n,dt=dt,regime=regime)
    out=ROOT/'data'/('screen' if screen else 'runs');out.mkdir(parents=True,exist_ok=True);path=out/f'{regime}_{label}_{seed}.npz';meta=path.with_suffix('.json')
    sources={k:digest(ROOT/k) for k in ('simulate.py','physics.py','vendor/numerics.py','common.py')}
    if path.exists():
        old=json.loads(meta.read_text(encoding='utf-8'))
        if old['config']!=config or old['sources']!=sources or old['sha256']!=digest(path):raise ValueError('Stale cache '+str(path))
        print('Cached '+path.name,flush=True);return path
    start=time.perf_counter();f=PilotFlow(n,c['nu']);f.force=initial(f,7788,kind='beltrami')*c['forcing_rms']
    diffusion=float(f.nu*f.k2[f.keep].max()*dt)
    if diffusion>p['gates']['diffusion']:raise ValueError('Diffusion gate')
    speed=np.random.default_rng(seed+840000).uniform(*c['initial_rms_range']);h=initial(f,seed,speed);x=markers(seed,c['marker_clusters'],c['marker_edge'])
    every=round(c['save_dt']/dt);steps=round(c['burn']/dt)
    if not all(np.isclose(round(t/dt)*dt,t) for t in (c['save_dt'],c['burn'],c['end'])):raise ValueError('Times must align with time step')
    maxcfl=0.;history=[];ht=[];hdiag=[];budget=np.zeros((1,2));burn_e0=energy(f,h);burn_residual=[]
    for j in range(steps+1):
        if j%every==0:
            obs,d=observe(f,h,x,c['marker_edge']);history.append(obs);ht.append(j*dt-c['burn']);hdiag.append(d);burn_residual.append(abs(energy(f,h)-burn_e0-budget[0,0]+budget[0,1])/burn_e0)
        if j<steps:
            hh,xx,b,pd,v=advance(f,[h],[x],dt);h,x=hh[0],xx[0];budget+=b;maxcfl=max(maxcfl,v*dt/f.dx)
            if maxcfl>p['gates']['cfl']:raise FloatingPointError('Burn CFL gate')
    burn_balance=float(max(burn_residual))
    h0=h.copy();x0=x.copy();e0=energy(f,h);unit=impulse(f,1,c['impulse_width']);amps=np.array(c['kick_rms_over_baseline_rms'])*np.sqrt(2*e0)
    hs=[h0.copy()]+[h0+a*unit for a in amps];xs=[x0.copy() for _ in hs];ed0=np.array([energy(f,q-h0) for q in hs[1:]]);energies0=np.array([energy(f,q) for q in hs])
    initial_rates=[];base_rhs=f.rhs_only(h0)
    for q in hs[1:]:
        delta=q-h0;dd=f.nu*f.inner(f.curl(delta),f.curl(delta))/L**3;net=f.inner(delta,f.rhs_only(q)-base_rhs)/L**3;initial_rates.append([net+dd,dd])
    budgets=np.zeros((3,2));dbudget=np.zeros((2,2));times=[];ed=[];es=[];observations=[];diagnostics=[];positions=[];radii=[];spread=[];transfer=[];budget_trace=[];energy_residual=[];difference_residual=[]
    state_hash=lambda a:hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
    initial_hashes=[state_hash(q) for q in hs]
    axis=np.arange(n)*f.dx-L/2;xyz=np.meshgrid(axis,axis,axis,indexing='ij',sparse=True);r2=sum(a*a for a in xyz)
    steps=round(c['end']/dt)
    for j in range(steps+1):
        if j%every==0:
            times.append(j*dt);deltas=[q-hs[0] for q in hs[1:]];ed.append([energy(f,d) for d in deltas]);es.append([energy(f,q) for q in hs]);positions.append(np.stack(xs))
            oo=[];dd=[]
            for q,xx in zip(hs,xs):o,d=observe(f,q,xx,c['marker_edge']);oo.append(o);dd.append(d)
            observations.append(oo);diagnostics.append(dd);r=[];part=[]
            for delta in deltas:
                du=f.real(delta);weight=.5*(du*du).sum(0);s=weight.sum();r.append(np.sqrt((weight*r2).sum()/max(s,1e-30)));part.append(s*s/max((weight*weight).sum()*n**3,1e-30))
            radii.append(r);spread.append(part);transfer.append(dbudget.copy());budget_trace.append(budgets.copy())
            energy_residual.append(abs(np.array(es[-1])-energies0-budgets[:,0]+budgets[:,1])/max(e0,1e-30));difference_residual.append(abs(np.array(ed[-1])-ed0-dbudget[:,0]+dbudget[:,1])/ed0)
        if j<steps:
            hs,xs,b,pd,v=advance(f,hs,xs,dt);budgets+=b;dbudget+=pd;maxcfl=max(maxcfl,v*dt/f.dx)
            if maxcfl>p['gates']['cfl'] or not np.isfinite(np.stack(hs)).all():raise FloatingPointError('Response finite/CFL gate')
    es=np.array(es);ed=np.array(ed);energy_balance=np.array(energy_residual);diff_balance=np.array(difference_residual)
    all_diag=hdiag+[q for triplet in diagnostics for q in triplet]
    np.savez_compressed(path,time=times,history_time=ht,history=history,observations=observations,energy=es,perturbation_energy=ed,initial_perturbation_energy=ed0,positions=positions,radius=radii,participation=spread,transfer=transfer,energy_budget=budgets,energy_budget_trace=budget_trace,difference_budget=dbudget,energy_residual=energy_residual,difference_residual=difference_residual,initial_markers=x0,amplitudes=amps,relative_amplitudes=c['kick_rms_over_baseline_rms'],initial_rates=initial_rates)
    save(str(meta.relative_to(ROOT)),{'seed':seed,'config':config,'sources':sources,'sha256':digest(path),'seconds':time.perf_counter()-start,'initial_rms':speed,'initial_state_sha256':initial_hashes,'final_state_sha256':[state_hash(q) for q in hs],'max_cfl':maxcfl,'diffusion_number':diffusion,'burn_energy_balance':burn_balance,'max_energy_balance':float(energy_balance.max()),'max_difference_balance':float(diff_balance.max()),'max_divergence':max(q['divergence'] for q in all_diag),'min_kcut_eta':min(q['kcut_eta'] for q in all_diag),'max_high_fraction':max(q['high_fraction'] for q in all_diag),'history_diagnostics':hdiag,'diagnostics':diagnostics})
    print(f'{path.stem}: {time.perf_counter()-start:.1f}s; CFL={maxcfl:.3f}; min kcut_eta={min(q["kcut_eta"] for q in all_diag):.3f}; high E={max(q["high_fraction"] for q in all_diag):.2g}; difference balance={diff_balance.max():.2g}',flush=True);return path
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);ap.add_argument('--regime',choices=['smooth','active'],required=True);ap.add_argument('--label',default='base');ap.add_argument('--screen',action='store_true');ap.add_argument('--n',type=int);ap.add_argument('--dt',type=float);args=ap.parse_args();run(**vars(args))
