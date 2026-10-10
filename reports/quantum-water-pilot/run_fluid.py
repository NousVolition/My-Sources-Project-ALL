"""Measured-property-to-flow connection, separate from the molecular model."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded
from scipy.integrate import simpson
from iapws import IAPWS95,D2O
from iapws._iapws import _D2O_Viscosity
from vendor.ns_solver import Solver,difference

def dump(path,data):
    path.write_text(json.dumps(data,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))

def properties(T):
    return {name:dict(rho=s.rho,mu=float(s.mu),nu=float(s.mu/s.rho)) for name,c in [('H2O',IAPWS95),('D2O',D2O)] for s in [c(T=T,P=.101325)]}

def ablations(p):
    h,d=p['H2O'],p['D2O']
    return {'H2O':h,'density_only':dict(rho=d['rho'],mu=h['mu'],nu=h['mu']/d['rho']),
            'viscosity_only':dict(rho=h['rho'],mu=d['mu'],nu=d['mu']/h['rho']),'D2O':d}

def exact_startup(y,t,rho,mu,G,H):
    m=np.arange(1,1000,2)
    steady=G*y*(H-y)/(2*mu)
    transient=(4*G*H*H/(mu*np.pi**3))*np.sum(np.sin(m[:,None]*np.pi*y/H)*np.exp(-mu/rho*(m*np.pi/H)**2*t)[:,None]/m[:,None]**3,axis=0)
    return steady-transient

def channel(rho,mu,G,n=129,dt=.0005,end=1.,H=.001):
    y=np.linspace(0,H,n);dy=y[1]-y[0];nu=mu/rho;r=nu*dt/dy**2
    ab=np.zeros((3,n-2));ab[0,1:]=-r/2;ab[1]=1+r;ab[2,:-1]=-r/2
    u=np.zeros(n);pert=np.sin(np.pi*y/H)*1e-5
    energy0=0.;work=0.;diss=0.;times=[0.];trace=[[0.,0.,0.,0.]];profiles=[u.copy()]
    for j in range(round(end/dt)):
        old=u.copy();oldp=pert.copy()
        rhs=(1-r)*u[1:-1]+r/2*(u[:-2]+u[2:])+dt*G/rho
        u[1:-1]=solve_banded((1,1),ab,rhs)
        pert[1:-1]=solve_banded((1,1),ab,(1-r)*pert[1:-1]+r/2*(pert[:-2]+pert[2:]))
        mid=(old+u)/2
        # Discrete energy balance is exact for CN with midpoint work/dissipation.
        work+=dt*G*np.sum(mid[1:-1])*dy
        diss+=dt*mu*np.sum((np.diff(mid)/dy)**2)*dy
        if (j+1)%max(1,round(.01/dt))==0:
            t=(j+1)*dt;E=.5*rho*np.sum(u*u)*dy;eps=mu*np.sum((np.diff(u)/dy)**2)*dy
            times.append(t);trace.append([simpson(u,x=y)/H,E,eps,float(np.linalg.norm(pert)/np.linalg.norm(np.sin(np.pi*y/H)*1e-5))]);profiles.append(u.copy())
    exact=exact_startup(y,end,rho,mu,G,H)
    E=.5*rho*np.sum(u*u)*dy
    result=dict(rho=rho,mu=mu,G_Pam=G,n=n,dt_s=dt,H_m=H,mean_velocity_ms=simpson(u,x=y)/H,steady_mean_ms=G*H**2/(12*mu),center_velocity_ms=u[n//2],kinetic_energy_Jm2=E,dissipation_Wm2=mu*np.sum((np.diff(u)/dy)**2)*dy,integrated_work_Jm2=work,integrated_dissipation_Jm2=diss,energy_relative_residual=abs(E-energy0-work+diss)/work,relative_L2_error=np.linalg.norm(u-exact)/np.linalg.norm(exact),perturbation_amplitude_ratio=np.linalg.norm(pert)/np.linalg.norm(np.sin(np.pi*y/H)*1e-5),exact_perturbation_amplitude_ratio=np.exp(-nu*(np.pi/H)**2*end),tau_s=H*H/(np.pi*np.pi*nu),Re_bulk=simpson(u,x=y)*rho/mu)
    return result,dict(y_m=y,times_s=times,profiles_ms=profiles,trace_columns=np.array(['mean_velocity_ms','kinetic_energy_Jm2','dissipation_Wm2','perturbation_amplitude_ratio']),trace=trace)

def vortex(p,U=.01,n=16,dt=.01,seed=0):
    L=.001;nu=p['nu']/(U*L);solver=Solver(n,nu,workers=1)
    h=solver.initial('taylor_green',epsilon=.01 if seed else 0.,seed=seed or 1);h0=h.copy();e0=.5*solver.inner(h,h)
    duration=.1*U/L;steps=int(np.ceil(duration/dt));dt=duration/steps;loss=0.;maxcfl=0
    for j in range(steps):
        h,dl,cfl=solver.step(h,dt);loss+=dl;maxcfl=max(maxcfl,cfl)
    obs,_,_=solver.observe(h,duration,loss,e0)
    return dict(**obs,n=n,U_ms=U,dt_dimensionless=dt,seed=seed,nu_dimensionless=nu,seconds_physical=.1,rms_velocity_ms=U*np.sqrt(2*obs['energy']),kinetic_energy_Jm3=p['rho']*U*U*obs['energy'],dissipation_Wm3=p['rho']*U**3/L*obs['dissipation'],Re_scale=U*L/p['nu'],cfl=maxcfl),h,h0,solver

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    props={str(T):properties(T) for T in [283.15,298.15,313.15]};dump(a.out/'properties.json',props)
    checks=[(298.15,0,10.035938),(298.15,1105,1092.6424),(298.15,1130,1088.3626),(373.15,1064,326.63791),(775,1,29.639474),(775,100,31.930085),(775,400,53.324172)]
    verification=[dict(T=T,rho=rho,reference_microPas=ref,computed_microPas=_D2O_Viscosity(rho,T)*1e6) for T,rho,ref in checks]
    assert max(abs(x['computed_microPas']/x['reference_microPas']-1) for x in verification)<1e-7
    dump(a.out/'property_verification.json',verification)
    channels=[]
    for Ts,p in props.items():
        for name,material in ablations(p).items():
            for G in [1.,10.,100.]:
                result,raw=channel(material['rho'],material['mu'],G)
                result.update(T_K=float(Ts),material=name);channels.append(result)
                np.savez_compressed(a.out/f'channel-T{Ts}-{name}-G{G}.npz',**raw)
    convergence=[]
    for name,p in props['298.15'].items():
        for n,dt in [(33,.001),(65,.001),(129,.001),(129,.0005),(129,.00025)]:
            result,_=channel(p['rho'],p['mu'],10.,n,dt,end=.1);result['material']=name;convergence.append(result)
    dump(a.out/'channel_results.json',channels);dump(a.out/'channel_convergence.json',convergence)
    vortices=[];fields={}
    for Ts,p in props.items():
        for name,material in ablations(p).items():
            for U in [.005,.01,.02]:
                record,h,h0,s=vortex(material,U);tag=f'T{Ts}-{name}-U{U}-n16'
                record.update(T_K=float(Ts),material=name,tag=tag);vortices.append(record);fields[tag]=h
                np.savez_compressed(a.out/f'vortex-{tag}.npz',final_hat=h,initial_hat=h0)
    vchecks=[]
    for name,p in props['298.15'].items():
        for n,dt,seed in [(24,.01,0),(32,.01,0),(24,.005,0),(24,.01,1),(24,.01,2),(24,.01,3)]:
            record,h,h0,s=vortex(p,.01,n,dt,seed);tag=f'T298.15-{name}-U0.01-n{n}-dt{dt}-s{seed}'
            record.update(T_K=298.15,material=name,tag=tag);vortices.append(record);fields[tag]=h
            np.savez_compressed(a.out/f'vortex-{tag}.npz',final_hat=h,initial_hat=h0)
        h16=fields[f'T298.15-{name}-U0.01-n16'];h24=fields[f'T298.15-{name}-U0.01-n24-dt0.01-s0'];h32=fields[f'T298.15-{name}-U0.01-n32-dt0.01-s0'];hh=fields[f'T298.15-{name}-U0.01-n24-dt0.005-s0']
        vchecks.append(dict(material=name,n16_n24=difference(h16,h24,Solver(16,p['nu']/(.01*.001))),n24_n32=difference(h24,h32,Solver(24,p['nu']/(.01*.001))),dt_half=difference(h24,hh,Solver(24,p['nu']/(.01*.001)))))
    dump(a.out/'vortex_results.json',vortices);dump(a.out/'vortex_convergence.json',vchecks)
    # Exact 2D Taylor-Green solution verifies the reused nonlinear solver.
    s=Solver(16,.1,workers=1);h0=s.initial('taylor_green_2d');h=h0.copy()
    for j in range(100):h,_,_=s.step(h,.01)
    analytic=h0*np.exp(-.2)
    exact_error=np.sqrt(s.inner(h-analytic,h-analytic)/s.inner(analytic,analytic))
    assert exact_error<1e-9
    dump(a.out/'exact_vortex_check.json',dict(relative_error=exact_error))
    bench=[]
    for temp,observed in [(5.,1.3052),(20.,1.2452),(40.,1.2017)]:
        p=properties(temp+273.15);pred=p['D2O']['mu']/p['H2O']['mu']
        bench.append(dict(T_C=temp,observed_viscosity_ratio=observed,predicted_viscosity_ratio=pred,relative_error_percent=100*(pred/observed-1),measurement_derived_conductance_ratio=1/observed,predicted_conductance_ratio=1/pred))
    dump(a.out/'historical_benchmark.json',bench)
    rng=np.random.default_rng(20261009);p=props['298.15'];samples=100000
    hmu=p['H2O']['mu']*rng.uniform(.995,1.005,samples);dmu=p['D2O']['mu']*rng.uniform(.99,1.01,samples)
    hr=p['H2O']['rho']*rng.uniform(.9998,1.0002,samples);dr=p['D2O']['rho']*rng.uniform(.9998,1.0002,samples)
    sensitivity=dict(assumptions='Independent rectangular property ranges, sensitivity only, not calibrated confidence.',steady_velocity_D_over_H_quantiles=np.quantile(hmu/dmu,[.025,.5,.975]),decay_time_D_over_H_quantiles=np.quantile((dr/dmu)/(hr/hmu),[.025,.5,.975]))
    dump(a.out/'sensitivity.json',sensitivity)
    print('Completed',len(channels),'channel cases,',len(vortices),'3D vortex cases and convergence/benchmark controls.',flush=True)

if __name__=='__main__': main()
