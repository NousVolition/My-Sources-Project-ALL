"""Validate the added Lorenz-driven mode and retain its finite-time stress data."""
from pathlib import Path
from dataclasses import asdict
from itertools import product
import hashlib,json
import numpy as np
from lorenz_hug import *
from stress_solver import energy_budget

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'lorenz-data'
def save(name,value):(DATA/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def main():
    DATA.mkdir(exist_ok=True)
    configs=list(product((.5,10.,24.,24.74,28.,100.,160.),(.05,.6,4.),(0.,1.5),(.35,.65)))
    save('protocol.json',dict(configurations=len(configs),rho=[.5,10,24,24.74,28,100,160],damping=[.05,.6,4],coupling=[0,1.5],amplitude=[.35,.65],
        mean_rule='mean = amplitude + 0.10',pressure_rule='mean + amplitude*tanh((Z-max(rho-1,0))/10)',
        warmup=50,end=20,grid=.01,seed='deterministic initial [1,1,1]',
        solver=dict(method='DOP853',rtol=1e-10,atol=1e-12,max_step=.02),
        independent='Radau at the same tolerances for each displayed preset on [0,5]; pressure and opening checked at max_step .01 too.',
        scope='One-way Lorenz-to-pressure drive. Hug does not feed back into Lorenz. No claim of physical calibration.',
        sources={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('lorenz_hug.py','run_lorenz.py')}))
    starts={rho:warm_start(Drive(rho=rho)) for rho in sorted({x[0] for x in configs})}
    rows=[]
    for i,(rho,damping,coupling,amplitude) in enumerate(configs):
        d=Drive(rho=rho,mean=amplitude+.1,amplitude=amplitude);p=model.Parameters(r=-.25,damping=damping,memory_coupling=coupling)
        sol,opened,_=integrate(d,p,end=20,initial=starts[rho]);t=np.linspace(0,20,2001);y=sol.sol(t).T;load=pressure(y[:,2],d)
        row=dict(name=f'lorenz-{i:03d}',drive=asdict(d),hug=asdict(p),opening=opened,energy=energy_budget(y[:,3:],p),
                 peak_lean=float(max(abs(y[:,3]))),pressure_min=float(min(load)),pressure_max=float(max(load)),
                 memory_min=float(min(y[:,5])),memory_max=float(max(y[:,5])),final_lean=float(y[-1,3]),nfev=sol.nfev)
        np.savez_compressed(DATA/(row['name']+'.npz'),time=t,state=y,pressure=load);rows.append(row)
        if (i+1)%16==0:print('Lorenz-driven stress',i+1,'/',len(configs),flush=True)
    save('screen.json',rows)
    presets=[];checks=[]
    for key,label,rho,amplitude,coupling in [
        ('quiet','Below convection onset',.5,.35,1.5),('steady','Steady Lorenz input',10,.35,1.5),
        ('chaotic','Chaotic Lorenz pressure',28,.35,1.5),('no-feedback','Same drive, lean coupling off',28,.35,0),
        ('opening','Stronger chaotic pressure: opening',28,.65,1.5)]:
        d=Drive(rho=rho,mean=amplitude+.1,amplitude=amplitude);p=model.Parameters(r=-.25,memory_coupling=coupling)
        sol,opened,start=integrate(d,p,end=40,initial=starts[rho]);t=np.linspace(0,40,4001);y=sol.sol(t).T;load=pressure(y[:,2],d)
        np.savez_compressed(DATA/(key+'.npz'),time=t,state=y,pressure=load)
        ref,op2,_=integrate(d,p,end=5,initial=start,method='Radau')
        fine,op3,_=integrate(d,p,end=5,initial=start,max_step=.01,rtol=1e-11)
        obs=np.linspace(0,5,1001);base=sol.sol(obs).T;scale=np.maximum(1,np.max(abs(ref.sol(obs).T[:,:6]),axis=0))
        e=float(np.max(abs(base[:,:6]-ref.sol(obs).T[:,:6])/scale))
        e2=float(np.max(abs(base[:,:6]-fine.sol(obs).T[:,:6])/scale))
        events=[x for x in (opened,op2,op3) if x is not None and x<=5]
        event_error=float(max(events)-min(events)) if events else 0.
        opening_consistent=(len(events) in (0,3))
        checks.append(dict(name=key,short_independent_error=e,short_refinement_error=e2,opening_error=event_error,
                           passed=e<1e-6 and e2<1e-6 and event_error<1e-6 and opening_consistent))
        np.savez_compressed(DATA/(key+'-checks.npz'),time=obs,base=base,radau=ref.sol(obs).T,refined=fine.sol(obs).T)
        presets.append(dict(key=key,label=label,drive=asdict(d),hug=asdict(p),opening=opened,
            samples=np.round(np.column_stack([t[::10],y[::10,:6],load[::10]]),8).tolist()))
    save('presets.json',presets);save('solver-checks.json',checks)
    lyaps=[]
    for initial in ((1.,1.,1.),(1.001,1.,1.)):
        for tol,step in ((1e-8,.05),(1e-10,.025)):
            info,t,value=lyapunov(Drive(),initial=initial,rtol=tol,max_step=step)
            np.savez_compressed(DATA/f'lyapunov-{len(lyaps)}.npz',time=t,estimate=value)
            lyaps.append(info);print('Lorenz largest Lyapunov estimate',info['estimate'],flush=True)
    save('lyapunov.json',lyaps)
    # Analytic fixed-point linearization and the standard Hopf threshold.
    d=Drive();hopf=d.sigma*(d.sigma+d.beta+3)/(d.sigma-d.beta-1)
    stable=float(max(np.linalg.eigvals(lorenz_jac([np.sqrt(d.beta*23)]*2+[23],Drive(rho=24))).real))
    unstable=float(max(np.linalg.eigvals(lorenz_jac([np.sqrt(d.beta*24)]*2+[24],Drive(rho=25))).real))
    origin_stable=float(max(np.linalg.eigvals(lorenz_jac([0,0,0],Drive(rho=.5))).real))
    origin_unstable=float(max(np.linalg.eigvals(lorenz_jac([0,0,0],Drive(rho=2))).real))
    nofb=next(p for p in presets if p['key']=='no-feedback');ts=np.array(nofb['samples'])[:,0]
    # A separate hug integration with zero forcing has identical q,v when coupling=0.
    p=model.Parameters(r=-.25,memory_coupling=0,pressure=0)
    direct=model.reference(p,ts,q0=.04)
    independence=float(np.max(abs(np.array(nofb['samples'])[:,4:6]-direct[:,:2])))
    controls=dict(hopf_threshold=hopf,stable_real_eigenvalue_at_24=stable,unstable_real_eigenvalue_at_25=unstable,
        origin_stable_below_one=origin_stable,origin_unstable_above_one=origin_unstable,
        coupling_off_lean_error=independence,
        all_pressure_bounds=all(r['pressure_min']>=r['drive']['mean']-r['drive']['amplitude']-1e-12 and r['pressure_max']<=r['drive']['mean']+r['drive']['amplitude']+1e-12 for r in rows),
        all_memory_bounds=all(r['memory_min']>=-1e-10 and r['memory_max']<=r['drive']['mean']+r['drive']['amplitude']+1e-10 for r in rows),
        max_scaled_energy_residual=max(r['energy']['scaled'] for r in rows))
    controls['passed']=bool(stable<0<unstable and origin_stable<0<origin_unstable and independence<1e-7 and controls['all_pressure_bounds'] and controls['all_memory_bounds'] and controls['max_scaled_energy_residual']<1e-6)
    save('controls.json',controls)
    summary=dict(configurations=len(rows),displayed_presets=len(presets),independent_checks=checks,
                 controls=controls,lyapunov_min=min(x['estimate'] for x in lyaps),lyapunov_max=max(x['estimate'] for x in lyaps),
                 all_checks_passed=bool(controls['passed'] and all(x['passed'] for x in checks) and all(x['estimate']>0 for x in lyaps)))
    save('summary.json',summary);print(json.dumps(summary,indent=2))
    if not summary['all_checks_passed']:raise SystemExit(1)

if __name__=='__main__':main()
