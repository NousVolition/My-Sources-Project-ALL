"""Run controls first, then matched two-way comparisons and finite stress tests."""
from pathlib import Path
from dataclasses import asdict
import json,hashlib,platform
import numpy as np
import scipy
from scipy.integrate import trapezoid
import two_way_hug as tw
from lorenz_hug import integrate as old_integrate, pressure, Drive, model

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'two-way-data'
def save(name,value):
    (DATA/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def sample(sol,d,end,spacing=.02):
    t=np.linspace(0,end,round(end/spacing)+1);y=sol.sol(t).T
    load=pressure(y[:,2],d)
    return t,y,load

def metrics(t,y,load,sol,d,p,g,burn):
    mask=t>=burn;tt=t[mask];yy=y[mask]
    mean=lambda value:float(trapezoid(value,tt)/(tt[-1]-tt[0]))
    kinetic=.5*np.sum(y[:,:3]**2,axis=1)
    residual=kinetic-kinetic[0]-y[:,9]+y[:,8]
    scale=max(1.,float(np.max(abs(kinetic))),float(np.max(abs(y[:,9]))),float(y[-1,8]))
    activity=np.sum(yy[:,:2]**2,axis=1)
    qmean=mean(yy[:,3]);zmean=mean(yy[:,2])
    loss=np.array([tw.resistance(a,b,g,sol.opening) for a,b in zip(t,y[:,5])])*np.sum(y[:,:2]**2,axis=1)
    return dict(activity_rms=float(np.sqrt(mean(activity))),mean_z=zmean,
        z_std=float(np.sqrt(max(0.,mean((yy[:,2]-zmean)**2)))),mean_pressure=mean(load[mask]),
        mean_imprint=mean(yy[:,5]),lean_rms=float(np.sqrt(mean(yy[:,3]**2))),mean_lean=qmean,
        opening=sol.opening,feedback_loss=float(y[-1,8]),max_feedback_power=float(max(loss)),
        scaled_activity_budget_error=float(max(abs(residual))/scale),
        pressure_bounds=bool(min(load)>=d.mean-d.amplitude-1e-12 and max(load)<=d.mean+d.amplitude+1e-12),
        memory_bounds=bool(min(y[:,5])>=-1e-9 and max(y[:,5])<=d.mean+d.amplitude+1e-9),
        passive_feedback=bool(min(loss)>=-1e-12),peak_state=float(np.max(abs(y[:,:6]))),nfev=sol.nfev)

def main():
    DATA.mkdir(exist_ok=True)
    d=Drive();p=model.Parameters(r=-.25,memory_coupling=1.5)
    initial=np.load(ROOT/'lorenz-data/chaotic.npz')['state'][0,:3]
    controls=[]
    # Independent numerical methods and refinement on the same short interval.
    for name,drive,g,params in [('off',d,0.,p),('weak',d,.5,p),('medium',d,2.,p),
        ('strong',d,8.,p),('extreme',d,100.,p),('opening',Drive(mean=.7,amplitude=.6),8.,p),
        ('stiff-hug',d,8.,model.Parameters(r=-.25,damping=500.,memory_coupling=1.5))]:
        start=(1.,1.,1.) if name=='opening' else initial
        base=tw.integrate(drive,params,g,end=5,initial=start)
        ref=tw.integrate(drive,params,g,end=5,initial=start,method='Radau')
        fine=tw.integrate(drive,params,g,end=5,initial=start,rtol=1e-11,max_step=.01)
        t,y,load=sample(base,drive,5,.005);r=ref.sol(t).T;f=fine.sol(t).T
        scale=np.maximum(1,np.max(abs(r[:,:6]),axis=0))
        e=float(np.max(abs(y[:,:6]-r[:,:6])/scale));ef=float(np.max(abs(y[:,:6]-f[:,:6])/scale))
        openings=[s.opening for s in (base,ref,fine)]
        consistent=all(x is None for x in openings) or all(x is not None for x in openings)
        ee=float(max(openings)-min(openings)) if all(x is not None for x in openings) else 0.
        row=dict(name=name,independent_error=e,refinement_error=ef,opening_error=ee,
            **metrics(t,y,load,base,drive,params,g,0))
        if g==0:
            old,_,_=old_integrate(drive,params,end=5,initial=start)
            row['original_one_way_error']=float(np.max(abs(y[:,:6]-old.sol(t).T[:,:6])/scale))
        if name=='opening':
            after=t>=base.opening+.8
            row['released_feedback_power']=float(max(tw.resistance(a,b,g,base.opening) for a,b in zip(t[after],y[after,5])))
        row['passed']=bool(e<1e-5 and ef<1e-5 and consistent and ee<1e-6 and
            row['scaled_activity_budget_error']<1e-7 and row['memory_bounds'] and row['passive_feedback'] and
            row.get('original_one_way_error',0)<1e-6 and row.get('released_feedback_power',0)==0)
        controls.append(row)
        np.savez_compressed(DATA/(name+'-checks.npz'),time=t,base=y,radau=r,refined=f)
        print('Control',name,row['passed'],e,flush=True)
    save('controls.json',controls)
    if not all(c['passed'] for c in controls):raise RuntimeError('Controls failed; sweep not started')
    # Suppressing the forward pressure modulation makes hug motion independent of Lorenz initial state.
    const=Drive(mean=.45,amplitude=0)
    a=tw.integrate(const,p,2,end=5,initial=initial);b=tw.integrate(const,p,2,end=5,initial=(1,1,1))
    ts=np.linspace(0,5,501)
    forward_off=float(np.max(abs(a.sol(ts)[3:6]-b.sol(ts)[3:6])))
    assert forward_off<1e-7
    sources={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in
        ['two_way_hug.py','run_two_way.py','lorenz_hug.py','stress_solver.py','hug_model.py'] if (ROOT/f).exists()}
    save('protocol.json',dict(rule='eta = g * bend(t)^2 * imprint/(1+imprint); add -eta*X and -eta*Y',
        forward_rule='Existing bounded Z-to-pressure mapping; pressure drives imprint and lean',
        sources=sources,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        method='DOP853',rtol=1e-10,max_step=.02,refined_rtol=1e-11,refined_max_step=.01,
        control_target=1e-5,end=80,burn=20,record_spacing=.02,feedback_strengths=[0,.5,2,8],
        scope='Proposed dissipative feedback; not a physical calibration or conservative energy exchange',
        constant_pressure_hug_independence_error=forward_off))
    rows=[]
    for seed,start in enumerate([initial,np.array([-8.,8.,27.]),np.array([1.,1.,1.])]):
        for g in [0.,.5,2.,8.]:
            name=f'paired-{seed}-g{g:g}'
            base=tw.integrate(d,p,g,initial=start);fine=tw.integrate(d,p,g,initial=start,rtol=1e-11,max_step=.01)
            t,y,load=sample(base,d,80);tf,yf,lf=sample(fine,d,80)
            row=dict(name=name,seed=seed,strength=g,initial=start.tolist(),drive=asdict(d),hug=asdict(p),
                metrics=metrics(t,y,load,base,d,p,g,20),refined=metrics(tf,yf,lf,fine,d,p,g,20))
            rows.append(row)
            np.savez_compressed(DATA/(name+'.npz'),time=t,state=y,pressure=load,refined=yf,refined_pressure=lf)
            print('Paired',name,'activity',row['metrics']['activity_rms'],flush=True)
    save('paired.json',rows)
    stress=[]
    for rho in [10.,28.,100.]:
        for g in [2.,20.]:
            for damping in [.05,50.]:
                drive=Drive(rho=rho);params=model.Parameters(r=-.25,damping=damping,memory_coupling=1.5)
                sol=tw.integrate(drive,params,g,end=20);t,y,load=sample(sol,drive,20)
                row=dict(name=f'stress-{len(stress):02}',rho=rho,strength=g,damping=damping,
                    metrics=metrics(t,y,load,sol,drive,params,g,5))
                stress.append(row);np.savez_compressed(DATA/(row['name']+'.npz'),time=t,state=y,pressure=load)
    save('stress.json',stress)
    summaries=[]
    for g in [.5,2.,8.]:
        paired=[r for r in rows if r['strength']==g]
        effect={}
        for metric in ['activity_rms','z_std','mean_pressure','mean_imprint','lean_rms']:
            changes=[];uncertainty=[]
            for r in paired:
                baseline=next(b for b in rows if b['seed']==r['seed'] and b['strength']==0)
                changes.append(100*(r['metrics'][metric]/baseline['metrics'][metric]-1))
                uncertainty.append(100*(r['refined'][metric]/baseline['refined'][metric]-1))
            effect[metric]=dict(percent_changes=changes,refined_percent_changes=uncertainty,
                mean_percent=float(np.mean(changes)),range_percent=[min(changes),max(changes)],
                largest_refinement_change_pp=float(np.max(abs(np.array(changes)-uncertainty))))
        summaries.append(dict(strength=g,effects=effect))
    allmetrics=[r['metrics'] for r in rows+stress]+[r['refined'] for r in rows]
    summary=dict(controls_passed=True,paired_configurations=len(rows),paired_solver_runs=2*len(rows),
        stress_configurations=len(stress),short_three_solver_controls=len(controls),effects=summaries,
        constant_pressure_hug_independence_error=forward_off,
        all_bounds_and_passivity=all(m['memory_bounds'] and m['pressure_bounds'] and m['passive_feedback'] for m in allmetrics),
        maximum_scaled_budget_residual=max(m['scaled_activity_budget_error'] for m in allmetrics),
        note='Ranges are finite-window sensitivity across three starts, not confidence intervals. Chaotic pointwise agreement is not required at long times.')
    save('summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    assert summary['all_bounds_and_passivity'] and summary['maximum_scaled_budget_residual']<1e-7

if __name__=='__main__':main()
