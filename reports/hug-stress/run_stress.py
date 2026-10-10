"""Predeclared extreme matrix, seeded parameter exploration and adversarial controls.

Run from any directory. Failed numerical trials are retained as results.
The baseline equations and all original published files remain unchanged.
"""
from pathlib import Path
from dataclasses import asdict
from itertools import product
import json,hashlib,platform,time
import numpy as np
import scipy
from scipy.stats import qmc
from stress_solver import *

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'

def save(name,value):
    (DATA/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def case(name,group,p,q0=.04,v0=0.,end=24.):
    return dict(name=name,group=group,parameters=asdict(p),q0=q0,v0=v0,end=round(np.ceil(end/.02)*.02,10))

def design():
    cases=[]
    for i,(r,damping,load,duration,coupling) in enumerate(product((-4.,0.,4.),(0.,.05,4.,500.),(.5,1.0001,100.),(.004,4.),(0.,8.))):
        cases.append(case(f'matrix-{i:03d}','matrix',model.Parameters(r=r,damping=damping,pressure=load,pulse_duration=duration,memory_coupling=coupling)))
    # 144 factorial settings; every factor combination is explicit.
    points=qmc.LatinHypercube(d=7,seed=20261010).random(64)
    for i,u in enumerate(points):
        duration=10**(-3+4.7*u[3])
        p=model.Parameters(r=-10+20*u[0],damping=10**(-3+5.7*u[1]),pressure=10**(-3+5*u[2]),
                           pulse_duration=duration,memory_coupling=10**(-3+4.5*u[4]))
        q0=(-1 if i%2 else 1)*10**(-12+12.4*u[5]);v0=(-1 if i%3 else 1)*10**(-8+9*u[6])
        cases.append(case(f'seeded-{i:03d}','seeded',p,q0,v0,end=max(24,2*duration+12)))
    cases += [case('missed-pulse','adversarial',model.Parameters(pressure=100,pulse_duration=.001,memory_coupling=8),q0=0),
              case('stiff-damping','adversarial',model.Parameters(pressure=0,damping=500),q0=2),
              case('boundary-escape','adversarial',model.Parameters(r=4,pressure=100,memory_coupling=8,pulse_duration=40),q0=.04,end=120),
              case('boundary-no-feedback','adversarial',model.Parameters(r=-1,pressure=100,memory_coupling=0,pulse_duration=40),q0=0,end=120)]
    for amplitude in (1-1e-8,1.,1+1e-8):
        cases.append(case(f'threshold-{amplitude:.8f}','threshold',model.Parameters(pressure=amplitude,memory_coupling=0),q0=0))
    for i,q0 in enumerate((0.,-1e-12,1e-12,-1e-8,1e-8,-.04,.04)):
        cases.append(case(f'bias-{i}','bias',model.Parameters(r=1,pressure=0,damping=1),q0=q0,end=200))
    cases += [case('long-conservative','long',model.Parameters(r=-1,damping=0,pressure=0),q0=2,v0=3,end=1000),
              case('long-weak-damping','long',model.Parameters(r=1,damping=.001,pressure=0),q0=.02,v0=.3,end=1000)]
    return cases


def compute(c):
    p=model.Parameters(**c['parameters']);q0=c['q0'];v0=c['v0'];end=c['end']
    evaluate,ref_info=resolved(p,q0,v0,end)
    obs=observation_grid(p,end);ref=evaluate(obs)
    ext,ext_info=extent(obs,ref,p)
    arrays={'reference_time':obs,'reference_state':ref,'extent':ext}
    fixed_results=[]
    for tag,dt in [('coarse',.02),('fine',.01)]:
        t,y,failure=fixed(p,q0,v0,end,dt)
        arrays[tag+'_time']=t;arrays[tag+'_state']=y
        error=None if failure else scaled_error(y,evaluate(t))
        fixed_results.append(dict(dt=dt,failure=failure,scaled_reference_error=error,
                           numerical_pass=failure is None and error<=1e-3,energy_budget=energy_budget(y,p)))
    if any(r['failure'] for r in fixed_results):refinement=None
    else:
        common=arrays['coarse_time'];coarse=arrays['coarse_state'];fine=arrays['fine_state'][::2]
        scale=np.maximum(1,np.max(abs(evaluate(common)[:,:3]),axis=0))
        refinement=float(np.max(abs(coarse[:,:3]-fine[:,:3])/scale))
    result=dict(**c,reference=ref_info,reference_energy=energy_budget(ref,p),
                reference_peak_abs_q=float(max(abs(ref[:,0]))),reference_peak_imprint=float(max(ref[:,2])),
                final_q=float(ref[-1,0]),extent=ext_info,fixed=fixed_results,
                refinement=refinement,false_convergence=bool(refinement is not None and refinement<=1e-4 and fixed_results[1]['scaled_reference_error']>1e-3),
                memory_bounds_pass=bool(min(ref[:,2])>=-1e-8 and max(ref[:,2])<=p.pressure+1e-8))
    np.savez_compressed(DATA/(c['name']+'.npz'),**arrays)
    save(c['name']+'.json',result)
    return result


def independent_checks(rows):
    # Every adversarial/threshold/bias/long case, 16 spread through the screen,
    # and the 12 worst fine-RK4 errors. Selection is declared in protocol.
    selected={r['name'] for r in rows if r['group'] not in ('matrix','seeded')}
    screen=[r for r in rows if r['group'] in ('matrix','seeded')]
    selected.update(screen[i]['name'] for i in np.linspace(0,len(screen)-1,16,dtype=int))
    worst=sorted(screen,key=lambda r:1e30 if r['fixed'][1]['failure'] else r['fixed'][1]['scaled_reference_error'],reverse=True)[:12]
    selected.update(r['name'] for r in worst)
    comparisons=[]
    for i,row in enumerate(r for r in rows if r['name'] in selected):
        p=model.Parameters(**row['parameters']);base=np.load(DATA/(row['name']+'.npz'));obs=base['reference_time']
        evaluate,info=resolved(p,row['q0'],row['v0'],row['end'],method='Radau',rtol=2e-10)
        y=evaluate(obs);error=scaled_error(y,base['reference_state'])
        comparisons.append(dict(name=row['name'],scaled_discrepancy=error,passed=error<2e-6,**info))
        np.savez_compressed(DATA/(row['name']+'-radau.npz'),time=obs,state=y)
        if (i+1)%8==0:print(json.dumps(dict(stage='independent implicit solver',completed=i+1)),flush=True)
    save('independent-solvers.json',comparisons)
    return comparisons


def controls(rows):
    controls=[]
    def check(name,ok,**details):controls.append(dict(name=name,passed=bool(ok),**details))
    # Resolve the missed pulse explicitly at two resolutions, without changing equations.
    row=next(r for r in rows if r['name']=='missed-pulse');p=model.Parameters(**row['parameters'])
    evaluate,_=resolved(p,0,0,24)
    repair=[]
    for dt in (.02,.01):
        t,y,f=fixed(p,0,0,24,dt,resolve_pulse=True);err=scaled_error(y,evaluate(t))
        repair.append(dict(dt=dt,error=err,failure=f))
        np.savez_compressed(DATA/f'pulse-resolved-{dt}.npz',time=t,state=y)
        check('explicit pulse resolution '+str(dt),f is None and err<1e-6,error=err)
    save('pulse-resolution.json',repair)
    check('short pulse false convergence detected',row['false_convergence'],peak_imprint=row['reference_peak_imprint'])
    stiff=next(r for r in rows if r['name']=='stiff-damping')
    check('stiff numerical failure distinguished from bounded reference',all(r['failure'] is not None for r in stiff['fixed']) and stiff['reference_peak_abs_q']<3)
    thresholds=sorted([r for r in rows if r['group']=='threshold'],key=lambda r:r['parameters']['pressure'])
    check('opening threshold preserved exactly',thresholds[0]['extent']['opening_time'] is None and all(r['extent']['opening_time'] is not None for r in thresholds[1:]))
    bias=[r for r in rows if r['group']=='bias']
    check('exactly balanced state stays balanced',bias[0]['reference_peak_abs_q']==0)
    check('tiny positive and negative starts reach opposite rests',all(abs(abs(r['final_q'])-1)<1e-7 and np.sign(r['final_q'])==np.sign(r['q0']) for r in bias[1:]))
    # RHS equality under (q,v)->(-q,-v), with memory/work/loss unchanged.
    p=model.Parameters(r=4,pressure=100,memory_coupling=8);y=np.array([2.,-3.,.7,.4,.8]);M=np.array([-1,-1,1,1,1])
    err=float(max(abs(model.rhs(1,y*M,p)-model.rhs(1,y,p)*M)))
    check('algebraic reflection including work and loss',err<1e-12,error=err)
    # Compare analytic arm extents with a dense direct geometry sampling.
    rng=np.random.default_rng(5001);max_diff=0.
    for _ in range(20):
        p=model.Parameters(pressure=1.2);t=float(rng.uniform(0,4));y=np.array([rng.uniform(-20,20),0,rng.uniform(0,100),0,0])
        exact,_=extent(np.array([t]),y[None,:],p);g=model.geometry(t,y,p,points=4097)
        sampled=np.max(abs(np.array(g['left']+g['right'])))
        max_diff=max(max_diff,float(abs(exact[0]-sampled)))
    check('analytic geometry extent verified against dense arms',max_diff<1e-5,max_difference=max_diff)
    # Exact energy derivative, independently calculated from its gradient.
    errors=[]
    for _ in range(50):
        p=model.Parameters(r=float(rng.uniform(-10,10)),damping=float(rng.uniform(0,500)),memory_coupling=float(rng.uniform(0,30)))
        y=np.array([rng.normal(),rng.normal(),rng.uniform(0,100),0,0]);f=model.rhs(.2,y,p)
        derivative=(-p.r*y[0]+y[0]**3)*f[0]+y[1]*f[1]
        balance=f[3]-f[4];errors.append(abs(derivative-balance)/max(1,abs(balance)))
    check('exact lean energy identity',max(errors)<1e-11,maximum=max(errors))
    check('all adaptive memory trajectories respect pressure bounds',all(r['memory_bounds_pass'] for r in rows))
    save('controls.json',controls)
    return controls


def main():
    DATA.mkdir(exist_ok=True)
    cases=design()
    files=[p for p in (ROOT/'baseline').rglob('*.py') if '__pycache__' not in p.parts]
    provenance={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    fingerprint=hashlib.sha256((ROOT/'run_stress.py').read_bytes()+(ROOT/'stress_solver.py').read_bytes()).hexdigest()
    protocol=dict(baseline_commit='991704f2c4d708d2ab3fd29532ce176426ef78c9',source_hashes=provenance,
       runner_fingerprint=fingerprint,cases=cases,seed=20261010,
       thresholds=dict(scaled_reference_error=1e-3,apparently_converged=1e-4,independent_solvers=2e-6,cube_extent=3,guard=1e12),
       scaling='Each q,v,m component divided by max(1, its reference peak absolute value). This is not percent relative error near zero.',
       independent_selection='All adversarial, threshold, bias and long cases; 16 evenly spaced screen cases; 12 worst fine-step errors, numerical failures ranked first.',
       references='Split at prescribed pulse end, max_step <= pulse_duration/32 within pulse; DOP853 everywhere, independent Radau subset; both rtol 2e-10 and atol 2e-12.',
       scope='Same proposed reduced model; no fluid boundary, transport, wall force, contact or fracture added.')
    path=DATA/'protocol.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['runner_fingerprint']!=fingerprint:raise RuntimeError('Runner changed; use a fresh data directory to avoid mixing results')
    save('protocol.json',protocol)
    start=time.perf_counter();rows=[]
    for i,c in enumerate(cases):
        existing=DATA/(c['name']+'.json')
        if existing.exists():row=json.loads(existing.read_text())
        else:row=compute(c)
        rows.append(row)
        if (i+1)%16==0 or c['group'] not in ('matrix','seeded'):
            print(json.dumps(dict(stage='stress screen',completed=i+1,total=len(cases),name=c['name'],
                  failed_fine=sum(not r['fixed'][1]['numerical_pass'] for r in rows),false_convergence=sum(r['false_convergence'] for r in rows),
                  outside_cube=sum(r['extent']['outside_cube'] for r in rows),elapsed_seconds=round(time.perf_counter()-start,1))),flush=True)
    save('screen.json',rows)
    comparisons=independent_checks(rows);checks=controls(rows)
    summary=dict(configurations=len(rows),groups={g:sum(r['group']==g for r in rows) for g in sorted({r['group'] for r in rows})},
        fixed_attempts=2*len(rows),reference_runs=len(rows),independent_radau_runs=len(comparisons),
        coarse_pass=sum(r['fixed'][0]['numerical_pass'] for r in rows),fine_pass=sum(r['fixed'][1]['numerical_pass'] for r in rows),
        coarse_stopped=sum(r['fixed'][0]['failure'] is not None for r in rows),fine_stopped=sum(r['fixed'][1]['failure'] is not None for r in rows),
        false_convergence=sum(r['false_convergence'] for r in rows),outside_cube=sum(r['extent']['outside_cube'] for r in rows),
        maximum_extent=max(r['extent']['maximum'] for r in rows),independent_failures=[r for r in comparisons if not r['passed']],
        controls_passed=sum(c['passed'] for c in checks),controls_total=len(checks),control_failures=[c for c in checks if not c['passed']],
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),elapsed_seconds=time.perf_counter()-start)
    save('summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    if summary['independent_failures'] or summary['control_failures']:raise SystemExit(1)

if __name__=='__main__':main()
