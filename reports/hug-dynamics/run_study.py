"""Run the fixed hug extension protocol, references and textbook controls."""
from pathlib import Path
from dataclasses import asdict
from itertools import product
import importlib.util
import json
import hashlib
import sys
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from hug_model import *

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
CHECKS=[]


def check(name, condition, value=None):
    CHECKS.append(dict(name=name,passed=bool(condition),value=value))


def save(name,value):
    (DATA/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def source_checks():
    expected={'hug_envelope.py':'27df38edb3c267e619b3685d6ff6a62a13abd6f4',
              'pressure_envelope.py':'5cf93cc02246d58e0a008d779617a84f2e4a99c5',
              'gated_fluid.py':'9af57c746a2e3aa872836626ea01660ad7a50662'}
    hashes={}
    for name,sha in expected.items():
        b=(HERE/'source'/name).read_bytes().replace(b'\r\n',b'\n')
        actual=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        check('pinned GitHub source '+name, actual==sha,actual)
        hashes[name]=dict(git_blob=actual,sha256_lf=hashlib.sha256(b).hexdigest())
    for filename in ('test_original_hug.py','test_original_pressure.py'):
        spec=importlib.util.spec_from_file_location(filename,HERE/'tests'/filename)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        for name in dir(module):
            if name.startswith('test_'):
                try:getattr(module,name)();check(filename+':'+name,True)
                except Exception as exc:check(filename+':'+name,False,repr(exc))
    from pressure_envelope import PressureEnvelope
    model=PressureEnvelope();state=model.step(.8,2.)
    check('original constant-load analytic memory',abs(state['imprint']-.8*(1-np.exp(-2/.8)))<1e-14)
    y=np.array([0,0,state['imprint'],0,0])
    g=geometry(2,y,Parameters(pressure=.8,memory_coupling=0))
    delta=float(np.max(np.abs(np.array(g['left'])-np.array(state['Norm']))))
    check('zero added lean recovers original arm geometry',delta<1e-14,delta)
    import gated_fluid as gf
    edge=float(gf.gate(np.array([2.8]),np.array([0.]),np.array([0.]))[0])
    core=float(gf.gate(np.array([0.]),np.array([0.]),np.array([0.]))[0])
    check('original corrected cube gate core and collar',core==1 and edge==0)
    save('sources.json',dict(repository='NousVolition/My-Sources-Project-ALL',
         commit='f5a2106f69bd1c4a0f8eb7a64f5291952c4c5e25',files=hashes))


def sweeps():
    rows=[]
    for i,(r,damping,load,q0) in enumerate(product((-1.,0.,1.),(.2,1.,4.),(0.,.8,1.2),(-.04,0.,.04))):
        p=Parameters(r=r,damping=damping,pressure=load)
        t,y=integrate(p,q0,dt=.02)
        tf,yf=integrate(p,q0,dt=.01)
        e=energy(y,p);residual=e-e[0]-y[:,3]+y[:,4]
        difference=float(np.max(abs(y[:,:3]-yf[::2,:3])))
        row=dict(id=i,**asdict(p),q0=q0,final_q=float(y[-1,0]),peak_abs_q=float(max(abs(y[:,0]))),
            max_step_difference=difference,energy_budget_error=float(max(abs(residual))),opening_time=opening_time(p))
        rows.append(row)
        check(f'sweep {i} time refinement',difference<2e-4,difference)
        check(f'sweep {i} oscillator energy budget',row['energy_budget_error']<2e-5,row['energy_budget_error'])
        np.savez_compressed(DATA/f'sweep-{i:03d}.npz',time=t,state=y,fine_time=tf,fine_state=yf)
    mirror=max(abs(rows[i]['final_q']+rows[i+2]['final_q']) for i in range(0,len(rows),3))
    zero=max(abs(r['final_q']) for r in rows if r['q0']==0)
    check('81-case reflected endpoints',mirror<1e-12,mirror)
    check('27 exactly balanced starts remain balanced',zero==0,zero)
    save('sweep.json',rows)
    return rows


def presets():
    visual=[];summary=[]
    for key,(label,p,q0) in PRESETS.items():
        t,y=integrate(p,q0,dt=.01)
        yr=reference(p,t,q0)
        yr2=reference(p,t,q0,tolerance=1e-12)
        err=float(np.max(abs(y[:,:3]-yr[:,:3])))
        referr=float(np.max(abs(yr-yr2)))
        check(key+' independent solver agreement',err<2e-5,err)
        check(key+' reference tolerance refinement',referr<2e-7,referr)
        maxgeom=0.;mirror=0.
        for ti,yi in zip(t[::10],y[::10]):
            g=geometry(ti,yi,p)
            maxgeom=max(maxgeom,float(np.max(np.abs(g['left']+g['right']))))
            ym=yi.copy();ym[:2]*=-1;gm=geometry(ti,ym,p)
            target=np.array(g['left'])*[-1,1]
            mirror=max(mirror,float(np.max(abs(target-np.array(gm['right'])))))
        check(key+' geometric reflection',mirror<1e-12,mirror)
        check(key+' stays inside chosen cube guide',maxgeom<3,maxgeom)
        check(key+' unit determinant for closed area',all(abs(geometry(ti,yi,p)['sx']*geometry(ti,yi,p)['sy']-1)<1e-14 for ti,yi in zip(t[::100],y[::100])))
        samples=[]
        for ti,yi in zip(t[::10],y[::10]):
            g=geometry(ti,yi,p)
            samples.append([round(float(v),7) for v in [ti,yi[0],yi[1],yi[2],pressure(ti,p),g['sx'],g['bend']]])
        visual.append(dict(key=key,label=label,parameters=asdict(p),q0=q0,opening_time=opening_time(p),samples=samples))
        summary.append(dict(key=key,label=label,final_q=float(y[-1,0]),final_memory=float(y[-1,2]),
                            peak_abs_q=float(max(abs(y[:,0]))),reference_error=err,reference_refinement=referr))
        np.savez_compressed(DATA/f'preset-{key}.npz',time=t,state=y,reference=yr)
    save('presets.json',visual);save('preset-checks.json',summary)
    # Verify that the zero-memory model's stable equilibria match its actual long runs.
    for sign in (-1,1):
        p=Parameters(r=1,damping=1,pressure=0)
        _,y=integrate(p,sign*.04,end=60)
        check(f'buckling equilibrium sign {sign}',abs(y[-1,0]-sign)<1e-9,float(y[-1,0]))
    return summary


def textbook_controls():
    numerics=[]
    for method,order in [('euler',1),('heun',2),('rk4',4)]:
        errors=[]
        for h in (.2,.1,.05,.025):
            y=np.array([1.])
            for i in range(round(1/h)):y=step(lambda t,z:-z,i*h,y,h,method)
            errors.append(float(abs(y[0]-np.exp(-1))))
        observed=float(np.log2(errors[-2]/errors[-1]))
        check(method+' global convergence order',abs(observed-order)<.12,observed)
        local=[]
        for h in (.4,.2,.1):
            local.append(float(abs(step(lambda t,z:-z,0,np.array([1.]),h,method)[0]-np.exp(-h))))
        localorder=float(np.log2(local[-2]/local[-1]))
        check(method+' one-step local error order',abs(localorder-order-1)<.15,localorder)
        numerics.append(dict(method=method,steps=[.2,.1,.05,.025],errors=errors,global_order=observed,
                             local_steps=[.4,.2,.1],local_errors=local,local_order=localorder))
    save('numerical-orders.json',numerics)
    # An independent hug comparison; switching memory can lower formal RK order.
    comparison=[];p=Parameters(r=-.25,memory_coupling=1.5)
    for method in ('euler','heun','rk4'):
        for dt in (.08,.04,.02,.01):
            t,y=integrate(p,dt=dt,method=method);yr=reference(p,t)
            comparison.append(dict(method=method,dt=dt,max_error=float(np.max(abs(y[:,:3]-yr[:,:3])))))
    save('hug-integrators.json',comparison)
    check('hug RK4 refinement improves agreement',comparison[-1]['max_error']<comparison[-4]['max_error'])
    for method in ('euler','heun','rk4'):
        errors=[r['max_error'] for r in comparison if r['method']==method]
        check(method+' hug refinement reduces every discrepancy',all(a>b for a,b in zip(errors,errors[1:])),errors)
    check('method comparison actually separates the methods',len({r['max_error'] for r in comparison if r['dt']==.04})==3)
    pend=[]
    for b in (0.,.2,2.):
        for v0 in (0.,2.2):
            ts=np.linspace(0,24,1201)
            sol=solve_ivp(lambda t,z:[z[1],-np.sin(z[0])-b*z[1],b*z[1]**2],(0,24),[.4,v0,0],
                          t_eval=ts,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05)
            assert sol.success
            e=.5*sol.y[1]**2-np.cos(sol.y[0]);err=float(max(abs(e-e[0]+sol.y[2])))
            check(f'pendulum energy identity b={b} v={v0}',err<1e-8,err)
            pend.append(dict(damping=b,v0=v0,initial_energy=float(e[0]),energy_budget_error=err,
                             angular_turns=float((sol.y[0,-1]-sol.y[0,0])/(2*np.pi))))
            np.savez_compressed(DATA/f'pendulum-{b}-{v0}.npz',time=ts,state=sol.y.T)
    save('pendulum.json',pend)
    stability=[]
    for r in (-1.,0.,1.):
        for q,eigen in equilibria(r,1):
            index=winding(lambda z:np.array([z[1],r*z[0]-z[0]**3-z[1]]),np.array([q,0.]))
            stability.append(dict(r=r,q=q,eigenvalues=[[float(v.real),float(v.imag)] for v in eigen],index=index,
                   classification='nonhyperbolic; nonlinear restoring potential' if r==0 else ('stable' if max(eigen.real)<0 else 'saddle')))
            check(f'index r={r} q={q}',abs(index-(-1 if r>0 and q==0 else 1))<1e-12,index)
    try:winding(lambda z:z,np.array([.1,0]),.1);check('index refuses a zero on curve',False)
    except ValueError:check('index refuses a zero on curve',True)
    save('stability.json',stability)
    # Nonhyperbolic scalar examples, evaluated on both sides of zero.
    scalar=[]
    for name,f,expected in [('negative cubic',lambda x:-x**3,'stable'),('positive cubic',lambda x:x**3,'unstable'),
                            ('square',lambda x:x*x,'one-sided'),('zero',lambda x:0*x,'neutral')]:
        scalar.append(dict(name=name,derivative_at_zero=0,left=float(f(-.1)),right=float(f(.1)),nonlinear_classification=expected))
    save('degenerate.json',scalar)
    # These candidate feedback laws are separate scalar tests, not water chemistry.
    growth=[]
    laws={'logistic':lambda x:x*(1-x),'autocatalytic':lambda x:1.*1.*x-1.*x*x,
          'gompertz':lambda x:-x*np.log(x),'allee':lambda x:x*(1-x)*(x-.3)}
    for name,f in laws.items():
        for x0 in (.1,.5,1.5):
            ts=np.linspace(0,30,301)
            sol=solve_ivp(lambda t,z:[f(z[0])],(0,30),[x0],t_eval=ts,rtol=1e-10,atol=1e-12,method='DOP853')
            assert sol.success
            if name in ('logistic','autocatalytic'):exact=1/(1+(1/x0-1)*np.exp(-ts))
            elif name=='gompertz':exact=np.exp(np.log(x0)*np.exp(-ts))
            else:exact=None
            err=float(max(abs(sol.y[0]-exact))) if exact is not None else None
            if err is not None:check(f'{name} analytic solution x0={x0}',err<2e-8,err)
            check(f'{name} positivity x0={x0}',min(sol.y[0])>0)
            growth.append(dict(law=name,x0=x0,final=float(sol.y[0,-1]),analytic_error=err,
                               samples=[[round(float(t),5),round(float(x),8)] for t,x in zip(ts,sol.y[0])]))
    save('growth.json',growth)
    # Infinite delayed-departure solutions; no solver is allowed to choose a unique one.
    nonunique=[]
    for delay in (0.,1.,2.):
        for sign in (-1,1):
            ts=np.linspace(0,4,401);elapsed=np.maximum(ts-delay,0)
            x=sign*(2*elapsed/3)**1.5;derivative=sign*np.sqrt(2*elapsed/3)
            err=float(max(abs(derivative-np.cbrt(x))))
            check(f'delayed cube-root solution {delay} {sign}',err<1e-14 and x[0]==0,err)
            nonunique.append(dict(delay=delay,sign=sign,residual=err))
    save('nonuniqueness.json',dict(family=nonunique,stationary_solution=0,
        conclusion='Same initial value permits stationary and arbitrarily delayed departures in either direction. Excluded from hug law.'))
    # Linear overdamped comparison with consistent initial slow velocity.
    slow=[]
    for mass in (1.,.2,.04):
        ts=np.linspace(0,5,501)
        sol=solve_ivp(lambda t,z:[z[1],(-z[1]-z[0])/mass],(0,5),[1,-1],t_eval=ts,rtol=1e-11,atol=1e-13)
        assert sol.success
        slow.append(dict(mass=mass,error=float(max(abs(sol.y[0]-np.exp(-ts))))))
    check('overdamped limit improves as mass decreases',slow[2]['error']<slow[1]['error']<slow[0]['error'])
    save('overdamped.json',slow)
    return numerics


def main():
    DATA.mkdir(exist_ok=True)
    # Protocol is saved before executing trials. Deterministic repeats are not independent samples.
    save('protocol.json',dict(model='reduced pressure-memory and lean extension',end=24,steps=[.02,.01],
       r=[-1,0,1],damping=[.2,1,4],pressure=[0,.8,1.2],starting_lean=[-.04,0,.04],
       thresholds=dict(sweep_step_difference=2e-4,energy_budget=2e-5,preset_reference=2e-5),
       limitations=['No fluid coupling, permeability or calibrated material law','No physical seconds, pressure or temperature calibration',
                   'Gate display is an initial support guide, not an active wall','Numerical agreement is not experimental validation']))
    source_checks();rows=sweeps();summary=presets();orders=textbook_controls()
    failures=[x for x in CHECKS if not x['passed']]
    report=dict(checks=len(CHECKS),passed=len(CHECKS)-len(failures),failed=len(failures),details=CHECKS,
       runs=dict(sweep_configurations=len(rows),sweep_integrations=2*len(rows),presets=8,preset_reference_integrations=16,
                 integrator_comparisons=12,integrator_reference_integrations=12,pendulum=6,growth=12,overdamped=3,buckling_long=2),
       maxima=dict(sweep_step_difference=max(r['max_step_difference'] for r in rows),
                   sweep_energy_budget=max(r['energy_budget_error'] for r in rows),
                   preset_reference=max(s['reference_error'] for s in summary)),
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    save('verification.json',report)
    print(json.dumps({k:report[k] for k in ['checks','passed','failed','runs','maxima','environment']},indent=2))
    if failures:print(json.dumps(failures,indent=2));sys.exit(1)


if __name__=='__main__':main()
