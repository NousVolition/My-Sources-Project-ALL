"""Reproducible bistability/hysteresis tests of the supplied Hill-feedback ODE.

Assumption: A is clamped, not a conserved pool with Ap. Then x=Ap/K,
tau=kd*t, s=kp*S*A/(kd*K), a=beta/(kd*K): x'=s+a*x**n/(1+x**n)-x.
Chosen dimensionless parameters illustrate mechanisms, not a fit to biology.
"""
from pathlib import Path
import argparse, json, platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent


def rhs(x, s, a, n=2):
    return s+a*x**n/(1+x**n)-x


def slope(x, a, n=2):
    return a*n*x**(n-1)/(1+x**n)**2-1


def equilibria(s, a, n=2):
    if n == 2:
        roots = np.roots([1, -(s+a), 1, -s])
    elif n == 1:
        roots = np.roots([1, 1-s-a, -s])
    else:
        raise ValueError('Equilibrium polynomial implemented for n=1 or n=2.')
    real = sorted(max(0., float(r.real)) for r in roots if abs(r.imag)<1e-8 and r.real>=-1e-10)
    unique = []
    for x in real:
        if not unique or abs(x-unique[-1])>1e-7:
            unique.append(x)
    return [{'x':x, 'slope':float(slope(x,a,n)), 'stable':bool(slope(x,a,n)<0)} for x in unique]


def folds(a):
    roots = np.roots([1, 0, 2, -2*a, 1])
    return sorted([{'x':float(z.real), 's':float(z.real-a*z.real**2/(1+z.real**2))}
                   for z in roots if abs(z.imag)<1e-8 and z.real>0], key=lambda z:z['s'])


def solve_segment(x0, duration, s, a, n=2, method='DOP853', rtol=1e-10, atol=1e-12, samples=None):
    sol = solve_ivp(lambda t,y: rhs(y,s,a,n), [0,duration], [x0], method=method,
                    rtol=rtol, atol=atol, t_eval=samples)
    assert sol.success, sol.message
    assert np.isfinite(sol.y).all() and np.min(sol.y)>-1e-9
    return sol.t, sol.y[0]


def sweep(a, n, dwell, levels):
    x=0.; up=[]; down=[]
    for s in levels:
        _, y=solve_segment(x,dwell,s,a,n)
        x=float(y[-1]);up.append(x)
    for s in levels[::-1]:
        _, y=solve_segment(x,dwell,s,a,n)
        x=float(y[-1]);down.append(x)
    return np.array(up), np.array(down[::-1])


def pulse(a,n,duration,s=.3,post=100.,method='DOP853',rtol=1e-10,atol=1e-12):
    t1,x1=solve_segment(0.,5.,0.,a,n,method,rtol,atol,np.linspace(0,5,101))
    t2,x2=solve_segment(x1[-1],duration,s,a,n,method,rtol,atol,np.linspace(0,duration,301))
    t3,x3=solve_segment(x2[-1],post,0.,a,n,method,rtol,atol,np.linspace(0,post,1001))
    return np.r_[t1,t2[1:]+5,t3[1:]+5+duration], np.r_[x1,x2[1:],x3[1:]], float(x2[-1])


def fixed_rk4(x0,duration,s,a,n,h):
    count=int(np.ceil(duration/h));h=duration/count;x=float(x0)
    for _ in range(count):
        f=lambda y:rhs(y,s,a,n)
        k1=f(x);k2=f(x+h*k1/2);k3=f(x+h*k2/2);k4=f(x+h*k3)
        x+=h*(k1+2*k2+2*k3+k4)/6
    return x


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    arrays={}; rows=[]
    reversible_folds=folds(1.8)
    assert len(reversible_folds)==2 and all(z['s']>0 for z in reversible_folds)
    low,high=[z['s'] for z in reversible_folds]
    def width(s):
        roots=equilibria(s,1.8)
        return roots[-1]['x']-roots[0]['x']
    static_area=quad(width,low,high,epsabs=1e-11,epsrel=1e-10)[0]
    levels=np.linspace(0,.3,601)
    threshold=np.mean([z['x'] for z in reversible_folds])
    for a,n,dwells in [(1.8,2,[5.,25.,125.]),(0.,2,[.5,2.,10.]),(1.4,2,[5.,25.,125.])]:
        for dwell in dwells:
            up,down=sweep(a,n,dwell,levels)
            key=f'a{a:g}_n{n}_d{dwell:g}'
            arrays[key+'_s'],arrays[key+'_up'],arrays[key+'_down']=levels,up,down
            row={'a':a,'n':n,'dwell_per_stimulus':dwell,
                 'loop_area':float(np.trapezoid(abs(down-up),levels)),
                 'max_up_down_difference':float(np.max(abs(down-up)))}
            if a==1.8:
                on=int(np.flatnonzero(up>threshold)[0])
                off=int(np.flatnonzero(down>threshold)[0])
                row.update(on_grid_bracket=[float(levels[on-1]),float(levels[on])],
                           off_grid_bracket=[float(levels[off-1]),float(levels[off])])
            if a==0:
                exact_up=[];exact_down=[];x=0
                for s in levels:
                    x=s+(x-s)*np.exp(-dwell);exact_up.append(x)
                for s in levels[::-1]:
                    x=s+(x-s)*np.exp(-dwell);exact_down.append(x)
                row['max_error_vs_linear_exact']=float(max(np.max(abs(up-exact_up)),np.max(abs(down-exact_down[::-1]))))
                assert row['max_error_vs_linear_exact']<1e-9
            rows.append(row)
            print(f'sweep a={a:g}, dwell={dwell:g}: area={row["loop_area"]:.8g}',flush=True)
    # Two stable states at the exact same stimulus, independent of sweep speed.
    interior=(low+high)/2
    basin_roots=equilibria(interior,1.8)
    assert len(basin_roots)==3 and [r['stable'] for r in basin_roots]==[True,False,True]
    basin_finals=[]
    for start in [0.,2.]:
        t,y=solve_segment(start,300,interior,1.8,samples=np.linspace(0,300,1201))
        arrays[f'basin_{start:g}_time'],arrays[f'basin_{start:g}_x']=t,y
        basin_finals.append(float(y[-1]))
    assert abs(basin_finals[0]-basin_roots[0]['x'])<1e-8
    assert abs(basin_finals[1]-basin_roots[-1]['x'])<1e-8

    # Persistent activation after a pulse, for a=3 at zero final input.
    separator=(3-np.sqrt(5))/2
    high_state=(3+np.sqrt(5))/2
    critical=quad(lambda x:1/rhs(x,.3,3.,2),0,separator,epsabs=1e-12,epsrel=1e-12)[0]
    numerical_critical=brentq(lambda duration:solve_segment(0,duration,.3,3.)[1][-1]-separator,
                             .5*critical,1.5*critical,xtol=1e-12)
    assert abs(critical-numerical_critical)<1e-8
    pulse_rows=[]
    for name,a,n,duration in [('short',3.,2,.95*critical),('long',3.,2,1.05*critical),
                              ('no_feedback',0.,2,1.05*critical),('weak_feedback',1.4,2,1.05*critical),
                              ('noncooperative',3.,1,1.05*critical)]:
        t,y,endpulse=pulse(a,n,duration)
        arrays[name+'_time'],arrays[name+'_x']=t,y
        pulse_rows.append({'name':name,'a':a,'n':n,'pulse_s':.3,'pulse_start':5.,'pulse_duration':duration,
                           'x_at_pulse_end':endpulse,'x_after_100_without_stimulus':float(y[-1])})
    assert pulse_rows[0]['x_after_100_without_stimulus']<1e-9
    assert abs(pulse_rows[1]['x_after_100_without_stimulus']-high_state)<1e-8
    assert abs(pulse_rows[2]['x_after_100_without_stimulus'])<1e-9
    assert abs(pulse_rows[3]['x_after_100_without_stimulus'])<1e-9
    assert abs(pulse_rows[4]['x_after_100_without_stimulus']-2)<1e-8
    # Persistence alone is insufficient evidence of bistability: n=1,a=3 has unstable zero.
    _,n1_seed=solve_segment(1e-6,100,0,3,1)
    _,n2_seed=solve_segment(1e-6,100,0,3,2)
    _,reset=solve_segment(high_state,30,0,0,2)
    assert abs(n1_seed[-1]-2)<1e-8 and abs(n2_seed[-1])<1e-9 and abs(reset[-1])<1e-9

    convergence=[]
    for factor in [.95,1.05]:
        duration=factor*critical
        t,reference,_=pulse(3,2,duration)
        _,independent,_=pulse(3,2,duration,method='Radau',rtol=1e-11,atol=1e-13)
        # Independent quadrature reference avoids adaptive-solver error masking
        # the very small RK4 differences. Report errors, not a fitted order.
        end=brentq(lambda x:quad(lambda y:1/rhs(y,.3,3.),0,x,epsabs=1e-13,epsrel=1e-13)[0]-duration,
                   0.,1.,xtol=1e-14)
        error=[]
        for h in [.05,.025]:
            error.append(abs(fixed_rk4(0,duration,.3,3.,2,h)-end))
        assert np.max(abs(reference-independent))<1e-7 and error[0]<1e-8 and error[1]<error[0]/2
        convergence.append({'duration_factor':factor,'max_DOP853_Radau_difference':float(np.max(abs(reference-independent))),
                            'RK4_endpoint_error_h005':error[0],'RK4_endpoint_error_h0025':error[1]})
    protocol={'source_equation':'dAp/dt = kp*S*A + beta*Ap^n/(K^n+Ap^n) - kd*Ap',
              'assumption':'A is constant/clamped; no total A+Ap conservation is specified or imposed.',
              'scaled_equation':'dx/dtau = s + a*x^n/(1+x^n) - x',
              'scaling':{'x':'Ap/K','tau':'kd*t','s':'kp*S*A/(kd*K)','a':'beta/(kd*K)'},
              'parameter_choice':'Illustrative dimensionless values, not measured biological parameters.',
              'solver':'solve_ivp DOP853, rtol=1e-10, atol=1e-12; integrate discontinuous pulse segments separately.',
              'sweep':'601 levels from s=0 to 0.3, then reverse, continuing from the previous final state. Each level held for a specified dwell. Up/down arrays aligned at equal s.',
              'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
              'randomness':'None. All calculations deterministic; exact final digits may vary by platform.',
              'scope':'Phenomenological concentration ODE; neither water MD nor an individual-molecule hierarchy model.'}
    result={'critical_feedback_a_for_n2':8/(3*np.sqrt(3)),
            'reversible':{'a':1.8,'n':2,'folds':reversible_folds,'static_loop_area':static_area,
                          'basin_stimulus':interior,'basin_roots':basin_roots,'basin_finals':basin_finals},
            'persistent':{'a':3,'n':2,'folds':folds(3),'equilibria_at_s0':equilibria(0,3),
                          'separator':separator,'high_state':high_state,'critical_pulse_duration_at_s03':critical,
                          'critical_duration_ODE_check':numerical_critical,'pulse_runs':pulse_rows},
            'sweeps':rows,'convergence':convergence,
            'controls':{'n1_a3_zero_slope':slope(0,3,1),'n1_tiny_seed_final':float(n1_seed[-1]),
                        'n2_tiny_seed_final':float(n2_seed[-1]),'feedback_removed_reset_after_30':float(reset[-1])},
            'checks_passed':True}
    np.savez_compressed(out/'trajectories.npz',**arrays)
    for filename,data in [('protocol.json',protocol),('results.json',result)]:
        (out/filename).write_text(json.dumps(data,indent=2),encoding='utf-8')
    print(json.dumps({'folds':reversible_folds,'critical_duration':critical,'pulse_runs':pulse_rows,'passed':True},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=HERE/'data')
    args=parser.parse_args();run(args.out)
