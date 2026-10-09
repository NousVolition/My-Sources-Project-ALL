"""Uniform lapping, nonuniform rotations, phase locking, and overdamped limits.

Supplied equations reduce to phi'=mu-sin(phi). Rotations live on a circle;
unwrapped phase is used for detecting slips and measuring complete periods.
"""
from pathlib import Path
import argparse, json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent


def rotation_period(omega,a):
    if not 0 <= a < omega:
        raise ValueError('Positive circulating case requires 0 <= a < omega.')
    exact=2*np.pi/np.sqrt(omega**2-a**2)
    quadrature=quad(lambda x:1/(omega-a*np.sin(x)),0,2*np.pi,
                    points=[np.pi/2],epsabs=1e-10,epsrel=1e-11)[0]
    def event(t,y):return y[0]-2*np.pi
    event.terminal=True;event.direction=1
    sol=solve_ivp(lambda t,y:omega-a*np.sin(y),[0,1.2*exact],[0.],method='DOP853',
                  rtol=1e-11,atol=1e-13,max_step=.1,events=event,dense_output=True)
    assert sol.success and len(sol.t_events[0])==1
    numeric=float(sol.t_events[0][0])
    t=np.linspace(0,numeric,1001)
    bottleneck=quad(lambda x:1/(omega-a*np.sin(x)),np.pi/2-.2,np.pi/2+.2,
                    points=[np.pi/2],epsabs=1e-10,epsrel=1e-11)[0]
    assert abs(numeric/exact-1)<1e-7 and abs(quadrature/exact-1)<1e-9
    return {'omega':omega,'a':a,'exact_period':exact,'ODE_period':numeric,'quadrature_period':quadrature,
            'ODE_relative_error':abs(numeric/exact-1),'fraction_time_in_04_rad_bottleneck':bottleneck/exact},t,sol.sol(t)[0]


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    arrays={};periods=[]
    for a in [0.,.5,.9,.99,.999,.9999]:
        record,t,x=rotation_period(1.,a)
        periods.append(record)
        arrays[f'rotation_a{a:g}_time'],arrays[f'rotation_a{a:g}_theta']=t,x
    delta=np.logspace(-7,-3,100)
    exact=2*np.pi/np.sqrt(1-(1-delta)**2)
    exponent=float(np.polyfit(np.log(delta),np.log(exact),1)[0])
    assert abs(exponent+.5)<.001
    arrays['scaling_delta'],arrays['scaling_period']=delta,exact

    locking=[]
    for mu in [-1.1,-.5,0.,.5,.99,1.01,1.1,1.5]:
        t=np.linspace(0,300,6001)
        sol=solve_ivp(lambda t,y:mu-np.sin(y),[0,300],[-.4],method='DOP853',
                      t_eval=t,rtol=1e-10,atol=1e-12,max_step=.2,dense_output=True)
        assert sol.success
        phi=sol.y[0];mean_rate=(phi[-1]-phi[4000])/(t[-1]-t[4000])
        record={'mu':mu,'final_unwrapped_phi':float(phi[-1]),'late_mean_phase_drift':float(mean_rate)}
        if abs(mu)<1:
            target=np.arcsin(mu)
            wrapped_error=float(abs(np.arctan2(np.sin(phi[-1]-target),np.cos(phi[-1]-target))))
            record.update(regime='locked',stable_phase=float(target),phase_error=wrapped_error,
                          stable_slope=float(-np.cos(target)))
            assert wrapped_error<1e-7 and abs(mean_rate)<1e-7
        else:
            slip=2*np.pi/np.sqrt(mu**2-1)
            record.update(regime='slipping',analytic_slip_period=float(slip),analytic_mean_drift=float(np.sign(mu)*np.sqrt(mu**2-1)))
            lower,upper=sorted([phi[0],phi[-1]])
            turns=np.arange(np.ceil(lower/(2*np.pi)),np.floor(upper/(2*np.pi))+1)*2*np.pi
            crossings=np.sort([brentq(lambda tt:sol.sol(tt)[0]-target,0,300,xtol=1e-11) for target in turns])
            intervals=np.diff(crossings)
            assert len(intervals)>1 and np.max(abs(intervals/slip-1))<1e-7
            record.update(measured_slip_period=float(np.mean(intervals)),
                          largest_full_cycle_relative_error=float(np.max(abs(intervals/slip-1))),
                          complete_cycle_count=len(intervals))
            arrays[f'mu{mu:g}_crossing_times']=crossings
        arrays[f'mu{mu:g}_time'],arrays[f'mu{mu:g}_phi']=t,phi
        locking.append(record)
    # Critical case: approach the degenerate equilibrium, no finite full turn.
    critical=solve_ivp(lambda t,y:1-np.sin(y),[0,300],[0.],rtol=1e-11,atol=1e-13,method='DOP853',max_step=.2)
    assert critical.success and critical.y[0,-1]<np.pi/2

    # First-order / full inertial pendulum at matched on-step initial conditions.
    pendulum=[]
    for mu in [.8,1.2]:
        t=np.linspace(0,40,2001)
        first=solve_ivp(lambda t,y:mu-np.sin(y),[0,40],[0.],t_eval=t,method='DOP853',rtol=1e-11,atol=1e-13).y[0]
        arrays[f'pendulum_mu{mu:g}_time'],arrays[f'pendulum_mu{mu:g}_first']=t,first
        errors=[]
        for eps in [.1,.01,.001]:
            sol=solve_ivp(lambda t,y:[y[1],(mu-np.sin(y[0])-y[1])/eps],
                          [0,40],[0.,0.],t_eval=t,method='Radau',rtol=1e-10,atol=1e-12)
            assert sol.success
            error=float(np.max(abs(sol.y[0]-first)))
            errors.append(error)
            arrays[f'pendulum_mu{mu:g}_eps{eps:g}_phi']=sol.y[0]
            pendulum.append({'mu':mu,'epsilon':eps,'max_unwrapped_angle_difference':error,
                             'final_full_angle':float(sol.y[0,-1]),'final_first_order_angle':float(first[-1])})
        assert errors[1]<errors[0]/5 and errors[2]<errors[1]/5

    T1,T2=60.,75.
    w1,w2=2*np.pi/T1,2*np.pi/T2
    lap=2*np.pi/(w1-w2)
    assert abs(lap-T1*T2/(T2-T1))<1e-10
    results={'lapping':{'T1_seconds':T1,'T2_seconds':T2,'first_lap_seconds':lap,
                        'relative_phase_at_lap':float((w1-w2)*lap)},
             'periods':periods,'near_threshold_log_log_exponent':exponent,
             'critical_mu1_final_angle_after_300':float(critical.y[0,-1]),
             'locking':locking,'pendulum':pendulum,'checks_passed':True}
    protocol={'uniform':'phi_dot=2*pi/T1-2*pi/T2; Tlap=T1*T2/(T2-T1) for T2>T1.',
              'nonuniform':'theta_dot=omega-a*sin(theta); positive rotation requires omega>a>=0.',
              'period':'T=2*pi/sqrt(omega^2-a^2). Near a=omega from below: T~pi*sqrt(2/omega)*(omega-a)^(-1/2).',
              'firefly':'Theta_dot=Omega; theta_dot=omega+A*sin(Theta-theta); phi=Theta-theta; tau=A*t; mu=(Omega-omega)/A; phi_prime=mu-sin(phi).',
              'locking':'Strict |mu|<1: stable phase arcsin(mu) modulo 2*pi. At |mu|=1 the equilibrium is degenerate; |mu|>1 gives persistent phase slips.',
              'boundary':'Exact unstable initial equilibria are exceptions to attraction. Locked phases need not coincide.',
              'pendulum':'I*theta_tt+b*theta_t+mgL*sin(theta)=Gamma, I=mL^2. tau=(mgL/b)*t, mu=Gamma/(mgL), epsilon=I*mgL/b^2=m^2*g*L^3/b^2.',
              'pendulum_initial_conditions':'theta=0, physical angular velocity=0; full and reduced angles compared at identical scaled times. Initial velocity layer is expected.',
              'scope':'Illustrative deterministic phase models; no firefly data fitted, no inference of water synchronization.',
              'numerics':'DOP853 for first-order dynamics; Radau for stiff inertial pendulum. Unwrapped phases used; independent quadrature validates periods.'}
    np.savez_compressed(out/'phase_trajectories.npz',**arrays)
    (out/'phase_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    (out/'phase_protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    print(json.dumps({'lapping':results['lapping'],'periods':periods,'scaling_exponent':exponent,'pendulum':pendulum,'passed':True},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=HERE/'data')
    args=parser.parse_args();run(args.out)
