"""Independent stress drivers for the unchanged, pinned reduced hug equations."""
from pathlib import Path
import sys
import math
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0,str(Path(__file__).resolve().parent/'baseline'))
import hug_model as model


def fixed(p,q0,v0,end,dt,resolve_pulse=False):
    """Same RK4 step; optional explicit subdivision of the pressure interval.

    A safety guard ends a numerical trial, not the mathematical solution.
    Only accepted finite states are retained. No clipping changes the ODE.
    """
    if resolve_pulse:
        cuts=[0,min(p.pulse_duration,end),end]
        intervals=[]
        for a,b in zip(cuts,cuts[1:]):
            if b<=a:continue
            h=min(dt,p.pulse_duration/64) if a==0 else dt
            n=int(math.ceil((b-a)/h));intervals.extend(np.linspace(a,b,n+1)[1:])
        target=np.array(intervals)
    else:
        n=int(round(end/dt));target=np.linspace(0,end,n+1)[1:]
    y=np.array([q0,v0,0.,0.,0.]);times=[0.];states=[y.copy()]
    failure=None
    for next_time in target:
        t=times[-1];h=next_time-t
        try:
            with np.errstate(over='raise',invalid='raise',divide='raise'):
                yn=model.step(lambda tt,yy:model.rhs(tt,yy,p),t,y,h,'rk4')
            if not np.isfinite(yn).all():raise FloatingPointError('nonfinite state')
            if np.max(np.abs(yn))>1e12:
                failure=dict(kind='guard_exceeded',time=float(next_time),largest_attempt=float(np.max(abs(yn))));break
        except (FloatingPointError,OverflowError) as error:
            failure=dict(kind='floating_point',time=float(next_time),message=str(error));break
        y=yn;times.append(float(next_time));states.append(y.copy())
    return np.array(times),np.array(states),failure


def jacobian(t,y,p):
    q,v,m,_,_=y;load=model.pressure(t,p);tau=.8 if load>m else 6.
    J=np.zeros((5,5));J[0,1]=1
    J[1,:3]=[p.r-3*q*q+p.memory_coupling*m,-p.damping,p.memory_coupling*q]
    J[2,2]=-1/tau
    J[3,:3]=[p.memory_coupling*m*v,p.memory_coupling*m*q,p.memory_coupling*q*v]
    J[4,1]=2*p.damping*v
    return J


def resolved(p,q0,v0,end,method='DOP853',rtol=2e-10):
    """Always resolve the prescribed pressure pulse before continuing.

    The solver is restarted at the known end of the pulse. Tight tolerances
    alone cannot require a solver to notice an unsampled forcing interval.
    """
    spans=[(0.,min(p.pulse_duration,end))]
    if p.pulse_duration<end:spans.append((p.pulse_duration,end))
    y=np.array([q0,v0,0.,0.,0.]);solutions=[];nfev=0
    for i,(a,b) in enumerate(spans):
        options={'jac':lambda t,z:jacobian(t,z,p)} if method=='Radau' else {}
        sol=solve_ivp(lambda t,z:model.rhs(t,z,p),(a,b),y,method=method,
                      rtol=rtol,atol=rtol*.01,max_step=min(.05,p.pulse_duration/32) if i==0 else .2,
                      dense_output=True,**options)
        if not sol.success or not np.isfinite(sol.y).all():raise RuntimeError(sol.message)
        y=sol.y[:,-1];solutions.append(sol);nfev+=sol.nfev
    def evaluate(times):
        times=np.asarray(times);out=np.empty((len(times),5))
        for i,sol in enumerate(solutions):
            mask=(times>=sol.t[0])&(times<=sol.t[-1])
            out[mask]=sol.sol(times[mask]).T
        return out
    return evaluate,dict(method=method,rtol=rtol,atol=rtol*.01,nfev=nfev,
                         accepted_steps=sum(len(s.t)-1 for s in solutions))


def observation_grid(p,end):
    # Resolve pulse and sample the full interval; peaks are still sampled peaks.
    return np.unique(np.r_[np.linspace(0,end,max(1201,int(end/.05)+1)),
                           np.linspace(0,min(p.pulse_duration,end),257)])


def scaled_error(values,reference):
    scale=np.maximum(1,np.max(abs(reference[:,:3]),axis=0))
    return float(np.max(abs(values[:,:3]-reference[:,:3])/scale))


def energy_budget(y,p):
    e=model.energy(y,p);residual=e-e[0]-y[:,3]+y[:,4]
    scale=max(1.,float(max(abs(e))),float(max(abs(y[:,3]))),float(max(abs(y[:,4]))))
    return dict(absolute=float(max(abs(residual))),scaled=float(max(abs(residual)))/scale)


def extent(times,values,p):
    """Exact arm extent at each sampled time, not a coarse vertex estimate."""
    sx=1+values[:,2]*(.04+.02*np.sin(2*np.pi*times/3));sy=1/sx
    opening=model.opening_time(p)
    if opening is None:bend=np.ones(len(times))
    else:bend=np.array([1-model.smoothstep((t-opening)/.8) for t in times])
    c=p.shear*values[:,0]*sy
    xmax=sx*(1-bend)+np.hypot(sx*bend,c)
    overall=np.maximum(xmax,sy)
    return overall,dict(maximum=float(max(overall)),outside_cube=bool(max(overall)>3),
        first_observed_outside=None if max(overall)<=3 else float(times[np.flatnonzero(overall>3)[0]]),
        opening_time=opening)
