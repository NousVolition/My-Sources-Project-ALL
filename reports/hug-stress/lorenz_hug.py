"""Optional one-way Lorenz pressure drive for the reduced hug.

This is an explicit proposed coupling, not a derived fluid/material law.
Lorenz rho is distinct from the hug's stiffness r. All units are model units.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from stress_solver import model

@dataclass(frozen=True)
class Drive:
    rho:float=28.
    sigma:float=10.
    beta:float=8/3
    mean:float=.45
    amplitude:float=.35
    scale:float=10.
    def __post_init__(self):
        if not all(np.isfinite(x) for x in asdict(self).values()):raise ValueError('Finite drive parameters required')
        if min(self.rho,self.sigma,self.beta,self.scale)<=0 or self.mean<self.amplitude or self.amplitude<0:raise ValueError('Invalid drive parameters')

def lorenz(s,d):
    x,y,z=s
    return np.array([d.sigma*(y-x),x*(d.rho-z)-y,x*y-d.beta*z])

def lorenz_jac(s,d):
    x,y,z=s
    return np.array([[-d.sigma,d.sigma,0],[d.rho-z,-1,-x],[y,x,-d.beta]])

def pressure(z,d):
    return d.mean+d.amplitude*np.tanh((z-max(d.rho-1,0))/d.scale)

def rhs(t,y,d,p):
    q,v,m,_,_=y[3:];load=float(pressure(y[2],d));tau=.8 if load>m else 6.
    feedback=p.memory_coupling*m*q
    return np.r_[lorenz(y[:3],d),v,p.r*q-q**3-p.damping*v+feedback,(load-m)/tau,feedback*v,p.damping*v*v]

def warm_start(d,initial=(1.,1.,1.)):
    sol=solve_ivp(lambda t,y:lorenz(y,d),(0,50),initial,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.025)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y[:,-1]

def integrate(d,p,end=40.,initial=None,q0=.04,method='DOP853',rtol=1e-10,max_step=.02):
    start=warm_start(d) if initial is None else np.asarray(initial)
    y0=np.r_[start,q0,0,0,0,0]
    def opening(t,y):return float(pressure(y[2],d))-1
    opening.direction=1;opening.terminal=False
    sol=solve_ivp(lambda t,y:rhs(t,y,d,p),(0,end),y0,method=method,rtol=rtol,atol=rtol*.01,
                  max_step=max_step,dense_output=True,events=opening)
    if not sol.success or not np.isfinite(sol.y).all():raise RuntimeError(sol.message)
    opened=0. if pressure(start[2],d)>=1 else (float(sol.t_events[0][0]) if len(sol.t_events[0]) else None)
    return sol,opened,start

def geometry(t,state,p,opening,points=65):
    q,_,m=state[3:6];sx=1+m*(.04+.02*np.sin(2*np.pi*t/3));sy=1/sx
    bend=1. if opening is None else 1-model.smoothstep((t-opening)/.8)
    left,right=model.arm_points(bend,sx,sy,points)
    for arm in (left,right):
        for pt in arm:pt[0]+=p.shear*q*pt[1]
    return dict(left=left,right=right,sx=sx,sy=sy,bend=bend)

def lyapunov(d,initial=(1.,1.,1.),rtol=1e-9,max_step=.03,duration=300.,burn=50.):
    # Evolve the tangent equation, renormalizing every quarter model-time unit.
    y=np.r_[initial,[1.,0,0]];logs=[];times=[]
    for i in range(int(round((duration+burn)/.25))):
        def f(t,s):return np.r_[lorenz(s[:3],d),lorenz_jac(s[:3],d)@s[3:]]
        sol=solve_ivp(f,(0,.25),y,method='DOP853',rtol=rtol,atol=rtol*.01,max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:,-1];length=np.linalg.norm(y[3:]);y[3:]/=length
        if (i+1)*.25>burn:logs.append(float(np.log(length)));times.append((i+1)*.25-burn)
    logs=np.array(logs);cumulative=np.cumsum(logs)/np.array(times)
    return dict(estimate=float(cumulative[-1]),duration=duration,burn=burn,rtol=rtol,max_step=max_step,
                initial=list(initial),block_estimates=[float(sum(x)/(len(x)*.25)) for x in np.array_split(logs,6)]),np.array(times),cumulative

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rho',type=float,default=28);parser.add_argument('--end',type=float,default=40)
    parser.add_argument('--output',type=Path,default=Path('lorenz-hug-run.npz'))
    args=parser.parse_args();d=Drive(rho=args.rho);p=model.Parameters(r=-.25,memory_coupling=1.5)
    sol,opened,start=integrate(d,p,end=args.end);t=np.linspace(0,args.end,int(args.end/.01)+1)
    np.savez_compressed(args.output,time=t,state=sol.sol(t).T,pressure=pressure(sol.sol(t)[2],d))
    print(json.dumps(dict(drive=asdict(d),hug=asdict(p),opening=opened,warm_start=start.tolist(),output=str(args.output))))
