"""Exploratory two-way closure; this is not a calibrated fluid interaction law.

State: X,Y,Z,q,v,imprint,hug_work,hug_loss,feedback_loss,Lorenz_native_work.
The new term removes Lorenz quadratic activity; it is not a conservative
transfer of physical energy to the hug. The original one-way code is retained.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from lorenz_hug import Drive, pressure, rhs as one_way_rhs, model


def bend(t, opening):
    return 1. if opening is None else 1.-model.smoothstep((t-opening)/.8)


def resistance(t, imprint, strength, opening):
    m=max(float(imprint),0.)
    return strength*bend(t,opening)**2*m/(1.+m)


def rhs(t,y,d,p,strength,opening):
    base=one_way_rhs(t,y[:8],d,p)
    x,z_y,z=y[:3]
    eta=resistance(t,y[5],strength,opening)
    base[:2]-=eta*y[:2]
    loss=eta*(x*x+z_y*z_y)
    native=(d.sigma+d.rho)*x*z_y-d.sigma*x*x-z_y*z_y-d.beta*z*z
    return np.r_[base,loss,native]


@dataclass
class Trajectory:
    pieces:list
    opening:float|None
    nfev:int
    def sol(self,t):
        times=np.atleast_1d(t).astype(float)
        result=np.empty((10,len(times)))
        for i,piece in enumerate(self.pieces):
            mask=(times>=piece.t[0])&(times<=piece.t[-1])
            result[:,mask]=piece.sol(times[mask])
        if np.any(times<self.pieces[0].t[0]) or np.any(times>self.pieces[-1].t[-1]):
            raise ValueError('Requested time outside completed integration')
        return result[:,0] if np.ndim(t)==0 else result


def integrate(d,p,strength,end=80.,initial=(1.,1.,1.),q0=.04,rtol=1e-10,max_step=.02,method='DOP853'):
    if not np.isfinite(strength) or strength<0:raise ValueError('Nonnegative finite feedback strength required')
    if end<=0:raise ValueError('Positive duration required')
    y0=np.r_[np.asarray(initial,dtype=float),q0,0.,0.,0.,0.,0.,0.]
    if y0.shape!=(10,) or not np.isfinite(y0).all():raise ValueError('Three finite Lorenz initial coordinates required')
    opened=0. if pressure(y0[2],d)>=1 else None
    def event(t,y):return float(pressure(y[2],d))-1
    event.direction=1;event.terminal=True
    options=dict(method=method,rtol=rtol,atol=rtol*.01,max_step=max_step,dense_output=True)
    first=solve_ivp(lambda t,y:rhs(t,y,d,p,strength,opened),(0,end),y0,
                    events=event if opened is None else None,**options)
    pieces=[first]
    if opened is None and len(first.t_events[0]):
        opened=float(first.t_events[0][0])
        # Stop at the first opening, then resume with an irreversible smooth release.
        if opened<end:
            pieces.append(solve_ivp(lambda t,y:rhs(t,y,d,p,strength,opened),
                                   (opened,end),first.y[:,-1],**options))
    for piece in pieces:
        if not piece.success or not np.isfinite(piece.y).all():raise RuntimeError(piece.message)
    return Trajectory(pieces,opened,sum(s.nfev for s in pieces))
