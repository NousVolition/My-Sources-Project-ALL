"""Dimensionless, noiseless, overdamped RSJ pair with a resistive load.

M phi' = i - sin(phi), M = I + alpha * ones((2,2)).
tau = (2e Ic r / hbar) t; V_total/(Ic r) = phi1' + phi2'.
This implements the same circuit as fluid_pilot/dynamics.py:circuit_rhs.
All phases remain unwrapped. Integrals use the RK4 stage quadrature.
"""
import numpy as np


def rhs(phi, bias, alpha):
    phi = np.asarray(phi, dtype=float)
    b = np.asarray(bias)[..., None] - np.sin(phi)
    a = np.asarray(alpha)[..., None]
    return b - a / (1 + 2 * a) * b.sum(axis=-1, keepdims=True)


def potential(phi):
    return -np.cos(phi).sum(axis=-1)


def dissipation(v, alpha):
    return (v * v).sum(axis=-1) + np.asarray(alpha) * v.sum(axis=-1)**2


def step(phi, bias, alpha, dt):
    k1 = rhs(phi, bias, alpha)
    k2 = rhs(phi + dt/2*k1, bias, alpha)
    k3 = rhs(phi + dt/2*k2, bias, alpha)
    k4 = rhs(phi + dt*k3, bias, alpha)
    p = phi + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    d = dt/6*(dissipation(k1, alpha) + 2*dissipation(k2, alpha)
              + 2*dissipation(k3, alpha) + dissipation(k4, alpha))
    return p, d


def integrate(phi0, bias, alpha, duration, dt, stride=1.0):
    """Batch independent trajectories; return phase and dissipated-energy samples."""
    if dt <= 0 or duration <= 0 or stride <= 0 or np.any(np.asarray(alpha) < 0):
        raise ValueError('Require positive times and nonnegative resistance ratio')
    n, every = round(duration/dt), round(stride/dt)
    if every < 1 or not np.isclose(n*dt, duration) or not np.isclose(every*dt, stride) or n % every:
        raise ValueError('duration and stride must be integral multiples of dt, and duration of stride')
    p = np.asarray(phi0, dtype=float).copy()
    if p.shape[-1] != 2 or not np.isfinite(p).all():
        raise ValueError('Expected finite phase pairs')
    acc = np.zeros(p.shape[:-1]); samples = [p.copy()]; energies = [acc.copy()]
    for j in range(n):
        p, d = step(p, bias, alpha, dt)
        acc += d
        if (j + 1) % every == 0:
            samples.append(p.copy()); energies.append(acc.copy())
    phases, dissipated = np.asarray(samples), np.asarray(energies)
    if not np.isfinite(phases).all():
        raise FloatingPointError('Nonfinite integration')
    return np.arange(len(samples))*stride, phases, dissipated


def voltage(start, end, window):
    if window <= 0:
        raise ValueError('Voltage window must be positive')
    return (np.asarray(end) - np.asarray(start)) / window


def power_residual(phases, dissipated, bias):
    work = np.asarray(bias)*(phases[-1] - phases[0]).sum(axis=-1)
    du = potential(phases[-1]) - potential(phases[0])
    residual = work - du - dissipated[-1]
    scale = 1 + abs(work) + abs(du) + dissipated[-1]
    return residual / scale


def synchronized_voltage(bias, alpha):
    """Infinite-time total voltage on the invariant synchronized branch only."""
    b = np.asarray(bias)
    return 2*np.sign(b)*np.sqrt(np.maximum(b*b - 1, 0))/(1 + 2*np.asarray(alpha))


def to_state(phi):
    """Return (p, log|tan(q/2)|), sign(q), with p=mean(phi), q=half difference.

    Initial phase differences must lie within (-2pi,2pi). These invariant
    strips suffice for the declared initial draws. Keeping log separation
    avoids artificial synchronization when q becomes smaller than an ulp of p.
    """
    phi=np.asarray(phi,float);p=phi.mean(-1);q=(phi[...,0]-phi[...,1])/2
    if np.any(abs(q)>=np.pi):raise ValueError('Initial half difference must be within (-pi,pi)')
    with np.errstate(divide='ignore'):
        z=np.log(abs(np.tan(q/2)))
    return np.stack([p,z],axis=-1),np.where(q<0,-1.,1.)


def from_state(state, sign):
    p,z=np.moveaxis(np.asarray(state),-1,0)
    q=2*np.asarray(sign)*np.arctan2(np.exp(np.minimum(z,0)),np.exp(np.minimum(-z,0)))
    return np.stack([p+q,p-q],axis=-1)


def stable_rhs(state,bias,alpha):
    p,z=np.moveaxis(np.asarray(state),-1,0)
    dp=(np.asarray(bias)+np.sin(p)*np.tanh(z))/(1+2*np.asarray(alpha))
    return np.stack([dp,-np.cos(p)],axis=-1)


def stable_velocity(state,sign,bias,alpha):
    p,z=np.moveaxis(np.asarray(state),-1,0)
    dp=(np.asarray(bias)+np.sin(p)*np.tanh(z))/(1+2*np.asarray(alpha))
    ez=np.exp(-abs(z));dq=-np.cos(p)*np.asarray(sign)*2*ez/(1+ez*ez)
    return np.stack([dp+dq,dp-dq],axis=-1)


def integrate_state(state0,sign,bias,alpha,duration,dt,stride=5.):
    """RK4 in analytically equivalent mean/log-separation coordinates."""
    n=round(duration/dt);every=round(stride/dt)
    if dt<=0 or every<1 or not np.isclose(n*dt,duration) or not np.isclose(every*dt,stride) or n%every:
        raise ValueError('Invalid time grid')
    state=np.asarray(state0,float).copy();acc=np.zeros(state.shape[:-1])
    samples=[state.copy()];energies=[acc.copy()]
    for j in range(n):
        k1=stable_rhs(state,bias,alpha);s2=state+dt/2*k1
        k2=stable_rhs(s2,bias,alpha);s3=state+dt/2*k2
        k3=stable_rhs(s3,bias,alpha);s4=state+dt*k3
        k4=stable_rhs(s4,bias,alpha)
        acc+=dt/6*(dissipation(stable_velocity(state,sign,bias,alpha),alpha)
                   +2*dissipation(stable_velocity(s2,sign,bias,alpha),alpha)
                   +2*dissipation(stable_velocity(s3,sign,bias,alpha),alpha)
                   +dissipation(stable_velocity(s4,sign,bias,alpha),alpha))
        state+=dt/6*(k1+2*k2+2*k3+k4)
        if (j+1)%every==0:samples.append(state.copy());energies.append(acc.copy())
    states=np.array(samples)
    if not np.isfinite(states[...,0]).all() or np.isnan(states[...,1]).any():
        raise FloatingPointError('Invalid transformed trajectory')
    return np.arange(len(samples))*stride,states,np.array(energies)
