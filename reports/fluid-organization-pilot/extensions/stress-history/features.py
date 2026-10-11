"""Prediction inputs stop at the impulse; future arrays are used only for targets."""
import numpy as np
def inputs(z,lags):
    now=np.asarray(z['history'])[-1];t=np.asarray(z['history_time']);past=[]
    if len(t)!=len(z['history']) or not np.isclose(t[-1],0) or np.any(t>1e-12) or np.any(np.diff(t)<=0):raise ValueError('Inputs must be ordered past observations ending at intervention time')
    for lag in lags:
        idx=np.flatnonzero(np.isclose(t,-lag))
        if len(idx)!=1:raise ValueError('Missing causal lag')
        past.extend(z['history'][idx[0]]-now)
    strength=np.asarray(z['relative_amplitudes']);rates=np.asarray(z['initial_rates'])/np.asarray(z['initial_perturbation_energy'])[:,None]
    conventional=np.c_[np.tile(now[:8],(len(strength),1)),strength]
    present=np.c_[np.tile(now,(len(strength),1)),strength,rates]
    history=np.tile(np.asarray(past),(len(strength),1))
    return conventional,present,history
def target(z):
    return np.trapezoid(z['perturbation_energy']/z['initial_perturbation_energy'][None,:],z['time'],axis=0)
def shuffle_groups(h,permutation):
    h=np.asarray(h)
    if len(h)%2:raise ValueError('Two disturbance strengths per independent field')
    return h.reshape(-1,2,h.shape[-1])[permutation].reshape(h.shape)
def derange(n,seed):
    rng=np.random.default_rng(seed)
    for _ in range(10000):
        p=rng.permutation(n)
        if np.all(p!=np.arange(n)):return p
    raise RuntimeError('Cannot derange')
