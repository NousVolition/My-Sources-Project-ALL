"""Abstract SIMS and continuously evolving directed relationship guards."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from scipy.special import expit

def make_cases(p):
    rng=np.random.default_rng(p['seed'])
    starts=[]
    for seed in range(p['base_pairs']):
        b=np.array([rng.uniform(*p['baseline_sampling']['A']),rng.uniform(*p['baseline_sampling']['B'])])
        x=b+rng.uniform(*p['baseline_sampling']['initial_offset'],size=2)
        starts.append((seed,b,x))
    pairs=[]
    for swapped in (False,True):
        for seed,b,x in starts:
            pairs.append(dict(seed=seed,swapped=swapped,b=(b[::-1] if swapped else b).tolist(),
                              x0=(x[::-1] if swapped else x).tolist()))
    cases=[]
    for pair_id,pair in enumerate(pairs):
        for mode in p['guard_modes']:
            for weight,k in p['weights'].items():
                for disturbance in p['disturbances']:
                    cases.append(dict(id=len(cases),pair=pair_id,**pair,mode=mode,weight=weight,k=k,
                                      disturbance=disturbance))
    return cases

def arrays(cases,p):
    b=np.array([c['b'] for c in cases]);k=np.array([c['k'] for c in cases])
    mode=np.array([list(p['guard_modes']).index(c['mode']) for c in cases])
    force=np.zeros_like(b)
    for i,c in enumerate(cases):
        if c['disturbance']!='none':
            force[i,0 if c['disturbance'].startswith('A') else 1]=1 if c['disturbance'].endswith('positive') else -1
    y=np.column_stack([np.array([c['x0'] for c in cases]),np.tile(p['guard_initial'],(len(cases),1))])
    return b,k,mode,force,y

def rhs(t,y,b,k,mode,force,p):
    d=y[:,1]-y[:,0]
    target=np.column_stack([expit(p['theta']+p['beta']*d),expit(p['theta']-p['beta']*d)])
    shared=mode==0
    target[shared]=np.mean(target[shared],axis=1)[:,None]
    dg=(target-y[:,2:])/p['tau'];dg[mode==2]=0
    dx=p['lambda']*(b-y[:,:2])+k*y[:,2:]*np.column_stack([d,-d])
    q=(t-p['pulse']['onset'])/p['pulse']['duration']
    if 0<q<1:dx=dx+p['pulse']['amplitude']*np.sin(np.pi*q)**2*force
    return np.column_stack([dx,dg])

def integrate(p,cases,dt):
    import time
    b,k,mode,force,y=arrays(cases,p)
    n=round(p['time']['end']/dt);stride=round(p['time']['output_dt']/dt)
    assert n<=p['resources']['maximum_step_count'] and abs(stride*dt-p['time']['output_dt'])<1e-14
    hist=[y.copy()];lo=1.;hi=0.;stages=0;start=time.monotonic()
    def check(z):
        nonlocal lo,hi,stages
        stages+=1
        if not np.isfinite(z).all():raise ArithmeticError('Nonfinite stage')
        lo=min(lo,float(z[:,2:].min()));hi=max(hi,float(z[:,2:].max()))
        if lo < -p['verification']['guard_interval_tolerance'] or hi > 1+p['verification']['guard_interval_tolerance']:
            raise ArithmeticError('Guard escaped interval; clipping forbidden')
    for step in range(n):
        t=step*dt
        check(y);a=rhs(t,y,b,k,mode,force,p)
        z=y+dt*a/2;check(z);bb=rhs(t+dt/2,z,b,k,mode,force,p)
        z=y+dt*bb/2;check(z);c=rhs(t+dt/2,z,b,k,mode,force,p)
        z=y+dt*c;check(z);d=rhs(t+dt,z,b,k,mode,force,p)
        y=y+dt*(a+2*bb+2*c+d)/6
        check(y)
        if (step+1)%stride==0:hist.append(y.copy())
        if step%500==0 and time.monotonic()-start>p['resources']['max_elapsed_seconds']:
            raise TimeoutError('Per-integration wall budget')
    return np.array(hist),dict(dt=dt,steps=n,rhs_evaluations=4*n,checked_stages=stages,
                               guard_min=lo,guard_max=hi,elapsed_seconds=time.monotonic()-start,
                               clipping=False)

def independent_rhs(p,case):
    """Scalar implementation, intentionally independent of vectorized rhs."""
    import math
    a,b=case['b'];ka,kb=case['k'];mode=case['mode'];pulse=case['disturbance']
    def fun(t,y):
        xa,xb,ga,gb=y
        ta=1/(1+math.exp(-(p['theta']+p['beta']*(xb-xa))))
        tb=1/(1+math.exp(-(p['theta']+p['beta']*(xa-xb))))
        if mode=='shared':ta=tb=(ta+tb)/2
        dga=(ta-ga)/p['tau'];dgb=(tb-gb)/p['tau']
        if mode=='frozen':dga=dgb=0.
        da=p['lambda']*(a-xa)+ka*ga*(xb-xa)
        db=p['lambda']*(b-xb)+kb*gb*(xa-xb)
        phase=(t-p['pulse']['onset'])/p['pulse']['duration']
        if 0<phase<1 and pulse!='none':
            amount=p['pulse']['amplitude']*math.sin(math.pi*phase)**2
            if pulse.endswith('negative'):amount=-amount
            if pulse.startswith('A'):da+=amount
            else:db+=amount
        return [da,db,dga,dgb]
    return fun

