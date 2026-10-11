"""Analytic, geometric, gradient and reciprocal-work checks before stress runs."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from clay_hug import *
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def main():
    DATA.mkdir(exist_ok=True);rows=[];rng=np.random.default_rng(20261010)
    p=Parameters(alpha=.5);mesh=Mesh(p)
    for name,scale in [('rest',1.),('contact',.68),('wall',1.6)]:
        x=mesh.initial*scale;rest=mesh.l0*.98;direction=rng.normal(size=x.shape);direction/=np.linalg.norm(direction)
        e,g,*_=mesh.energy_force(x,rest);h=1e-6
        numerical=(mesh.energy_force(x+h*direction,rest)[0]-mesh.energy_force(x-h*direction,rest)[0])/(2*h)
        error=abs(numerical-np.sum(g*direction))/max(1,abs(numerical))
        rows.append(dict(name='energy-gradient-'+name,error=float(error),passed=bool(error<2e-7)))
        state=mesh.pack();state[:mesh.position_size]=x.ravel();state[mesh.position_size:mesh.position_size+mesh.internal_size]=rest
        state[-7:-4]=[3.,4.,20.];rhs=mesh.rhs(0,state)
        energy=lambda s:mesh.energy_force(*mesh.unpack(s)[:2])[0]+.5*p.lorenz_capacity*np.sum(mesh.unpack(s)[2]**2)
        directional=(energy(state+h*rhs)-energy(state-h*rhs))/(2*h)
        expected=rhs[-4]-rhs[-3]-rhs[-2]
        error=abs(directional-expected)/max(1,abs(expected))
        rows.append(dict(name='reciprocal-work-'+name,error=float(error),passed=bool(error<2e-6)))
    for shape in ['circle','ellipse','lobed']:
        params=replace(p,object_shape=shape);points=np.array([[1.2,.4],[-.8,.3],[.2,-.7]])
        gap,normal=object_gap(points,params);v=rng.normal(size=points.shape)
        finite=(object_gap(points+1e-6*v,params)[0]-object_gap(points-1e-6*v,params)[0])/2e-6
        error=float(max(abs(finite-np.sum(normal*v,axis=1))))
        rows.append(dict(name='obstacle-gradient-'+shape,error=error,passed=bool(error<1e-7)))
    # Independent analytic solutions for a single yielded bond and fading imprint.
    H=p.hardening;tm=p.memory_time;tp=p.yield_time;limit=p.yield_strain
    def rate(t,r,length):
        z=length-(1+H)*r[0]+H
        return [z/tm+soft_threshold(z,limit)/tp]
    length=1.5;lam=(1+H)*(1/tm+1/tp);eq=(length+H-limit*tm/(tm+tp))/(1+H)
    t=np.linspace(0,.1,201);sol=solve_ivp(lambda t,r:rate(t,r,length),(0,.1),[1.],t_eval=t,rtol=1e-12,atol=1e-14)
    exact=eq+(1-eq)*np.exp(-lam*t);error=float(max(abs(sol.y[0]-exact)))
    rows.append(dict(name='analytic-overstress-yield',error=error,passed=bool(error<1e-10)))
    t=np.linspace(0,20,201);sol=solve_ivp(lambda t,r:[-H*(r[0]-1)/tm],(0,20),[1.1],t_eval=t,rtol=1e-12,atol=1e-14)
    exact=1+.1*np.exp(-H*t/tm);error=float(max(abs(sol.y[0]-exact)))
    rows.append(dict(name='analytic-free-bond-imprint-decay',error=error,passed=bool(error<1e-10),decay_time=tm/H))
    off=replace(p,alpha=0);m,t,y,failure=integrate(off,end=2,dt=.001)
    def lorenz(t,l):
        x,z_y,z=l;return [10*(z_y-x),x*(28-z)-z_y,x*z_y-(8/3)*z]
    ref=solve_ivp(lorenz,(0,2),[1,1,1],method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
    error=float(np.max(abs(y[:,-7:-4]-ref.sol(t).T)))
    rows.append(dict(name='uncoupled-Lorenz-independent-solver',error=error,passed=bool(error<2e-6 and failure is None)))
    # Rotation must rotate the nodal motion without changing the reservoir response.
    phi=.4;R=np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    other=Mesh(replace(p,object_angle=p.object_angle+phi));state=mesh.pack();state[:mesh.position_size]=(mesh.initial*.75).ravel()
    rotated=state.copy();rotated[:mesh.position_size]=(mesh.unpack(state)[0]@R.T).ravel()
    a=mesh.rhs(0,state);b=other.rhs(0,rotated)
    error=max(float(np.max(abs(a[:mesh.position_size].reshape(-1,2)@R.T-b[:mesh.position_size].reshape(-1,2)))),float(np.max(abs(a[mesh.position_size:]-b[mesh.position_size:]))))
    rows.append(dict(name='rotational-equivariance',error=error,passed=bool(error<1e-8)))
    result=dict(checks=rows,all_passed=all(r['passed'] for r in rows))
    (DATA/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if not result['all_passed']:raise SystemExit(1)
if __name__=='__main__':main()
