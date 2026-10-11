"""Coupled Jacobian: spatial finite differences, exact rest/port/work blocks.

Bookkeeping variables have exactly zero derivative columns. Supplying these
known blocks avoids adaptive finite-difference factors growing on zero columns.
No force, constraint, or continuous material equation is changed here.
"""
import numpy as np
def coupled_jacobian(m,rhs,t,y):
    p=m.p;x,rest,l=m.unpack(y);k=m.position_size;e=m.internal_size;N=m.size
    f=rhs(t,y);J=np.zeros((N,N))
    # Only spatial columns need finite differences. Reverse a probe if it leaves
    # the barrier domain. The integration error controls remain independent.
    for col in range(k):
        h=1e-7*max(1,abs(y[col]));probe=y.copy();probe[col]+=h
        try:value=rhs(t,probe)
        except FloatingPointError:
            h=-h;probe=y.copy();probe[col]+=h;value=rhs(t,probe)
        J[:k,col]=(value[:k]-f[:k])/h
    mass=np.repeat(m.mass,2);delta=x[m.i]-x[m.j];length=np.linalg.norm(delta,axis=1);unit=delta/length[:,None]
    s=length-rest-p.hardening*(rest-m.l0)
    C=(1/p.memory_time+(abs(s)>p.yield_strain*m.l0)/p.yield_time) if p.plastic else np.zeros(e)
    for bond,(i,j) in enumerate(zip(m.i,m.j)):
        J[2*i:2*i+2,k+bond]=m.k[bond]*unit[bond]/(p.viscosity*m.mass[i])
        J[2*j:2*j+2,k+bond]=-m.k[bond]*unit[bond]/(p.viscosity*m.mass[j])
        J[k+bond,2*i:2*i+2]=C[bond]*unit[bond];J[k+bond,2*j:2*j+2]=-C[bond]*unit[bond]
        J[k+bond,k+bond]=-(1+p.hardening)*C[bond]
    radius=np.linalg.norm(x,axis=1);u=x/radius[:,None];lam=p.alpha/10;v=f[:k].reshape(-1,2)
    J[:k,k+e]=-lam*u.ravel()/p.viscosity
    gA=(m.mass[:,None]*u).ravel();dadot=gA@J[:k]
    dadot[:k]+=(m.mass[:,None]/radius[:,None]*(v-u*np.sum(u*v,axis=1)[:,None])).ravel()
    X,Y,Z=l;ii=k+e
    J[ii,ii]=-10;J[ii,ii+1]=10
    if p.feedback:J[ii]+=lam/p.lorenz_capacity*dadot
    J[ii+1,ii]=p.rho-Z;J[ii+1,ii+1]=-1;J[ii+1,ii+2]=-X
    J[ii+2,ii]=Y;J[ii+2,ii+1]=X;J[ii+2,ii+2]=-8/3
    J[-4,ii:ii+3]=p.lorenz_capacity*np.array([(10+p.rho)*Y-20*X,(10+p.rho)*X-2*Y,-16/3*Z])
    J[-3]=2*p.viscosity*(mass*f[:k])@J[:k]
    rate=f[k:k+e];loss_derivative=m.k*(rate+s*C)
    for bond,(i,j) in enumerate(zip(m.i,m.j)):
        J[-2,2*i:2*i+2]+=loss_derivative[bond]*unit[bond]
        J[-2,2*j:2*j+2]-=loss_derivative[bond]*unit[bond]
    J[-2,k:k+e]=-loss_derivative*(1+p.hardening)
    J[-1]=-lam*X*dadot;J[-1,ii]-=lam*float(np.dot(gA,f[:k]))
    return J
