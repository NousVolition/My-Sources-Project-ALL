"""Coupled Jacobian: spatial finite differences, exact rest/port/work blocks.

Bookkeeping variables have exactly zero derivative columns. Supplying these
known blocks avoids adaptive finite-difference factors growing on zero columns.
No force, constraint, or continuous material equation is changed here.
"""
import numpy as np
from scipy.sparse import csc_matrix
from analytic_spatial import hessian
from numerical_jacobian import coupled_jacobian as numerical
def spatial_groups(m,x):
    """Conservative dependency coloring; all boundary nodes may self-contact."""
    k=m.position_size;pattern=np.zeros((k,k),dtype=bool)
    for tri in m.tri:
        ids=np.ravel(np.column_stack([2*tri,2*tri+1]));pattern[np.ix_(ids,ids)]=True
    # A point-segment barrier is identically zero outside its support. Include
    # every pair that a finite-difference probe can activate. Two nodes move
    # by at most sqrt(2)*hmax, so point-segment distance changes by < 4*hmax.
    hmax=1e-7*max(1,float(np.max(abs(x))))
    active=m.self_gap(x,False) < m.p.self_width+8*hmax
    for ids0 in zip(m.si[active],m.sa[active],m.sb[active]):
        ids=np.ravel(np.column_stack([2*np.array(ids0),2*np.array(ids0)+1]))
        pattern[np.ix_(ids,ids)]=True
    colors=np.full(k,-1,dtype=int)
    for col in np.argsort(-pattern.sum(axis=0),kind='stable'):
        conflict=np.any(pattern[pattern[:,col]],axis=0);used=set(colors[conflict]);color=0
        while color in used:color+=1
        colors[col]=color
    groups=[np.flatnonzero(colors==c) for c in range(colors.max()+1)]
    assert all(np.max(pattern[:,g].sum(axis=1))<=1 for g in groups)
    return pattern,groups
def coupled_jacobian(m,rhs,t,y):
    if not m.p.surface_contact or not m.p.continuum_shape:return numerical(m,rhs,t,y)
    p=m.p;x,rest,l=m.unpack(y);k=m.position_size;e=m.internal_size;N=m.size
    f=rhs(t,y);J=np.zeros((N,N))
    J[:k,:k]=-hessian(m,x,rest,l[0])/(p.viscosity*np.repeat(m.mass,2)[:,None])
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
    return csc_matrix(J)
