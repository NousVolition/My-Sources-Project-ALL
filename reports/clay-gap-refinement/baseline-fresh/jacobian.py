"""Sparse assembly of the same material, rest, reciprocal-port and work blocks.

No material law or contact force changes. Dense arrays are restricted to the
three Lorenz and four work rows; local material blocks stay sparse throughout.
"""
import numpy as np
from scipy.sparse import coo_matrix,csc_matrix,hstack,vstack
from analytic_spatial import hessian
from numerical_jacobian import coupled_jacobian as numerical

def coupled_jacobian(m,rhs,t,y):
 if not m.p.surface_contact or not m.p.continuum_shape:return numerical(m,rhs,t,y)
 p=m.p;x,rest,l=m.unpack(y);k=m.position_size;e=m.internal_size;N=m.size;ii=k+e
 f=rhs(t,y);mass=np.repeat(m.mass,2)
 H=hessian(m,x,rest,l[0]);spatial=H.multiply((-1/(p.viscosity*mass))[:,None])
 delta=x[m.i]-x[m.j];length=np.linalg.norm(delta,axis=1);unit=delta/length[:,None]
 nodes=np.column_stack([m.i,m.j]);ids=(2*nodes[:,:,None]+np.arange(2)).reshape(e,4)
 signs=np.array([1.,1.,-1.,-1.]);directions=np.tile(unit,(1,2))*signs
 rr=np.broadcast_to(np.arange(e)[:,None],ids.shape)
 B=coo_matrix(((m.k[:,None]*directions/(p.viscosity*mass[ids])).ravel(),(ids.ravel(),rr.ravel())),shape=(k,e)).tocsc()
 s=length-rest-p.hardening*(rest-m.l0)
 C=(1/p.memory_time+(abs(s)>p.yield_strain*m.l0)/p.yield_time) if p.plastic else np.zeros(e)
 Rpos=coo_matrix(((C[:,None]*directions).ravel(),(rr.ravel(),ids.ravel())),shape=(e,k)).tocsc()
 Rrest=coo_matrix((-(1+p.hardening)*C,(np.arange(e),np.arange(e))),shape=(e,e)).tocsc()
 radius=np.linalg.norm(x,axis=1);u=x/radius[:,None];lam=p.alpha/10;v=f[:k].reshape(-1,2)
 port=csc_matrix(np.column_stack([-lam*u.ravel()/p.viscosity,np.zeros((k,6))]))
 V=hstack([spatial,B,port],format='csc')
 R=hstack([Rpos,Rrest,csc_matrix((e,7))],format='csc')
 gA=(m.mass[:,None]*u).ravel();dadot=np.asarray(gA@V).ravel()
 dadot[:k]+=(m.mass[:,None]/radius[:,None]*(v-u*np.sum(u*v,axis=1)[:,None])).ravel()
 X,Y,Z=l;L=np.zeros((3,N));L[:,ii:ii+3]=[[-10,10,0],[p.rho-Z,-1,-X],[Y,X,-8/3]]
 if p.feedback:L[0]+=lam/p.lorenz_capacity*dadot
 W=np.zeros((4,N));W[0,ii:ii+3]=p.lorenz_capacity*np.array([(10+p.rho)*Y-20*X,(10+p.rho)*X-2*Y,-16/3*Z])
 W[1]=2*p.viscosity*np.asarray((mass*f[:k])@V).ravel()
 rate=f[k:k+e];loss=m.k*(rate+s*C)
 np.add.at(W[2],ids.ravel(),(loss[:,None]*directions).ravel())
 W[2,k:k+e]=-loss*(1+p.hardening)
 W[3]=-lam*X*dadot;W[3,ii]-=lam*float(np.dot(gA,f[:k]))
 J=vstack([V,R,csc_matrix(L),csc_matrix(W)],format='csc');J.eliminate_zeros()
 return J
