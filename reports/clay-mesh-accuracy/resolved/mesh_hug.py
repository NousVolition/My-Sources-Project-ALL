"""Mesh-consistent contact integration and deformation-based cell protection.

This is a revised dimensionless constitutive model, not a calibration of clay.
Historical material/contact code is retained for explicit ablations.
"""
from dataclasses import dataclass
import numpy as np
from clay_hug import Mesh as SpringMesh, object_gap, object_outline
from contact_hug import barrier
from cell_hug import Mesh as CellMesh, Parameters as CellParameters, shape_values

@dataclass(frozen=True)
class Parameters(CellParameters):
    surface_contact: bool = True
    continuum_shape: bool = True
    quadrature: int = 32
    self_neighborhood: float = .4
    def __post_init__(self):
        super().__post_init__()
        if self.quadrature < 2: raise ValueError('At least two quadrature points')
        if not np.isfinite(self.self_neighborhood) or self.self_neighborhood<=0:raise ValueError('Positive reference neighborhood')

class Mesh(CellMesh):
    def __init__(self,p):
        super().__init__(p)
        X=self.initial[self.tri]
        self.invD=np.linalg.inv(np.stack([X[:,1]-X[:,0],X[:,2]-X[:,0]],axis=2))
        self.G=np.concatenate([-self.invD.sum(axis=1)[:,None,:],self.invD],axis=1)
        q,w=np.polynomial.legendre.leggauss(p.quadrature)
        self.surface_t=(q+1)/2; self.surface_w=w/2
        self.inner_i=np.arange(p.angles);self.inner_j=np.roll(self.inner_i,-1)
        self.surface_length=np.linalg.norm(self.initial[self.inner_j]-self.initial[self.inner_i],axis=1)
        # Nonlocal surface contact has a fixed material neighborhood, independent
        # of mesh edge count. Local deformation is handled by the material energy.
        n=p.angles
        arc=lambda a,b:np.minimum(abs(a%n-b%n),n-abs(a%n-b%n))*2*np.pi/n
        same=self.si//n==self.sa//n
        distance=np.minimum(arc(self.si,self.sa),arc(self.si,self.sb))*np.linalg.norm(self.initial[self.si],axis=1)
        self.nonlocal_pair=(~same)|(distance>p.self_neighborhood)
        radii=np.linalg.norm(self.initial,axis=1)
        vertex_weight=2*radii*np.sin(np.pi/n)
        segment_weight=np.linalg.norm(self.initial[self.sb]-self.initial[self.sa],axis=1)
        # Double boundary integral, with a half factor for the two directions.
        self.surface_pair_weight=.5*vertex_weight[self.si]*segment_weight/p.self_width

    def self_force(self,x):
        gap,normal,t=self.self_gap(x) # Includes the independent crossing-domain guard.
        active=self.nonlocal_pair & (gap<self.p.self_width)
        if not np.any(active):return 0.,np.zeros_like(x)
        e,d=barrier(gap[active],self.p.self_width,self.p.self_strength)
        w=self.surface_pair_weight[active];i=self.si[active];a=self.sa[active];b=self.sb[active];s=t[active]
        force=(w*d)[:,None]*normal[active]
        grad=self.scatter(i,force)-self.scatter(a,(1-s[:,None])*force)-self.scatter(b,s[:,None]*force)
        return float(np.dot(w,e)),grad

    def deformation(self,x):
        return np.einsum('tai,taj->tij',x[self.tri],self.G)

    def material(self,x,rest):
        if not self.p.continuum_shape:return super().material(x,rest)
        energy,grad,length,ratio=SpringMesh.material(self,x,rest)
        p=self.p;area,ag=self.areas(x);J=area/self.a0
        F=self.deformation(x);S=np.sum(F*F,axis=(1,2));Q=2*J/S
        e,d=barrier(J-p.area_floor,p.area_activation-p.area_floor,p.area_strength)
        energy+=np.dot(self.a0,e)
        grad+=self.scatter(self.tri.ravel(),(d[:,None,None]*ag).reshape(-1,2))
        e,d=barrier(Q-p.shape_floor,p.shape_activation-p.shape_floor,p.shape_strength)
        energy+=np.dot(self.a0,e)
        dS=2*np.einsum('tij,taj->tai',F,self.G)
        dQ=2*(ag/self.a0[:,None,None]/S[:,None,None]-J[:,None,None]*dS/S[:,None,None]**2)
        grad+=self.scatter(self.tri.ravel(),((self.a0*d)[:,None,None]*dQ).reshape(-1,2))
        return float(energy),grad,length,ratio

    def energy_force(self,x,rest):
        if not self.p.surface_contact:return super().energy_force(x,rest)
        p=self.p;energy,grad,length,ratio=self.material(x,rest)
        if np.min(rest)<=0:raise FloatingPointError('Nonpositive internal rest length')
        # The whole-edge certificate guards the domain; it is not a quadrature rule.
        if np.min(self.clearance(x,False))<=0:raise FloatingPointError('Object edge clearance')
        radius=np.linalg.norm(x,axis=1);unit=x/np.maximum(radius[:,None],1e-12)
        actuation=float(np.dot(self.mass,radius));grad+=p.preload*self.mass[:,None]*unit;energy+=p.preload*actuation
        i,j=self.inner_i,self.inner_j;t=self.surface_t
        pts=(1-t[None,:,None])*x[i,None]+t[None,:,None]*x[j,None]
        gap,normal=object_gap(pts.reshape(-1,2),p)
        weights=(self.surface_length[:,None]*self.surface_w).ravel()
        e,d=barrier(gap,p.barrier_width,p.contact);energy+=np.dot(weights,e)
        cg=((weights*d)[:,None]*normal).reshape(len(i),len(t),2)
        grad+=self.scatter(i,((1-t[None,:,None])*cg).sum(axis=1))+self.scatter(j,(t[None,:,None]*cg).sum(axis=1))
        for sign in [-1,1]:
            e,d=barrier(p.wall_half-sign*x,p.barrier_width,p.wall)
            energy+=np.sum(self.mass[:,None]*e);grad-=sign*self.mass[:,None]*d
        e,g=self.self_force(x);energy+=e;grad+=g
        return float(energy),grad,length,ratio,unit,actuation

    def diagnostics(self,y):
        d=super().diagnostics(y);x=self.unpack(y)[0]
        F=self.deformation(x);J=self.areas(x)[0]/self.a0
        d['min_deformation_quality']=float(np.min(2*J/np.sum(F*F,axis=(1,2))))
        # Independent physical arclength weighting, not equal weight per mesh edge.
        t,w=np.polynomial.legendre.leggauss(32);t=(t+1)/2;w=w/2
        i,j=self.inner_i,self.inner_j
        points=(1-t[None,:,None])*x[i,None]+t[None,:,None]*x[j,None]
        g=object_gap(points.reshape(-1,2),self.p)[0].reshape(len(i),-1)
        lengths=np.linalg.norm(x[j]-x[i],axis=1)
        d['gap_arclength_rms']=float(np.sqrt(np.sum(lengths[:,None]*w*g*g)/lengths.sum()))
        return d

def swept_shape(m,a,b):
    from cell_hug import swept_shape as old_swept
    J,q,h=old_swept(m,a,b)
    if not m.p.continuum_shape:return J,q,h
    xa=m.unpack(a)[0];xb=m.unpack(b)[0]
    a0=m.areas(xa)[0];a1=m.areas(xb)[0];am=m.areas((xa+xb)/2)[0]
    quad=2*(a1+a0-2*am);linear=a1-a0-quad
    tau=np.clip(-linear/np.where(abs(quad)>1e-24,2*quad,1),0,1)
    amin=np.minimum(np.minimum(a0,a1),a0+linear*tau+quad*tau*tau)
    s0=np.sum(m.deformation(xa)**2,axis=(1,2));s1=np.sum(m.deformation(xb)**2,axis=(1,2))
    return J,float(np.min(2*amin/m.a0/np.maximum(s0,s1))),h
