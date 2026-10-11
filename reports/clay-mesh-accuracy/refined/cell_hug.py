"""Yielding hug with barriers against cell compression and needle-like shape.

These forces change the abstract material law; they are not position clipping,
remeshing or experimentally calibrated clay behavior. The earlier contact model
is preserved verbatim in contact_hug.py and clay_hug.py.
"""
from dataclasses import dataclass
import numpy as np
from contact_hug import Parameters as ContactParameters, Mesh as ContactMesh, barrier, object_outline, object_gap

@dataclass(frozen=True)
class Parameters(ContactParameters):
    area_floor: float = .25
    area_activation: float = .65
    area_strength: float = 100.
    shape_floor: float = .35
    shape_activation: float = .65
    shape_strength: float = 100.
    def __post_init__(self):
        super().__post_init__()
        for prefix in ['area','shape']:
            floor,active,strength=[getattr(self,prefix+'_'+s) for s in ['floor','activation','strength']]
            if not (np.isfinite([floor,active,strength]).all() and 0<floor<active<1 and strength>0):
                raise ValueError('Require 0 < floor < activation < 1 and positive strength: '+prefix)

def shape_values(m,x):
    points=x[m.tri];a,b,c=points[:,0],points[:,1],points[:,2]
    lengths2=np.stack([np.sum((b-a)**2,axis=1),np.sum((c-b)**2,axis=1),np.sum((a-c)**2,axis=1)],axis=1)
    total=lengths2.sum(axis=1);area,area_grad=m.areas(x)
    quality=4*np.sqrt(3)*area/total
    altitude=2*area/np.sqrt(lengths2.max(axis=1))
    return area,area_grad,total,quality,altitude

def geometry(m,x):
    area,_,_,q,h=shape_values(m,x);_,_,_,q0,h0=shape_values(m,m.initial)
    return dict(min_area_ratio=float(np.min(area/m.a0)),min_quality=float(np.min(q)),
        min_quality_ratio=float(np.min(q/q0)),min_altitude_ratio=float(np.min(h/h0)))

class Mesh(ContactMesh):
    def __init__(self,p):
        super().__init__(p)
        _,_,_,self.q0,self.h0=shape_values(self,self.initial)

    def material(self,x,rest):
        energy,grad,length,ratio=super().material(x,rest)
        area,ag,total,q,_=shape_values(self,x);p=self.p
        e,d=barrier(area/self.a0-p.area_floor,p.area_activation-p.area_floor,p.area_strength)
        energy+=float(np.dot(self.a0,e))
        grad+=self.scatter(self.tri.ravel(),(d[:,None,None]*ag).reshape(-1,2))
        e,d=barrier(q/self.q0-p.shape_floor,p.shape_activation-p.shape_floor,p.shape_strength)
        energy+=float(np.dot(self.a0,e))
        a,b,c=x[self.tri[:,0]],x[self.tri[:,1]],x[self.tri[:,2]]
        ds=2*np.stack([2*a-b-c,2*b-a-c,2*c-a-b],axis=1)
        dq=(4*np.sqrt(3)/self.q0)[:,None,None]*(ag/total[:,None,None]-area[:,None,None]*ds/total[:,None,None]**2)
        grad+=self.scatter(self.tri.ravel(),((self.a0*d)[:,None,None]*dq).reshape(-1,2))
        return float(energy),grad,length,ratio

    def diagnostics(self,y):
        result=super().diagnostics(y);result.update(geometry(self,self.unpack(y)[0]));return result

def swept_shape(m,a,b):
    """Bounds on area, quality and altitude along a straight accepted step.

    Area is quadratic and minimized exactly. Sum of squared edge lengths and
    squared longest-edge length are convex, hence bounded by endpoint maxima.
    Dividing the area lower bound by these maxima gives conservative bounds.
    """
    xa=m.unpack(a)[0];xb=m.unpack(b)[0]
    area0,_,s0,_,h0=shape_values(m,xa);area1,_,s1,_,h1=shape_values(m,xb)
    amid=m.areas((xa+xb)/2)[0];quad=2*(area1+area0-2*amid);linear=area1-area0-quad
    tau=np.clip(-linear/np.where(abs(quad)>1e-24,2*quad,1),0,1)
    amin=np.minimum(np.minimum(area0,area1),area0+linear*tau+quad*tau*tau)
    _,_,_,qref,href=shape_values(m,m.initial)
    qlower=4*np.sqrt(3)*amin/np.maximum(s0,s1)/qref
    longest=np.maximum(2*area0/h0,2*area1/h1)
    hlower=2*amin/longest/href
    return float(np.min(amin/m.a0)),float(np.min(qlower)),float(np.min(hlower))
