"""A 2-D filled deformable hug with yielding bonds, fading memory and contact.

Dimensionless abstract material, not calibrated wet clay. The cohesive annular
mesh is a cross-section, not a rendered prescribed curve. A reciprocal mechanical
port cancels work exchanged between a weighted Lorenz reservoir and the material.
"""
from dataclasses import dataclass,asdict,replace
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    angles:int=24
    layers:int=3
    modulus:float=2.
    bulk:float=5.
    viscosity:float=1.5
    hardening:float=.1
    yield_strain:float=.015
    yield_time:float=.1
    memory_time:float=1.
    contact:float=100.
    wall:float=100.
    preload:float=.25
    alpha:float=.5
    feedback:bool=True
    plastic:bool=True
    rho:float=28.
    lorenz_capacity:float=.01
    object_shape:str='ellipse'
    object_angle:float=.35
    object_x:float=0.
    object_y:float=0.
    wall_half:float=2.2
    def __post_init__(self):
        if self.angles<12 or self.layers<2:raise ValueError('Mesh too small')
        for key in ('modulus','bulk','viscosity','yield_time','memory_time','contact','wall','lorenz_capacity','wall_half'):
            if not np.isfinite(getattr(self,key)) or getattr(self,key)<=0:raise ValueError(key)
        if min(self.hardening,self.yield_strain,self.alpha)<0:raise ValueError('Nonnegative material and port parameters required')
        if self.object_shape not in ('circle','ellipse','lobed'):raise ValueError('Unknown object')


def soft_threshold(x,limit):return np.sign(x)*np.maximum(abs(x)-limit,0.)


def object_gap(points,p):
    """Radial gap and its exact gradient for the specified star-shaped obstacle.

    A negative gap means overlap. Gap is exact distance for a circle; for other
    shapes it is a radial gap, not an exact closest-point distance.
    """
    points=points-np.array([p.object_x,p.object_y])
    radius=np.linalg.norm(points,axis=1);safe=np.maximum(radius,1e-12)
    angle=np.arctan2(points[:,1],points[:,0])-p.object_angle
    if p.object_shape=='circle':r=np.full(len(points),.85);dr=np.zeros(len(points))
    elif p.object_shape=='ellipse':
        a,b=1.05,.65;den=b*b*np.cos(angle)**2+a*a*np.sin(angle)**2
        r=a*b/np.sqrt(den);dr=-r*(a*a-b*b)*np.sin(angle)*np.cos(angle)/den
    else:r=.9*(1+.22*np.cos(3*angle));dr=-.594*np.sin(3*angle)
    radial=points/safe[:,None];tangent=np.column_stack([-radial[:,1],radial[:,0]])
    return radius-r,radial-(dr/safe)[:,None]*tangent


def object_outline(p,count=241):
    theta=np.linspace(-np.pi,np.pi,count);unit=np.column_stack([np.cos(theta),np.sin(theta)])
    gap,_=object_gap(unit+np.array([p.object_x,p.object_y]),p)
    return (1-gap)[:,None]*unit+np.array([p.object_x,p.object_y])


class Mesh:
    def __init__(self,p):
        self.p=p;n=p.angles
        theta=np.arange(n)*2*np.pi/n
        self.initial=np.vstack([r*np.column_stack([np.cos(theta),np.sin(theta)]) for r in np.linspace(1.25,1.65,p.layers)])
        self.n=len(self.initial)
        tris=[]
        for layer in range(p.layers-1):
            for j in range(n):
                a=layer*n+j;b=(layer+1)*n+j;c=(layer+1)*n+(j+1)%n;d=layer*n+(j+1)%n
                tris.extend([(a,b,c),(a,c,d)])
        self.tri=np.array(tris);self.a0,self.area_grad=self.areas(self.initial)
        assert np.all(self.a0>0)
        self.mass=np.bincount(self.tri.ravel(),weights=np.repeat(self.a0/3,3),minlength=self.n)
        edge_weights={}
        for tri,area in zip(self.tri,self.a0):
            for i,j in ((tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])):
                edge=tuple(sorted((i,j)));edge_weights[edge]=edge_weights.get(edge,0)+area/3
        self.edges=np.array(list(edge_weights));self.edge_weight=np.array(list(edge_weights.values()))
        self.i,self.j=self.edges.T;self.l0=np.linalg.norm(self.initial[self.i]-self.initial[self.j],axis=1)
        self.k=p.modulus*self.edge_weight/self.l0**2
        # Three-point quadrature on each inner surface segment, including curvature contact.
        self.bi=np.repeat(np.arange(n),3);self.bj=np.repeat((np.arange(n)+1)%n,3)
        self.bt=np.tile([.1127016653792583,.5,.8872983346207417],n)
        length=np.linalg.norm(self.initial[(np.arange(n)+1)%n]-self.initial[np.arange(n)],axis=1)
        self.bw=np.repeat(length,3)*np.tile([5/18,8/18,5/18],n)
        self.position_size=2*self.n;self.internal_size=len(self.edges);self.size=self.position_size+self.internal_size+7

    def areas(self,x):
        a,b,c=x[self.tri[:,0]],x[self.tri[:,1]],x[self.tri[:,2]]
        area=.5*((b[:,0]-a[:,0])*(c[:,1]-a[:,1])-(b[:,1]-a[:,1])*(c[:,0]-a[:,0]))
        grad=.5*np.stack([np.column_stack([b[:,1]-c[:,1],c[:,0]-b[:,0]]),
            np.column_stack([c[:,1]-a[:,1],a[:,0]-c[:,0]]),np.column_stack([a[:,1]-b[:,1],b[:,0]-a[:,0]])],axis=1)
        return area,grad

    def pack(self,initial=(1.,1.,1.)):
        return np.r_[self.initial.ravel(),self.l0,initial,np.zeros(4)]

    def unpack(self,y):
        k=self.position_size;e=self.internal_size
        return y[:k].reshape(self.n,2),y[k:k+e],y[k+e:k+e+3]

    def scatter(self,index,values):
        return np.column_stack([np.bincount(index,weights=values[:,j],minlength=self.n) for j in range(2)])

    def material(self,x,rest):
        p=self.p;delta=x[self.i]-x[self.j];length=np.linalg.norm(delta,axis=1)
        if np.min(length)<1e-8:raise FloatingPointError('Collapsed material bond')
        strain=length-rest;force=self.k*strain
        grad=self.scatter(self.i,force[:,None]*delta/length[:,None])-self.scatter(self.j,force[:,None]*delta/length[:,None])
        energy=.5*np.sum(self.k*(strain**2+p.hardening*(rest-self.l0)**2))
        area,ag=self.areas(x);ratio=area/self.a0
        if np.min(ratio)<=.01:raise FloatingPointError('Material cell reached inversion guard')
        energy+=p.bulk*np.sum(self.a0*(ratio-1-np.log(ratio)))
        grad+=self.scatter(self.tri.ravel(),((p.bulk*(1-1/ratio))[:,None,None]*ag).reshape(-1,2))
        return float(energy),grad,length,float(np.min(ratio))

    def energy_force(self,x,rest):
        p=self.p;energy,grad,length,ratio=self.material(x,rest)
        # Radial actuation has a potential, A=sum(node reference area * radius).
        radius=np.linalg.norm(x,axis=1);unit=x/np.maximum(radius[:,None],1e-12)
        actuation=float(np.dot(self.mass,radius));grad+=p.preload*self.mass[:,None]*unit;energy+=p.preload*actuation
        for points,weights,indices,blend in [(x,self.mass,None,None),
            ((1-self.bt[:,None])*x[self.bi]+self.bt[:,None]*x[self.bj],self.bw,(self.bi,self.bj),self.bt)]:
            gap,normal=object_gap(points,p);overlap=np.minimum(gap,0.)
            energy+=.5*p.contact*np.dot(weights,overlap**2)
            cg=p.contact*(weights*overlap)[:,None]*normal
            if indices is None:grad+=cg
            else:grad+=self.scatter(indices[0],(1-blend[:,None])*cg)+self.scatter(indices[1],blend[:,None]*cg)
        outside=np.maximum(abs(x)-p.wall_half,0.)
        energy+=.5*p.wall*np.sum(self.mass[:,None]*outside**2)
        grad+=p.wall*self.mass[:,None]*outside*np.sign(x)
        return float(energy),grad,length,ratio,unit,actuation

    def rhs(self,t,y):
        p=self.p;x,rest,lorenz=self.unpack(y);energy,gradient,length,ratio,unit,actuation=self.energy_force(x,rest)
        X,Y,Z=lorenz;lam=p.alpha/10
        velocity=(-gradient-lam*X*self.mass[:,None]*unit)/(p.viscosity*self.mass[:,None])
        adot=float(np.sum(self.mass[:,None]*unit*velocity))
        internal_force=self.k*(length-rest-p.hardening*(rest-self.l0))
        rest_rate=(internal_force/self.k/p.memory_time+soft_threshold(internal_force/self.k,p.yield_strain*self.l0)/p.yield_time) if p.plastic else np.zeros_like(rest)
        dissipation=float(p.viscosity*np.sum(self.mass[:,None]*velocity**2))
        plastic_loss=float(np.dot(internal_force,rest_rate))
        port=-lam*X*adot
        ld=np.array([10*(Y-X),X*(p.rho-Z)-Y,X*Y-(8/3)*Z])
        if p.feedback:ld[0]+=lam*adot/p.lorenz_capacity
        native=p.lorenz_capacity*((10+p.rho)*X*Y-10*X*X-Y*Y-(8/3)*Z*Z)
        return np.r_[velocity.ravel(),rest_rate,ld,native,dissipation,plastic_loss,port]

    def diagnostics(self,y):
        p=self.p;x,rest,lorenz=self.unpack(y);energy,_,length,ratio,_,_=self.energy_force(x,rest)
        fractions=np.linspace(0,1,17)
        surface=((1-fractions[None,:,None])*x[np.arange(p.angles),None]+fractions[None,:,None]*x[(np.arange(p.angles)+1)%p.angles,None]).reshape(-1,2)
        gap,_=object_gap(surface,p)
        return dict(material_energy=energy,lorenz_energy=float(.5*p.lorenz_capacity*np.dot(lorenz,lorenz)),
            imprint=float(np.sqrt(np.dot(self.edge_weight,((rest-self.l0)/self.l0)**2)/sum(self.edge_weight))),
            min_area_ratio=ratio,area_ratio=float(sum(self.areas(x)[0])/sum(self.a0)),
            penetration=float(max(0,-np.min(gap))),surface_gap_rms=float(np.sqrt(np.mean(np.maximum(gap,0)**2))),
            contact_fraction=float(np.mean(abs(gap)<.03)),wall_overlap=float(np.max(np.maximum(abs(x)-p.wall_half,0))))


def integrate(p,end=12,dt=.001,initial=(1.,1.,1.),state=None,record_step=.04,method='RK4'):
    mesh=Mesh(p);y=mesh.pack(initial) if state is None else np.array(state,copy=True);y[-4:]=0
    times=[0.];states=[y.copy()];failure=None
    if method!='RK4':
        sol=solve_ivp(mesh.rhs,(0,end),y,method=method,rtol=1e-8,atol=1e-10,max_step=dt,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        times=np.linspace(0,end,round(end/record_step)+1);states=sol.sol(times).T
        return mesh,np.array(times),states,None
    steps=int(np.ceil(end/dt));h=end/steps;stride=max(1,round(record_step/h))
    for step in range(1,steps+1):
        t=(step-1)*h
        try:
            a=mesh.rhs(t,y);b=mesh.rhs(t+h/2,y+h*a/2);c=mesh.rhs(t+h/2,y+h*b/2);d=mesh.rhs(t+h,y+h*c)
            yn=y+h*(a+2*b+2*c+d)/6
            if not np.isfinite(yn).all():raise FloatingPointError('Nonfinite state')
            mesh.energy_force(*mesh.unpack(yn)[:2]);y=yn
        except FloatingPointError as e:
            failure=dict(time=t,message=str(e),dt=h);break
        if step%stride==0 or step==steps:times.append(step*h);states.append(y.copy())
    if times[-1]<(step-1)*h:times.append((step-1)*h);states.append(y.copy())
    return mesh,np.array(times),np.array(states),failure
