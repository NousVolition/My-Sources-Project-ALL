"""Version 2: conservative whole-edge barriers and guarded adaptive integration.

Uses the version-1 material and reciprocal work port, replacing soft penalties.
This is a small custom overdamped solver, not the full IPC algorithm.
"""
from dataclasses import dataclass
import numpy as np
from clay_hug import Parameters as OldParameters, Mesh as OldMesh, object_gap, object_outline


@dataclass(frozen=True)
class Parameters(OldParameters):
    barrier_width: float = .04
    contact_subdivisions: int = 8
    def __post_init__(self):
        super().__post_init__()
        if self.barrier_width <= 0 or self.contact_subdivisions < 2:
            raise ValueError('Positive barrier width and at least two subdivisions required')


def barrier(gap,width,strength):
    """Localized C2 logarithmic energy; no finite energy at nonpositive gap."""
    if np.min(gap) <= 0:
        raise FloatingPointError('Contact or wall clearance reached zero')
    d=np.minimum(gap,width);z=d-width;log=np.log(d/width)
    return -strength*z*z*log, -strength*(2*z*log+z*z/d)


class Mesh(OldMesh):
    def __init__(self,p):
        super().__init__(p)
        s=p.contact_subdivisions
        self.ci=np.repeat(self.i,s);self.cj=np.repeat(self.j,s)
        self.ct=np.tile((np.arange(s)+.5)/s,len(self.edges))
        self.cw=np.repeat(self.edge_weight/s,s)
        self.lipschitz=1.32 # sqrt(1+(.594/.702)^2) < 1.32 for lobed radial gap.

    def clearance(self,x,gradient=True):
        """Lower bounds valid on entire edges, with exact energy derivatives.

        Circle/ellipse: minimize transformed squared radius along each edge.
        Lobed: midpoint gap minus a proven Lipschitz bound on each subsegment.
        The latter is conservative, so subdivision changes the modeled offset.
        """
        p=self.p;a=x[self.i];b=x[self.j];v=b-a
        if p.object_shape!='lobed':
            c=np.cos(p.object_angle);s=np.sin(p.object_angle)
            rot=np.array([[c,-s],[s,c]])
            axes=np.array([.85,.85]) if p.object_shape=='circle' else np.array([1.05,.65])
            center=np.array([p.object_x,p.object_y]);u=(a-center)@rot/axes;w=v@rot/axes
            t=np.clip(-np.sum(u*w,axis=1)/np.maximum(np.sum(w*w,axis=1),1e-24),0,1)
            z=u+t[:,None]*w;r=np.linalg.norm(z,axis=1);scale=min(axes)
            gap=scale*(r-1)
            if not gradient:return gap
            normal=(scale*z/np.maximum(r[:,None],1e-24)/axes)@rot.T
            return gap,self.i,self.j,(1-t[:,None])*normal,t[:,None]*normal,self.edge_weight
        # Radial gradient bound applies wherever radius >= minimum object radius.
        center=np.array([p.object_x,p.object_y]);u=a-center
        t=np.clip(-np.sum(u*v,axis=1)/np.maximum(np.sum(v*v,axis=1),1e-24),0,1)
        if np.min(np.linalg.norm(u+t[:,None]*v,axis=1))<.702:
            raise FloatingPointError('Edge entered minimum-radius obstacle disk')
        points=(1-self.ct[:,None])*x[self.ci]+self.ct[:,None]*x[self.cj]
        g,n=object_gap(points,p);d=x[self.cj]-x[self.ci];length=np.linalg.norm(d,axis=1)
        k=self.lipschitz/(2*p.contact_subdivisions)
        gap=g-k*length
        if not gradient:return gap
        tangent=d/np.maximum(length[:,None],1e-24)
        return gap,self.ci,self.cj,(1-self.ct[:,None])*n+k*tangent,self.ct[:,None]*n-k*tangent,self.cw

    def energy_force(self,x,rest):
        p=self.p;energy,grad,length,ratio=self.material(x,rest)
        if np.min(rest)<=0:raise FloatingPointError('Nonpositive internal rest length')
        radius=np.linalg.norm(x,axis=1);unit=x/np.maximum(radius[:,None],1e-12)
        actuation=float(np.dot(self.mass,radius));grad+=p.preload*self.mass[:,None]*unit;energy+=p.preload*actuation
        gap,i,j,ga,gb,weight=self.clearance(x)
        e,d=barrier(gap,p.barrier_width,p.contact);energy+=float(np.dot(weight,e))
        grad+=self.scatter(i,(weight*d)[:,None]*ga)+self.scatter(j,(weight*d)[:,None]*gb)
        # Each coordinate has a separate positive/negative square-wall clearance.
        for sign in [-1,1]:
            e,d=barrier(p.wall_half-sign*x,p.barrier_width,p.wall)
            energy+=float(np.sum(self.mass[:,None]*e));grad-=sign*self.mass[:,None]*d
        return float(energy),grad,length,ratio,unit,actuation

    def swept_clearance(self,a,b):
        """Certify straight interpolation of each accepted pair, not dense RK output."""
        xa=self.unpack(a)[0];xb=self.unpack(b)[0]
        ga=self.clearance(xa,False);gb=self.clearance(xb,False)
        move=np.linalg.norm(xb-xa,axis=1)
        if self.p.object_shape=='lobed':
            bound=self.lipschitz*(1+1/self.p.contact_subdivisions)*np.maximum(move[self.ci],move[self.cj])/2
        else:bound=np.maximum(move[self.i],move[self.j])/2
        safe=float(np.min(np.minimum(ga,gb)-bound))
        # Triangle signed area is quadratic in interpolation time: check its minimum.
        a0=self.areas(xa)[0];a1=self.areas(xb)[0];amid=self.areas((xa+xb)/2)[0]
        q=2*(a1+a0-2*amid);l=a1-a0-q
        tau=np.clip(-l/np.where(abs(q)>1e-24,2*q,1),0,1)
        area=np.minimum(np.minimum(a0,a1),a0+l*tau+q*tau*tau)
        return safe,float(np.min(area/self.a0))

    def diagnostics(self,y):
        d=super().diagnostics(y);x=self.unpack(y)[0]
        d['certified_clearance']=float(np.min(self.clearance(x,False)))
        d['wall_clearance']=float(np.min(self.p.wall_half-abs(x)))
        return d


def integrate(p,end=12,dt=.004,initial=(1.,1.,1.),state=None,record_step=.04,
              rtol=2e-6,atol=2e-8,method='RK4'):
    """Step-doubled RK4, domain rejection, swept-edge and area certificates.

    Accepted path consists of two half steps. Output lies exactly on accepted
    states; integration stops at requested recording times. No position clipping.
    """
    mesh=Mesh(p);y=mesh.pack(initial) if state is None else np.array(state,copy=True);y[-4:]=0
    mesh.energy_force(*mesh.unpack(y)[:2])
    times=[0.];states=[y.copy()];t=0.;h=dt;accepted=rejected=domain_rejected=0
    minimum_clearance=float(np.min(mesh.clearance(mesh.unpack(y)[0],False)))
    minimum_swept=minimum_clearance;minimum_area=1.;min_h=dt;failure=None
    forward=returned=0.;previous_power=mesh.rhs(0,y)[-1]
    targets=np.r_[np.arange(record_step,end,record_step),end]
    def step(t,y,h):
        a=mesh.rhs(t,y);b=mesh.rhs(t+h/2,y+h*a/2);c=mesh.rhs(t+h/2,y+h*b/2);d=mesh.rhs(t+h,y+h*c)
        z=y+h*(a+2*b+2*c+d)/6
        if not np.isfinite(z).all():raise FloatingPointError('Nonfinite state')
        mesh.energy_force(*mesh.unpack(z)[:2]);return z
    for target in targets:
        while t<target-1e-12:
            h=min(h,dt,target-t)
            if h<1e-11:
                failure={'time':t,'message':'Adaptive step below minimum','step':h};break
            try:
                full=step(t,y,h);half=step(t,y,h/2);fine=step(t+h/2,half,h/2)
                scale=atol+rtol*np.maximum(abs(y),abs(fine))
                error=float(np.max(abs(fine-full)/scale)/15)
                c1,a1=mesh.swept_clearance(y,half);c2,a2=mesh.swept_clearance(half,fine)
                if min(c1,c2)<=0 or min(a1,a2)<=.01:raise FloatingPointError('Swept geometry guard')
            except FloatingPointError:
                rejected+=1;domain_rejected+=1;h*=.5;continue
            if error>1:
                rejected+=1;h*=max(.2,.9*error**(-.2));continue
            minimum_clearance=min(minimum_clearance,float(np.min(mesh.clearance(mesh.unpack(half)[0],False))),float(np.min(mesh.clearance(mesh.unpack(fine)[0],False))))
            minimum_swept=min(minimum_swept,c1,c2);minimum_area=min(minimum_area,a1,a2)
            powers=[previous_power,mesh.rhs(t+h/2,half)[-1],mesh.rhs(t+h,fine)[-1]]
            forward+=h/4*(max(powers[0],0)+2*max(powers[1],0)+max(powers[2],0))
            returned+=h/4*(max(-powers[0],0)+2*max(-powers[1],0)+max(-powers[2],0))
            previous_power=powers[2]
            min_h=min(min_h,h/2);accepted+=2;y=fine;t+=h
            h=min(dt,h*min(2.,max(.5,.9*max(error,1e-12)**(-.2))))
        times.append(t);states.append(y.copy())
        if failure:break
    audit=dict(accepted_half_steps=accepted,rejected_trials=rejected,domain_rejections=domain_rejected,
               minimum_step=min_h,minimum_accepted_clearance=minimum_clearance,minimum_swept_clearance=minimum_swept,
               minimum_swept_area_ratio=minimum_area,rtol=rtol,atol=atol,max_step=dt,
               forward_work_step_trapezoid=float(forward),return_work_step_trapezoid=float(returned))
    return mesh,np.array(times),np.array(states),failure,audit
