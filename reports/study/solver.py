"""New test harness: strict dealiased Fourier and fourth-order finite differences.

FD4 uses physical-space stencils for momentum transport and viscosity. Its
periodic discrete pressure Poisson system is diagonalized with FFTs; this shared
linear-solver backend is disclosed, not described as an entirely FFT-free code.
All pressure projection/filtering is numerical; no body force or feedback is added.
"""
import numpy as np
from scipy import fft, ndimage
from initial_design import geometry, initial_parts, L


class Solver:
    def __init__(self,n,method='fourier',nu=.001,workers=4):
        self.n,self.method,self.nu,self.workers=n,method,nu,workers
        self.dx=L/n
        _,self.k,self.k2,self.keep=geometry(n)
        if method=='fourier':
            self.q=self.k;self.lap=-self.k2
        elif method=='fd4':
            self.q=tuple((8*np.sin(k*self.dx)-np.sin(2*k*self.dx))/(6*self.dx) for k in self.k)
            self.q=tuple(np.where(np.abs(q)<1e-12,0,q) for q in self.q)
            self.lap=sum((-30+32*np.cos(k*self.dx)-2*np.cos(2*k*self.dx))/(12*self.dx**2) for k in self.k)
        else:raise ValueError(method)
        q2=sum(q*q for q in self.q)
        self.inv=np.zeros_like(q2);np.divide(1,q2,out=self.inv,where=q2>0)
        self.weights=np.full((1,1,n//2+1),2.);self.weights[...,0]=1;self.weights[...,-1]=1
        self.scale=L**3/n**6

    def real(self,h):
        return fft.irfftn(h,s=(self.n,)*3,axes=(-3,-2,-1),workers=self.workers)

    def hat(self,u):
        return fft.rfftn(u,axes=(-3,-2,-1),workers=self.workers)

    def project(self,h):
        h=h*self.keep
        dot=sum(q*c for q,c in zip(self.q,h))*self.inv
        out=np.stack([c-q*dot for q,c in zip(self.q,h)])
        out[:,0,0,0]=0
        return out

    def curl(self,h,physical=False):
        q=self.k if physical else self.q
        x,y,z=q
        return 1j*np.stack((y*h[2]-z*h[1],z*h[0]-x*h[2],x*h[1]-y*h[0]))

    def inner(self,a,b):
        return float(np.sum((np.conj(a)*b).real*self.weights)*self.scale)

    def derivative(self,u,j):
        return ndimage.correlate1d(u,np.array([1,-8,0,8,-1])/(12*self.dx),axis=j,mode='wrap')

    def initial(self,case,gamma):
        core,bg,positions=initial_parts(self.n,case,gamma,self.workers)
        raw=core+bg
        projected=self.project(raw)
        change=np.sqrt(self.inner(projected-raw,projected-raw)/self.inner(raw,raw))
        return projected,positions,float(change)

    def rhs(self,h):
        h=self.project(h)
        u=self.real(h)
        wh=self.curl(h)
        w=self.real(wh)
        if self.method=='fourier':
            force=np.stack((u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]))
            rh=self.project(self.hat(force))+self.nu*self.lap*h
        else:
            force=np.zeros_like(u)
            for i in range(3):
                for j in range(3):
                    force[i]-=.5*(u[j]*self.derivative(u[i],j)+self.derivative(u[j]*u[i],j))
                    force[i]+=self.nu*ndimage.correlate1d(u[i],np.array([-1,16,-30,16,-1])/(12*self.dx**2),axis=j,mode='wrap')
            rh=self.project(self.hat(force))
        dE=self.inner(h,rh)
        dissE=-self.nu*self.inner(h,self.lap*h)
        dZ=self.inner(wh,self.curl(rh))
        dissZ=-self.nu*self.inner(wh,self.lap*wh)
        # Native-vorticity budgets plus a common physical Fourier-vorticity diagnostic.
        physical=w if self.method=='fourier' else self.real(self.curl(h,physical=True))
        W=float(np.sqrt(np.sum(physical*physical,axis=0)).max())
        nativeW=float(np.sqrt(np.sum(w*w,axis=0)).max())
        rates=np.array([W,dE,dissE,dZ,dissZ,dZ+dissZ,nativeW])
        umaxsum=float(np.sum(np.abs(u),axis=0).max())
        return rh,rates,umaxsum

    def step(self,h,dt):
        # SSP RK3; the same Butcher weights integrate diagnostic rates.
        r1,a,c1=self.rhs(h)
        y=self.project(h+dt*r1)
        r2,b,c2=self.rhs(y)
        z=self.project(.75*h+.25*(y+dt*r2))
        r3,c,c3=self.rhs(z)
        out=self.project(h/3+2*(z+dt*r3)/3)
        integral=dt*(a/6+b/6+2*c/3)
        return out,integral,max(c1,c2,c3)*dt/self.dx

    def sample(self,u,points):
        coords=((points+L/2)%L/self.dx).T
        return np.stack([ndimage.map_coordinates(v,coords,order=1,mode='grid-wrap') for v in u],axis=1)

    def move_markers(self,points,before,after,dt):
        # Passive diagnostic labels, not a force or vortex identity assumption.
        a=self.sample(before,points)
        mid=points+dt*a
        b=self.sample(after,mid)
        return (points+.5*dt*(a+b)+L/2)%L-L/2

    def width(self,mag,index,direction):
        d=direction/np.linalg.norm(direction)
        seed=np.eye(3)[np.argmin(np.abs(d))]
        e=np.cross(d,seed);e/=np.linalg.norm(e);f=np.cross(d,e)
        origin=np.asarray(index)*self.dx-L/2
        peak=float(mag[index]);level=peak/2
        radii=np.arange(0,L/2+self.dx/8,self.dx/4)
        widths=[];censored=0
        for angle in np.linspace(0,np.pi,12,endpoint=False):
            axis=np.cos(angle)*e+np.sin(angle)*f
            crossings=[]
            for sign in (-1,1):
                coords=(((origin+sign*radii[:,None]*axis)+L/2)%L/self.dx).T
                values=ndimage.map_coordinates(mag,coords,order=1,mode='grid-wrap')
                hits=np.flatnonzero(values<level)
                if not len(hits):crossings.append(None);continue
                j=int(hits[0]);assert j>0
                crossing=radii[j-1]+(level-values[j-1])*(radii[j]-radii[j-1])/(values[j]-values[j-1])
                crossings.append(float(crossing))
            if None in crossings:censored+=1
            else:widths.append(sum(crossings))
        return {'minimum_half_peak_chord':min(widths) if widths else None,
                'maximum_half_peak_chord':max(widths) if widths else None,
                'minimum_chord_cells':min(widths)/self.dx if widths else None,
                'censored_chord_directions':censored,'sample_spacing':self.dx/4}

    def observe(self,h,time,markers,accum,initial):
        u=self.real(h);wh=self.curl(h,physical=True);w=self.real(wh)
        mag=np.sqrt(np.sum(w*w,axis=0));idx=np.unravel_index(np.argmax(mag),mag.shape)
        W=float(mag[idx]);rms=float(np.sqrt(np.mean(mag*mag)))
        native=self.curl(h);En=.5*self.inner(h,h);Z=.5*self.inner(native,native)
        div=float(np.abs(self.real(1j*sum(q*c for q,c in zip(self.q,h)))).max())
        physical_div=float(np.abs(self.real(1j*sum(q*c for q,c in zip(self.k,h)))).max())
        grads=[[self.real(1j*self.k[j]*h[i]) for j in range(3)] for i in range(3)]
        x,y,z=np.meshgrid(np.arange(self.n)*self.dx-L/2,np.arange(self.n)*self.dx-L/2,np.arange(self.n)*self.dx-L/2,indexing='ij',sparse=True)
        # A fixed central-axis tube is an explicit ROI; it is not an automatic identity tracker.
        roi=np.broadcast_to((x*x+y*y)<1.0**2,mag.shape)
        cidx=np.unravel_index(np.argmax(np.where(roi,mag,-1)),mag.shape)
        cw=w[(slice(None),)+cidx]
        jac=np.array([[grads[i][j][cidx] for j in range(3)] for i in range(3)])
        S=(jac+jac.T)/2;ev,evec=np.linalg.eigh(S)
        stretching=float(cw@S@cw/(cw@cw))
        # Material points initially on the central tube complement the fixed ROI.
        # With viscosity, these points need not remain vorticity maxima.
        central_points=markers[:64]
        mw=self.sample(w,central_points)
        mj=np.stack([self.sample(np.stack(grads[i]),central_points) for i in range(3)],axis=1)
        mnorm=np.sum(mw*mw,axis=1)
        mstretch=np.einsum('ni,nij,nj->n',mw,mj,mw)/np.maximum(mnorm,1e-300)
        production=float(sum(np.mean(w[i]*grads[i][j]*w[j])*L**3 for i in range(3) for j in range(3)))
        center=np.array([g[(self.n//2,)*3] for row in grads for g in row]).reshape(3,3)
        marker_distance=[]
        core=markers[:64]
        for group in range(4):
            outer=markers[64+group*32:64+(group+1)*32]
            delta=(outer[:,None,:]-core[None,:,:]+L/2)%L-L/2
            near=np.sqrt(np.sum(delta*delta,axis=2)).min(axis=1)
            marker_distance.append({'minimum':float(near.min()),'median':float(np.median(near))})
        top=np.maximum.reduce(np.broadcast_arrays(*[np.abs(k)/(2*np.pi/L) for k in self.k]))
        high=top>=.8*(self.n/3)
        tail=self.inner(wh*high,wh*high)/max(self.inner(wh,wh),1e-300)
        row={'t':time,'Wmax':W,'W_central_roi':float(mag[cidx]),'rms_spin':rms,'spin_ratio':W/rms,
             'peak_location':[float(i*self.dx-L/2) for i in idx],
             'central_peak_location':[float(i*self.dx-L/2) for i in cidx],
             'energy':En,'native_enstrophy':Z,'spectral_enstrophy':.5*self.inner(wh,wh),
             'divergence_native':div,'divergence_physical':physical_div,
             'mean_velocity':[float(c.mean()) for c in u],
             'central_parallel_stretch':stretching,'central_strain_eigenvalues':ev.tolist(),
             'central_alignment':float(abs(cw@evec[:,-1])/np.linalg.norm(cw)),
             'central_marker_parallel_stretch_mean':float(np.mean(mstretch)),
             'central_marker_parallel_stretch_spin_weighted':float(np.sum(mstretch*mnorm)/max(np.sum(mnorm),1e-300)),
             'central_marker_parallel_stretch_min':float(np.min(mstretch)),
             'central_marker_parallel_stretch_max':float(np.max(mstretch)),
             'central_marker_Wmax':float(np.sqrt(mnorm).max()),
             'origin_axial_strain':float(center[2,2]),'direct_spectral_enstrophy_production':production,
             'outer_marker_distances':marker_distance,'high_band_enstrophy_fraction':tail,
             'I':float(accum[0]),'I_native':float(accum[6]),
             'energy_balance_residual':En-initial[0]+accum[2],
             'energy_time_integral_residual':En-initial[0]-accum[1],
             'native_enstrophy_balance_residual':Z-initial[1]-accum[5]+accum[4],
             'integrated_viscous_energy_loss':float(accum[2]),
             'integrated_native_production':float(accum[5]),'integrated_native_enstrophy_loss':float(accum[4]),
             'core_width':self.width(mag,cidx,cw),'global_width':self.width(mag,idx,w[(slice(None),)+idx])}
        for cutoff in (50,100,200,400):
            count=int(np.count_nonzero(mag>=cutoff));row[f'cells_{cutoff}']=count;row[f'volume_{cutoff}']=count*self.dx**3
        row['half_peak_cells']=int(np.count_nonzero(mag>=.5*W))
        row['half_peak_volume']=row['half_peak_cells']*self.dx**3
        return row


def seed_markers(positions):
    core=np.zeros((64,3));core[:,2]=np.linspace(-L/2,L/2,64,endpoint=False)
    groups=[core]
    for cy,cz,_ in positions:
        p=np.zeros((32,3));p[:,0]=np.linspace(-L/2,L/2,32,endpoint=False);p[:,1]=cy;p[:,2]=cz;groups.append(p)
    return np.concatenate(groups)


def validate():
    import time
    results=[]
    for method in ('fourier','fd4'):
        s=Solver(24,method);n=s.n
        rng=np.random.default_rng(42)
        h=s.project(s.hat(rng.normal(size=(3,n,n,n))))
        div=np.max(np.abs(s.real(1j*sum(q*c for q,c in zip(s.q,h)))))
        assert div<1e-11
        x=np.arange(n)*L/n
        u=np.zeros((3,n,n,n));u[0]=np.sin(2*np.pi*x[None,:,None]/L)
        h=s.hat(u);dt=.001
        new,_,_=s.step(h,dt)
        rate=(2*np.pi/L)**2 if method=='fourier' else (30-32*np.cos(2*np.pi/n)+2*np.cos(4*np.pi/n))/(12*s.dx**2)
        exact=u*np.exp(-s.nu*rate*dt)
        err=float(np.max(np.abs(s.real(new)-exact)));assert err<1e-11
        # Taylor-Green vortex: its convective acceleration is a pure gradient.
        a=2*np.pi*x/L;xx,yy,zz=np.meshgrid(a,a,a,indexing='ij')
        tg=np.stack((np.sin(xx)*np.cos(yy),-np.cos(xx)*np.sin(yy),np.zeros_like(xx)))
        ht=s.project(s.hat(tg));rhs,rates,_=s.rhs(ht)
        assert abs(rates[1]+rates[2])<1e-10
        if method=='fourier':assert np.max(np.abs(s.real(rhs)+2*s.nu*(2*np.pi/L)**2*tg))<1e-10
        results.append({'method':method,'projection_divergence':float(div),'shear_decay_error':err,'energy_semidiscrete_residual':float(rates[1]+rates[2])})
    for method in ('fourier','fd4'):
        from protocol import GAMMA
        s=Solver(80,method);h,_,change=s.initial('aligned',GAMMA)
        start=time.perf_counter();h,a,cfl=s.step(h,.00025);elapsed=time.perf_counter()-start
        results.append({'benchmark_method':method,'n':80,'step_seconds':elapsed,'cfl':cfl,'initial_projection_change':change})
    return results


if __name__=='__main__':
    import json
    print(json.dumps(validate(),indent=2))
