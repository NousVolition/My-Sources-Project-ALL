"""Periodic 3D incompressible NS; strict 2/3 Galerkin cutoff, double precision RK4.

Fourier projection/rotational advection adapted from the pinned stress-history
study in My-Sources-Project-ALL (f9def0b). No eddy-viscosity or particle forces.
Large-scale negative damping supplies constant power separately to each copy.
The difference-energy budget therefore explicitly includes forcing work Q.
"""
import numpy as np
from scipy import fft

L = 6.0

class Flow:
    def __init__(self, n=48, nu=.025, power=.2):
        self.n,self.nu,self.power=n,nu,power
        self.dx=L/n
        m=fft.fftfreq(n)*n; p=fft.rfftfreq(n)*n
        self.m=(m[:,None,None],m[None,:,None],p[None,None,:])
        self.k=tuple(2*np.pi/L*a for a in self.m)
        self.s2=sum(a*a for a in self.m); self.k2=sum(a*a for a in self.k)
        self.inv=np.divide(1.,self.k2,out=np.zeros_like(self.k2),where=self.k2>0)
        self.keep=(abs(self.m[0])<n/3)&(abs(self.m[1])<n/3)&(self.m[2]<n/3)
        self.low=(self.s2>0)&(self.s2<=4)
        cut=int(np.ceil(n/3)-1); self.cut=cut*2*np.pi/L
        self.high=self.keep&((abs(self.m[0])>=.8*cut)|(abs(self.m[1])>=.8*cut)|(self.m[2]>=.8*cut))
        self.weights=np.full((1,1,n//2+1),2.); self.weights[...,0]=self.weights[...,-1]=1.
        self.bins=np.floor(np.sqrt(self.s2)+.5).astype(int)
    def hat(self,u):return fft.rfftn(u,axes=(-3,-2,-1),workers=1)
    def real(self,h):return fft.irfftn(h,s=(self.n,)*3,axes=(-3,-2,-1),workers=1)
    def project(self,h):
        d=sum(self.k[j]*h[...,j,:,:,:] for j in range(3))*self.inv
        return np.stack([(h[...,j,:,:,:]-self.k[j]*d)*self.keep for j in range(3)],axis=-4)
    def curl(self,h):
        x,y,z=self.k; a,b,c=[h[...,j,:,:,:] for j in range(3)]
        return 1j*np.stack((y*c-z*b,z*a-x*c,x*b-y*a),axis=-4)
    def inner(self,a,b):return np.sum((a.conj()*b).real*self.weights,axis=(-4,-3,-2,-1))/self.n**6
    def force(self,h):
        low=h*self.low; norm=self.inner(low,low)
        if self.power==0:return np.zeros_like(h)
        if np.any(norm<1e-15):raise FloatingPointError('Empty forcing band')
        return low*(self.power/norm)[...,None,None,None,None]
    def parts(self,h):
        u=self.real(h); w=self.real(self.curl(h))
        a,b,c=[u[...,j,:,:,:] for j in range(3)]; x,y,z=[w[...,j,:,:,:] for j in range(3)]
        adv=self.project(self.hat(np.stack((b*z-c*y,c*x-a*z,a*y-b*x),axis=-4)))
        force=self.force(h); rhs=adv-self.nu*self.k2*h+force
        rates=np.stack((self.inner(h,force),self.nu*self.inner(h,self.k2*h)),axis=-1)
        if len(h)>1:
            d=h[1:]-h[0]; P=self.inner(d,adv[1:]-adv[0]); D=self.nu*self.inner(d,self.k2*d); Q=self.inner(d,force[1:]-force[0])
            difference=np.stack((P,D,Q),axis=-1)
        else:difference=np.zeros((0,3))
        speed=float(np.sqrt(np.sum(u*u,axis=-4)).max())
        return rhs,rates,difference,speed
    def step(self,h,dt):
        stages=[]; budgets=[]; differences=[]; speeds=[]
        for a in (0,.5,.5,1):
            z=h if not stages else h+a*dt*stages[-1]
            r,b,d,s=self.parts(z); stages.append(r);budgets.append(b);differences.append(d);speeds.append(s)
        weights=(1,2,2,1)
        out=self.project(h+dt/6*sum(a*b for a,b in zip(weights,stages)))
        return out,dt/6*sum(a*b for a,b in zip(weights,budgets)),dt/6*sum(a*b for a,b in zip(weights,differences)),max(speeds)*dt/self.dx

def resize(h,n):
    """Exact embedding/truncation of a Fourier polynomial, excluding Nyquist modes."""
    old=h.shape[-3]; out=np.zeros(h.shape[:-3]+(n,n,n//2+1),complex)
    c=min((old-1)//2,(n-1)//2)
    ix=np.arange(-c,c+1); iz=np.arange(c+1)
    out[...,ix[:,None,None]%n,ix[None,:,None]%n,iz[None,None,:]]=h[...,ix[:,None,None]%old,ix[None,:,None]%old,iz[None,None,:]]*(n/old)**3
    return out

def initial(f,seed,speed=1.4):
    ref=Flow(32,f.nu,0); rng=np.random.default_rng(seed)
    h=ref.project(ref.hat(rng.normal(size=(3,32,32,32))))
    h*=np.exp(-ref.s2/8)*(ref.s2<=16);h[:,0,0,0]=0
    h*=speed/np.sqrt(ref.inner(h,h))
    return f.project(resize(h,f.n))

def impulse(f,amplitude,width=1.):
    ref=Flow(64,f.nu,0); a=np.arange(64)*L/64-L/2
    x,y,z=np.meshgrid(a,a,a,indexing='ij',sparse=True); kappa=(L/(2*np.pi*width))**2
    potential=np.exp(kappa*(np.cos(2*np.pi*x/L)+np.cos(2*np.pi*y/L)+np.cos(2*np.pi*z/L)-3))
    ph=ref.hat(potential)*((abs(ref.m[0])<=5)&(abs(ref.m[1])<=5)&(ref.m[2]<=5))
    h=np.stack((1j*ref.k[1]*ph,-1j*ref.k[0]*ph,np.zeros_like(ph)))
    h=f.project(resize(h,f.n));return h*(amplitude/np.sqrt(f.inner(h,h)))

def observe(f,h):
    """Scalar budgets, spectrum, isotropic-approximation Reynolds numbers and local data."""
    density=np.sum(abs(h)**2,axis=0)*f.weights/f.n**6/2
    E=float(density.sum());eps=float(2*f.nu*np.sum(density*f.k2));u=f.real(h);w=f.real(f.curl(h))
    speed=np.sqrt(np.sum(u*u,axis=0)); adv=f.parts(h[None])[0][0]+f.nu*f.k2*h-f.force(h)
    flux=float(-f.inner(h*(f.s2<16),adv)); eta=(f.nu**3/eps)**.25
    spectrum=np.bincount(f.bins.ravel(),weights=density.ravel(),minlength=f.n)
    k=np.arange(1,len(spectrum))*2*np.pi/L
    ell=float(3*np.pi/4*np.sum(spectrum[1:]/k)/E)
    # Fourier evaluation at the fixed impulse centre, x=y=z=0 in [-L/2,L/2).
    phase=np.exp(1j*np.pi*sum(f.m)); ev=lambda a:np.sum((a*phase).real*f.weights,axis=(-3,-2,-1))/f.n**3
    localu=ev(h); grad=np.stack([ev(1j*k*h) for k in f.k],axis=1)
    s=(grad+grad.T)/2; localw=ev(f.curl(h))
    components=np.mean(u*u,axis=(1,2,3))
    global_values=[E,eps,float(np.mean(np.sum(u*w,axis=0))),float(np.mean(speed**4)),flux,
                   float(density[f.s2<=4].sum()/E),float(density[f.s2>=16].sum()/E),ell]
    local_values=[*localu,*grad.ravel(),float(localw@s@localw),float(np.linalg.norm(s))]
    stats=dict(energy=E,dissipation=eps,re_lambda=float(np.sqrt(15)*(2*E/3)/np.sqrt(f.nu*eps)),
               integral_re=float(np.sqrt(2*E)*ell/f.nu),integral_length=ell,turnover=float(ell/np.sqrt(2*E/3)),
               keta=f.cut*eta,tail_energy=float(density[f.high].sum()/E),flux_k4=flux,
               component_energies=components.tolist(),peak_speed=float(speed.max()),
               divergence=float(abs(f.real(1j*sum(k*a for k,a in zip(f.k,h)))).max()))
    return np.array(global_values+local_values),stats,spectrum
