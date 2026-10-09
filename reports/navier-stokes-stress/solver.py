"""Unforced, constant-density 3D incompressible Navier-Stokes on [0,2*pi)^3.
Strict 2/3 Fourier truncation, pressure projection, rotational nonlinearity,
classical RK4. No molecular masses, forcing, hyperviscosity or clipping.
"""
import numpy as np
from scipy import fft


class Solver:
    def __init__(self, n, nu, workers=2):
        if nu <= 0: raise ValueError('This experiment requires positive viscosity')
        self.n, self.nu, self.workers = n, nu, workers
        self.dx = 2*np.pi/n
        m = fft.fftfreq(n)*n; p = fft.rfftfreq(n)*n
        self.k = (m[:,None,None], m[None,:,None], p[None,None,:])
        self.k2 = sum(k*k for k in self.k)
        self.inv = np.divide(1., self.k2, out=np.zeros_like(self.k2), where=self.k2>0)
        self.keep = np.maximum.reduce(np.broadcast_arrays(*[abs(k) for k in self.k])) < n/3
        self.high = np.maximum.reduce(np.broadcast_arrays(*[abs(k) for k in self.k])) >= .8*n/3
        self.weights = np.full((1,1,n//2+1), 2.); self.weights[...,0] = self.weights[...,-1] = 1.
        self.shell = np.rint(np.sqrt(self.k2)).astype(int)

    def hat(self, u): return fft.rfftn(u, axes=(-3,-2,-1), workers=self.workers)
    def real(self, h): return fft.irfftn(h, s=(self.n,)*3, axes=(-3,-2,-1), workers=self.workers)
    def inner(self, a, b): return float(np.sum((a.conj()*b).real*self.weights)/self.n**6)
    def project(self, h):
        dot = sum(k*c for k,c in zip(self.k,h))*self.inv
        return np.stack([(c-k*dot)*self.keep for k,c in zip(self.k,h)])
    def curl(self, h):
        x,y,z = self.k
        return 1j*np.stack([y*h[2]-z*h[1], z*h[0]-x*h[2], x*h[1]-y*h[0]])

    def initial(self, family, epsilon=0., seed=1):
        a=np.arange(self.n)*self.dx; x,y,z=np.meshgrid(a,a,a,indexing='ij',sparse=True)
        if family=='kida':
            u=np.stack(np.broadcast_arrays(np.sin(x)*(np.cos(3*y)*np.cos(z)-np.cos(y)*np.cos(3*z)),
                np.sin(y)*(np.cos(3*z)*np.cos(x)-np.cos(z)*np.cos(3*x)),
                np.sin(z)*(np.cos(3*x)*np.cos(y)-np.cos(x)*np.cos(3*y))))
        elif family=='taylor_green':
            u=np.stack(np.broadcast_arrays(np.sin(x)*np.cos(y)*np.cos(z),-np.cos(x)*np.sin(y)*np.cos(z),np.zeros_like(x+y+z)))
        elif family=='abc':
            u=np.stack(np.broadcast_arrays(np.sin(z)+np.cos(y),np.sin(x)+np.cos(z),np.sin(y)+np.cos(x)))
        elif family=='taylor_green_2d':
            u=np.stack(np.broadcast_arrays(np.sin(x)*np.cos(y)+0*z,-np.cos(x)*np.sin(y)+0*z,np.zeros_like(x+y+z)))
        else: raise ValueError(family)
        h=self.project(self.hat(u)); h[:,0,0,0]=0
        if epsilon:
            rng=np.random.default_rng(seed); raw=np.zeros_like(u)
            for mode in [(1,1,0),(1,0,2),(0,2,1),(2,-1,1),(1,2,-2),(2,2,1)]:
                phase=mode[0]*x+mode[1]*y+mode[2]*z+rng.uniform(0,2*np.pi)
                raw+=rng.normal(size=3)[:,None,None,None]*np.cos(phase)
            d=self.project(self.hat(raw));d[:,0,0,0]=0
            d-=h*self.inner(h,d)/self.inner(h,h)
            d*=np.sqrt(self.inner(h,h)/self.inner(d,d))
            h=np.sqrt(1-epsilon**2)*h+epsilon*d
        return h

    def rhs(self,h):
        u=self.real(h); w=self.real(self.curl(h))
        cross=np.stack([u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]])
        r=self.project(self.hat(cross))-self.nu*self.k2*h
        diss=self.nu*self.inner(h,self.k2*h)
        return r,diss,float(np.sum(abs(u),axis=0).max())

    def step(self,h,dt):
        k1,d1,u1=self.rhs(h); k2,d2,u2=self.rhs(h+dt*k1/2)
        k3,d3,u3=self.rhs(h+dt*k2/2);k4,d4,u4=self.rhs(h+dt*k3)
        out=self.project(h+dt*(k1+2*k2+2*k3+k4)/6)
        return out,dt*(d1+2*d2+2*d3+d4)/6,max(u1,u2,u3,u4)*dt/self.dx

    def observe(self,h,t,loss,e0):
        wh=self.curl(h);w=self.real(wh);u=self.real(h)
        mag=np.sqrt(np.sum(w*w,axis=0));E=.5*self.inner(h,h);Z=.5*self.inner(wh,wh)
        grad=[[self.real(1j*self.k[j]*h[i]) for j in range(3)] for i in range(3)]
        production=sum(float(np.mean(w[i]*grad[i][j]*w[j])) for i in range(3) for j in range(3))
        dz=self.nu*self.inner(wh,self.k2*wh)
        spectral=.5*np.sum(abs(h)**2,axis=0)*self.weights/self.n**6
        spectrum=np.bincount(self.shell.ravel(),weights=spectral.ravel())
        r=dict(t=float(t),energy=E,enstrophy=Z,Wmax=float(mag.max()),W_rms=float(np.sqrt(2*Z)),
               dissipation=2*self.nu*Z,stretch_production=production,enstrophy_dissipation=dz,
               enstrophy_net=production-dz,energy_residual=(E-e0+loss)/e0,
               tail_Z=self.inner(wh*self.high,wh*self.high)/max(2*Z,1e-300),
               divergence_rms=np.sqrt(self.inner(1j*sum(k*c for k,c in zip(self.k,h)),1j*sum(k*c for k,c in zip(self.k,h)))),
               mean_speed=float(np.linalg.norm(h[:,0,0,0].real/self.n**3)),
               integrated_dissipation=float(loss),peak_location=(np.array(np.unravel_index(mag.argmax(),mag.shape))*self.dx).tolist(),
               half_peak_volume_fraction=float(np.mean(mag>=.5*mag.max())),
               max_velocity=float(np.sqrt(np.sum(u*u,axis=0)).max()))
        return r,spectrum,mag[:,:,self.n//2]


def restrict(h,n):
    """Common Fourier coefficients with FFT amplitude normalization; no interpolation."""
    old=h.shape[1];m=np.rint(fft.fftfreq(n)*n).astype(int)
    return h[:,m%old][:,:,m%old][:,:,:,:n//2+1]*(n/old)**3


def difference(coarse, fine, coarse_solver):
    n=coarse_solver.n; nf=fine.shape[1]; common=restrict(fine,n)
    coarse_energy=coarse_solver.inner(coarse,coarse)
    weights=np.full((1,1,nf//2+1),2.);weights[...,0]=weights[...,-1]=1.
    fine_energy=float(np.sum(abs(fine)**2*weights)/nf**6)
    # Coarse has no retained Nyquist modes. Cross term is exact on common modes.
    delta2=max(0.,coarse_energy+fine_energy-2*coarse_solver.inner(coarse,common))
    return float(np.sqrt(delta2/max(fine_energy,1e-300)))
