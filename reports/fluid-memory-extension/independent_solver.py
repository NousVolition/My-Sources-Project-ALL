"""Independent CPU implementation: full complex modes, padded conservative flux,
separately written pressure projection/observations and adaptive SciPy DOP853.
No reference solver code is imported. Coefficients are normalized Fourier series.
"""
import numpy as np
from scipy import fft
from scipy.integrate import solve_ivp

class IndependentSolver:
    def __init__(self,modes,nu,injection=0):
        positive=np.asarray(modes,int)
        signed={tuple(k) for k in positive}
        signed.update(tuple(-k) for k in positive if k[2]>0)
        self.modes=np.array(sorted(signed),int);self.nu=nu;self.injection=injection
        self.lookup={tuple(k):j for j,k in enumerate(self.modes)}
        self.positive=positive;self.to_positive=np.array([self.lookup[tuple(k)] for k in positive])
        self.k2=np.sum(self.modes**2,axis=1)
        self.band=(np.max(abs(self.modes),axis=1)<=2)&(self.k2>0)
        # >3*kmax per dimension is sufficient for alias-free retained products.
        self.pad=fft.next_fast_len(3*int(abs(self.modes).max())+1)
        self.indices=tuple((self.modes[:,j]%self.pad) for j in range(3))
    def expand(self,h):
        result=np.zeros(h.shape[:-1]+(len(self.modes),),complex)
        for j,k in enumerate(self.positive):
            result[...,self.lookup[tuple(k)]]=h[...,j]
            if k[2]>0:result[...,self.lookup[tuple(-k)]]=h[...,j].conj()
        return result
    def contract(self,h):return h[...,self.to_positive]
    def project(self,h):
        dot=np.sum(h*self.modes.T[None],axis=1)
        quotient=np.divide(dot,self.k2,out=np.zeros_like(dot),where=self.k2>0)
        return h-quotient[:,None]*self.modes.T[None]
    def inner(self,a,b):return np.sum((a.conj()*b).real,axis=(-2,-1))
    def rhs(self,h):
        b=len(h);grid=np.zeros((b,3,self.pad,self.pad,self.pad),complex)
        grid[:,:,self.indices[0],self.indices[1],self.indices[2]]=h*self.pad**3
        u=fft.ifftn(grid,axes=(-3,-2,-1),workers=2).real
        nonlinear=np.zeros_like(h)
        for i in range(3):
            for j in range(i,3):
                product=fft.fftn(u[:,i]*u[:,j],axes=(-3,-2,-1),workers=2)
                coefficients=product[:,self.indices[0],self.indices[1],self.indices[2]]/self.pad**3
                nonlinear[:,i]-=1j*self.modes[:,j]*coefficients
                if i!=j:nonlinear[:,j]-=1j*self.modes[:,i]*coefficients
        force=np.zeros_like(h)
        if self.injection:
            low=h*self.band
            force=low*(self.injection/self.inner(low,low))[:,None,None]
        return self.project(nonlinear)-self.nu*self.k2*h+force, self.inner(h,force)-self.nu*self.inner(h,h*self.k2)
    def sample(self,h,x):
        phase=np.exp(1j*np.einsum('bci,mi->bcm',x,self.modes))
        velocity=np.einsum('bcm,bim->bci',phase,h).real
        gradient=np.einsum('bcm,bim,mj->bcij',phase,h,1j*self.modes).real
        return velocity,gradient
    def integrate(self,coefficients,positions,horizon,rtol=2e-10,atol=2e-12,max_step=.02):
        h=self.expand(coefficients);b=len(h);x=np.asarray(positions,float);c=x.shape[1]
        f=np.tile(np.eye(3),(b,c,1,1));hs=h.size;xs=x.size;fs=f.size
        initial=np.concatenate([h.ravel(),x.ravel(),f.ravel(),np.zeros(b)]).astype(complex)
        def unpack(y):return y[:hs].reshape(h.shape),y[hs:hs+xs].real.reshape(x.shape),y[hs+xs:hs+xs+fs].real.reshape(f.shape)
        def rhs(t,y):
            hh,xx,ff=unpack(y);dh,db=self.rhs(hh);dx,j=self.sample(hh,xx)
            return np.concatenate([dh.ravel(),dx.ravel(),(j@ff).ravel(),db])
        result=solve_ivp(rhs,(0,horizon),initial,method='DOP853',rtol=rtol,atol=atol,max_step=max_step)
        if not result.success:raise RuntimeError(result.message)
        hh,xx,ff=unpack(result.y[:,-1]);energy0=.5*self.inner(h,h);budget=result.y[-b:,-1].real
        return dict(end_fields=self.contract(hh),positions=xx,tangent=ff,nfev=result.nfev,accepted_steps=len(result.t)-1,
                    divergence=float(abs(np.sum(hh*self.modes.T[None],axis=1)).max()),
                    energy_budget=float(abs((.5*self.inner(hh,hh)-energy0-budget)/energy0).max()),
                    volume=float(abs(np.linalg.det(ff)-1).max()),times=result.t,rtol=rtol,atol=atol)
