"""GPU diagnostics for large refinements; identical definitions to CPU observe.
Only diagnostics move to GPU. Integration remains the previously verified code.
"""
import cupy as cp
import numpy as np
from gpu_solver import GPUSolver

class RefinementSolver(GPUSolver):
    def observe(self,h,t,loss,e0):
        wh=self.curl(h);w=self.real(wh);u=self.real(h)
        mag=cp.sqrt(cp.sum(w*w,axis=0));E=.5*self.inner(h,h);Z=.5*self.inner(wh,wh)
        # Compute one gradient at a time to limit GPU memory.
        production=0.
        for i in range(3):
            for j in range(3):production+=float(cp.mean(w[i]*self.real(1j*self.k[j]*h[i])*w[j]))
        dz=self.nu*self.inner(wh,self.k2*wh)
        spectral=.5*cp.sum(abs(h)**2,axis=0)*self.weights/self.n**6
        spectrum=cp.asnumpy(cp.bincount(self.shell.ravel(),weights=spectral.ravel()))
        div=1j*sum(k*c for k,c in zip(self.k,h))
        r=dict(t=float(t),energy=E,enstrophy=Z,Wmax=float(mag.max()),W_rms=float(np.sqrt(2*Z)),
            dissipation=2*self.nu*Z,stretch_production=production,enstrophy_dissipation=dz,
            enstrophy_net=production-dz,energy_residual=(E-e0+loss)/e0,
            tail_Z=self.inner(wh*self.high,wh*self.high)/max(2*Z,1e-300),divergence_rms=np.sqrt(self.inner(div,div)),
            mean_speed=float(cp.linalg.norm(h[:,0,0,0].real/self.n**3)),integrated_dissipation=float(loss),
            peak_location=(np.array(np.unravel_index(int(mag.argmax()),mag.shape))*self.dx).tolist(),
            half_peak_volume_fraction=float(cp.mean(mag>=.5*mag.max())),max_velocity=float(cp.sqrt(cp.sum(u*u,axis=0)).max()))
        return r,spectrum,cp.asnumpy(mag[:,:,self.n//2])

def validate_diagnostics():
    checks=[]
    for n in [32,64]:
        s=RefinementSolver(n,.005);h=s.initial('kida',.001,7)
        for _ in range(10):h,_,_=s.step(h,.01)
        e0=.5*s.inner(h,h);a,sp,sl=s.cpu.observe(cp.asnumpy(h),.1,0,e0);b,gp,gl=s.observe(h,.1,0,e0)
        err=max(abs(a[k]-b[k])/max(1.,abs(a[k])) for k in a if k!='peak_location')
        # Symmetry can tie maxima; scalar magnitude must agree, chosen index need not.
        se=float(np.max(abs(sp-gp)));le=float(np.max(abs(sl-gl)))
        checks.append(dict(n=n,max_scaled_scalar_difference=err,spectrum_max_difference=se,slice_max_difference=le,passed=err<1e-11 and se<1e-12 and le<1e-11))
    return checks
