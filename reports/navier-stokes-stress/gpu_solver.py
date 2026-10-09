"""Optional CuPy backend of the same complex128 Fourier/RK4 algorithm.
The CPU solver is retained as the numerical reference and for diagnostics.
"""
import numpy as np
import cupy as cp
from solver import Solver

class GPUSolver(Solver):
    def __init__(self,n,nu,workers=2):
        self.cpu=Solver(n,nu,workers)
        self.n,self.nu,self.workers,self.dx=n,nu,workers,self.cpu.dx
        for key in ['k2','inv','keep','high','weights','shell']:setattr(self,key,cp.asarray(getattr(self.cpu,key)))
        self.k=tuple(cp.asarray(k) for k in self.cpu.k)
    def hat(self,u):return cp.fft.rfftn(cp.asarray(u),axes=(-3,-2,-1))
    def real(self,h):return cp.fft.irfftn(h,s=(self.n,)*3,axes=(-3,-2,-1))
    def inner(self,a,b):return float(cp.sum((a.conj()*b).real*self.weights)/self.n**6)
    def project(self,h):
        h=cp.asarray(h);dot=sum(k*c for k,c in zip(self.k,h))*self.inv
        return cp.stack([(c-k*dot)*self.keep for k,c in zip(self.k,h)])
    def curl(self,h):
        x,y,z=self.k
        return 1j*cp.stack([y*h[2]-z*h[1],z*h[0]-x*h[2],x*h[1]-y*h[0]])
    def initial(self,family,epsilon=0.,seed=1):
        # Exactly the same CPU-built analytic start, including random variates.
        return cp.asarray(self.cpu.initial(family,epsilon,seed))
    def rhs(self,h):
        u=self.real(h);w=self.real(self.curl(h))
        cross=cp.stack([u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]])
        r=self.project(self.hat(cross))-self.nu*self.k2*h
        return r,self.nu*self.inner(h,self.k2*h),float(cp.sum(abs(u),axis=0).max())
    def observe(self,h,t,loss,e0):return self.cpu.observe(cp.asnumpy(h),t,loss,e0)

def validate():
    checks=[]
    for n in [16,32]:
        for family in ['abc','kida','taylor_green']:
            c=Solver(n,.005);g=GPUSolver(n,.005);a=c.initial(family,.001,7);b=cp.asarray(a)
            ra,da,_=c.rhs(a);rb,db,_=g.rhs(b)
            rhs=float(np.sqrt(c.inner(ra-cp.asnumpy(rb),ra-cp.asnumpy(rb))/c.inner(ra,ra)))
            for _ in range(20):a,_,_=c.step(a,.005);b,_,_=g.step(b,.005)
            delta=a-cp.asnumpy(b);err=float(np.sqrt(c.inner(delta,delta)/c.inner(a,a)))
            checks.append(dict(n=n,family=family,rhs_relative_L2=rhs,trajectory_relative_L2=err,dissipation_difference=abs(da-db),pass_=bool(rhs<1e-11 and err<1e-11)))
    return checks

if __name__=='__main__':
    import json,time
    from pathlib import Path
    checks=validate();r=dict(cupy=cp.__version__,gpu=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),checks=checks,passed=all(c['pass_'] for c in checks))
    Path('gpu-validation.json').write_text(json.dumps(r,indent=2))
    print(json.dumps(r,indent=2));raise SystemExit(not r['passed'])
