"""Two-dimensional incompressible Navier-Stokes with evolving SIMS guards.

All quantities are nondimensional. Pressure is enforced by the Fourier Leray
projection. Viscosity enters through divergence of the full symmetric stress.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.fft import fft2,ifft2
from scipy.special import expit

def cases(p):
    return [dict(id=i,mode=m,weight=w,k=k,pulse=q)
            for i,(m,w,k,q) in enumerate((m,w,k,q) for m in p['modes']
                for w,k in p['weights'].items() for q in p['pulses'])]

class Model:
    def __init__(self,n,p,cs,swapped=False):
        self.n=n;self.p=p;self.cs=cs;self.swapped=swapped
        modes=np.rint(np.fft.fftfreq(n)*n).astype(int)
        self.kx=modes[:,None];self.ky=modes[None,:];self.k2=self.kx**2+self.ky**2
        self.inv=np.divide(1.,self.k2,out=np.zeros_like(self.k2,dtype=float),where=self.k2!=0)
        self.keep=(3*np.abs(self.kx)<n)&(3*np.abs(self.ky)<n)
        self.x,self.y=np.meshgrid(np.arange(n)*2*np.pi/n,np.arange(n)*2*np.pi/n,indexing='ij')
        self.k=np.array([c['k'] for c in cs])[:,:,None,None]
        self.shared=np.array([c['mode']=='shared' for c in cs])
        self.frozen=np.array([c['mode']=='frozen' for c in cs])
        self.feedback=np.array([c['mode']!='no_feedback' for c in cs])[:,None,None]
        self.chi=np.array([0. if c['mode']=='no_strain' else p['chi'] for c in cs])[:,None,None]
        self.b=np.array(p['baseline'][::-1] if swapped else p['baseline'])[None,:,None,None]
        self.force=np.zeros((len(cs),2,1,1))
        for i,c in enumerate(cs):
            if c['pulse']!='none':self.force[i,0 if c['pulse']=='A' else 1]=1
        self.profile=((1+np.cos(self.x-np.pi/2))/2)**2*((1+np.cos(self.y-np.pi/3))/2)**2
    def fft(self,f):return fft2(f,axes=(-2,-1),norm='forward',workers=1)
    def real(self,h):return ifft2(h,axes=(-2,-1),norm='forward',workers=1).real
    def project(self,h):
        z=h*self.keep
        dot=self.kx*z[:,0]+self.ky*z[:,1]
        z[:,0]-=self.kx*self.inv*dot;z[:,1]-=self.ky*self.inv*dot
        return z
    def start(self):
        x,y=self.x,self.y
        u=np.sin(x)*np.cos(y)-.25*np.sin(2*x+y)-.3*np.cos(x-2*y)
        v=-np.cos(x)*np.sin(y)+.5*np.sin(2*x+y)-.15*np.cos(x-2*y)
        a=self.p['baseline'][0]+.15*np.cos(x)+.1*np.sin(y)
        b=self.p['baseline'][1]+.1*np.sin(x)-.15*np.cos(y)
        if self.swapped:a,b=b,a
        f=np.stack([u,v,a,b,np.full_like(x,.5),np.full_like(x,.5)])
        return self.project(self.fft(np.broadcast_to(f,(len(self.cs),)+f.shape).copy()))
    def physical(self,h):
        f=self.real(h);dx=self.real(1j*self.kx*h);dy=self.real(1j*self.ky*h)
        return f,dx,dy
    def viscosity(self,f):
        active=self.p['mu']+self.p['a']*(f[:,4]+f[:,5])
        return np.where(self.feedback,active,self.p['mu']+self.p['a'])
    def rhs(self,t,h):
        f,dx,dy=self.physical(h)
        sxx=2*dx[:,0];syy=2*dy[:,1];sxy=dy[:,0]+dx[:,1]
        strain=dx[:,0]**2+dy[:,1]**2+.5*sxy**2
        nu=self.viscosity(f)
        stress=self.fft(np.stack([nu*sxx,nu*syy,nu*sxy],axis=1))*self.keep
        rhs=-self.fft(f[:,0,None]*dx+f[:,1,None]*dy)*self.keep
        rhs[:,0]+=1j*self.kx*stress[:,0]+1j*self.ky*stress[:,2]
        rhs[:,1]+=1j*self.kx*stress[:,2]+1j*self.ky*stress[:,1]
        d=f[:,3]-f[:,2]
        local=self.p['lambda']*(self.b-f[:,2:4])+self.k*f[:,4:6]*np.stack([d,-d],axis=1)
        phase=(t-self.p['pulse']['onset'])/self.p['pulse']['duration']
        if 0<phase<1:local+=self.force*self.p['pulse']['amplitude']*np.sin(np.pi*phase)**2*self.profile
        target=expit(self.p['theta']+self.p['beta']*np.stack([d,-d],axis=1)+self.chi[:,None]*strain[:,None])
        target[self.shared]=target[self.shared].mean(axis=1,keepdims=True)
        reaction=(target-f[:,4:6])/self.p['tau']
        rhs[:,2:4]+=self.fft(local)*self.keep-self.p['diffusivity_x']*self.k2*h[:,2:4]
        rhs[:,4:6]+=self.fft(reaction)*self.keep-self.p['diffusivity_g']*self.k2*h[:,4:6]
        rhs[self.frozen,4:6]=0
        diss=np.mean(2*nu*strain,axis=(-2,-1))
        aux=dict(dissipation=diss,guard_min=float(f[:,4:6].min()),guard_max=float(f[:,4:6].max()),
                 speed=float((np.abs(f[:,0])+np.abs(f[:,1])).max()),
                 divergence=float(np.abs(dx[:,0]+dy[:,1]).max()),finite=bool(np.isfinite(f).all()))
        return self.project(rhs),aux
    def independent_rhs(self,t,h):
        """Separate curl-form momentum and conservative scalar transport check."""
        out=np.zeros_like(h)
        for i,c in enumerate(self.cs):
            f=self.real(h[i]);dx=self.real(1j*self.kx*h[i]);dy=self.real(1j*self.ky*h[i])
            nu=self.p['mu']+self.p['a']*(f[4]+f[5]) if c['mode']!='no_feedback' else np.full_like(f[0],self.p['mu']+self.p['a'])
            fx=1j*self.kx*self.fft(2*nu*dx[0])+1j*self.ky*self.fft(nu*(dy[0]+dx[1]))
            fy=1j*self.kx*self.fft(nu*(dy[0]+dx[1]))+1j*self.ky*self.fft(2*nu*dy[1])
            wh=1j*self.kx*h[i,1]-1j*self.ky*h[i,0]
            wd=-self.fft(f[0]*self.real(1j*self.kx*wh)+f[1]*self.real(1j*self.ky*wh))+1j*self.kx*fy-1j*self.ky*fx
            out[i,0]=1j*self.ky*self.inv*wd;out[i,1]=-1j*self.kx*self.inv*wd
            for j in range(2,6):
                kap=self.p['diffusivity_x'] if j<4 else self.p['diffusivity_g']
                out[i,j]=-1j*self.kx*self.fft(f[0]*f[j])-1j*self.ky*self.fft(f[1]*f[j])-kap*self.k2*h[i,j]
            d=f[3]-f[2];s=dx[0]**2+dy[1]**2+.5*(dy[0]+dx[1])**2
            theta=self.p['theta'];beta=self.p['beta'];chi=0 if c['mode']=='no_strain' else self.p['chi']
            qa=1/(1+np.exp(-(theta+beta*d+chi*s)));qb=1/(1+np.exp(-(theta-beta*d+chi*s)))
            if c['mode']=='shared':qa=qb=(qa+qb)/2
            a=self.p['lambda']*(self.b[0,0,0,0]-f[2])+c['k'][0]*f[4]*d
            b=self.p['lambda']*(self.b[0,1,0,0]-f[3])-c['k'][1]*f[5]*d
            phase=(t-self.p['pulse']['onset'])/self.p['pulse']['duration']
            if 0<phase<1 and c['pulse']!='none':
                kick=self.p['pulse']['amplitude']*np.sin(np.pi*phase)**2*self.profile
                if c['pulse']=='A':a+=kick
                else:b+=kick
            out[i,2]+=self.fft(a);out[i,3]+=self.fft(b)
            out[i,4]+=self.fft((qa-f[4])/self.p['tau']);out[i,5]+=self.fft((qb-f[5])/self.p['tau'])
            if c['mode']=='frozen':out[i,4:6]=0
        return out*self.keep
    def diagnostics(self,h):
        f,dx,dy=self.physical(h);nu=self.viscosity(f)
        strain=dx[:,0]**2+dy[:,1]**2+.5*(dy[:,0]+dx[:,1])**2
        return dict(energy=.5*np.mean(f[:,0]**2+f[:,1]**2,axis=(-2,-1)),
                    dissipation=np.mean(2*nu*strain,axis=(-2,-1)),
                    enstrophy=.5*np.mean((dx[:,1]-dy[:,0])**2,axis=(-2,-1)),
                    mean_velocity=f[:,:2].mean(axis=(-2,-1)),
                    mean_guard=f[:,4:6].mean(axis=(-2,-1)),
                    mean_viscosity=nu.mean(axis=(-2,-1)),
                    gap_rms=np.sqrt(np.mean((f[:,3]-f[:,2])**2,axis=(-2,-1))),
                    gain_gap_rms=np.sqrt(np.mean((self.k[:,0]*f[:,4]-self.k[:,1]*f[:,5])**2,axis=(-2,-1))))

def run(p,cs,n,dt,progress=None,swapped=False):
    import time
    m=Model(n,p,cs,swapped);h=m.start();steps=round(p['time']['end']/dt);stride=round(p['time']['output_dt']/dt)
    assert abs(steps*dt-p['time']['end'])<1e-12 and abs(stride*dt-p['time']['output_dt'])<1e-12
    initial=m.diagnostics(h);loss=np.zeros(len(cs));records=[];snapshots=[];snap_times=[]
    saved_steps={round(t/dt) for t in p['time']['field_times']}
    look={(c['mode'],c['weight'],c['pulse']):i for i,c in enumerate(cs)}
    baselines=np.array([look.get((c['mode'],c['weight'],'none'),i) for i,c in enumerate(cs)])
    bounds=dict(guard_min=1.,guard_max=0.,max_CFL=0.,max_divergence=0.,max_energy_residual=0.,max_mean_drift=0.,checked_stages=0)
    started=time.monotonic()
    def check(aux):
        bounds['checked_stages']+=1
        bounds['guard_min']=min(bounds['guard_min'],aux['guard_min']);bounds['guard_max']=max(bounds['guard_max'],aux['guard_max'])
        bounds['max_CFL']=max(bounds['max_CFL'],aux['speed']*dt*n/(2*np.pi));bounds['max_divergence']=max(bounds['max_divergence'],aux['divergence'])
        if not aux['finite'] or aux['guard_min'] < -p['hard_limits']['guard_tolerance'] or aux['guard_max']>1+p['hard_limits']['guard_tolerance']:raise ArithmeticError('Nonfinite state or guard bound failure; no clipping allowed')
        if bounds['max_CFL']>p['hard_limits']['CFL'] or bounds['max_divergence']>p['hard_limits']['divergence']:raise ArithmeticError('CFL or incompressibility limit')
    for step in range(steps+1):
        t=step*dt
        if step%stride==0:
            d=m.diagnostics(h);f=m.real(h);delta=f-f[baselines]
            d.update(t=t,response_rms=np.sqrt(np.mean(delta**2,axis=(-2,-1))),
                     flow_response_rms=np.sqrt(np.mean(np.sum(delta[:,:2]**2,axis=1),axis=(-2,-1))),
                     energy_residual=(d['energy']-initial['energy']+loss)/initial['energy'],loss=loss.copy())
            bounds['max_energy_residual']=max(bounds['max_energy_residual'],float(np.abs(d['energy_residual']).max()))
            bounds['max_mean_drift']=max(bounds['max_mean_drift'],float(np.abs(d['mean_velocity']-initial['mean_velocity']).max()))
            if bounds['max_energy_residual']>p['hard_limits']['energy_relative'] or bounds['max_mean_drift']>p['hard_limits']['mean_drift']:raise ArithmeticError('Energy or momentum budget failure')
            records.append(d)
        if step in saved_steps:snapshots.append(m.real(h));snap_times.append(t)
        if step==steps:break
        a,aa=m.rhs(t,h);check(aa)
        b,bb=m.rhs(t+dt/2,h+dt*a/2);check(bb)
        c,cc=m.rhs(t+dt/2,h+dt*b/2);check(cc)
        dd,ee=m.rhs(t+dt,h+dt*c);check(ee)
        loss+=dt/6*(aa['dissipation']+2*bb['dissipation']+2*cc['dissipation']+ee['dissipation'])
        h=m.project(h+dt*(a+2*b+2*c+dd)/6)
        if step%50==0:
            if time.monotonic()-started>p['resources']['max_batch_seconds']:raise TimeoutError('Batch time budget')
            if progress:progress(step*dt,bounds)
    # The accepted endpoint is checked as well as every derivative stage.
    _,aux=m.rhs(p['time']['end'],h);check(aux)
    return dict(records=records,bounds=bounds,elapsed_seconds=time.monotonic()-started,n=n,dt=dt),np.array(snapshots),np.array(snap_times)
