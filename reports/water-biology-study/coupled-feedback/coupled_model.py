"""Two-way passive suspension/flow model on a two-dimensional unit torus.

Velocity changes material transport; the transported volume fraction changes
viscous stress. This is an assumed dilute constitutive model, not bacteria data.
"""
import numpy as np


class CoupledFlow:
    def __init__(self, n=32, nu0=.01, diffusivity=.001, mean_fraction=.02, feedback=True, contrast=.015):
        self.n=n; self.nu0=nu0; self.D=diffusivity
        self.mean_fraction=mean_fraction; self.feedback=feedback; self.contrast=contrast
        mode=np.rint(np.fft.fftfreq(n)*n).astype(int)
        self.kx=2*np.pi*mode[:,None]; self.ky=2*np.pi*mode[None,:]
        self.k2=self.kx**2+self.ky**2
        self.inv=np.zeros_like(self.k2);np.divide(1,self.k2,out=self.inv,where=self.k2!=0)
        self.keep=(3*abs(mode[:,None])<n)&(3*abs(mode[None,:])<n)
        self.x,self.y=np.meshgrid(np.arange(n)/n,np.arange(n)/n,indexing='ij')

    def project(self,h):
        out=h.copy();dot=self.kx*out[0]+self.ky*out[1]
        out[0]-=self.kx*dot*self.inv;out[1]-=self.ky*dot*self.inv
        return out*self.keep

    def start(self, seed=100, organization='x', disturbed=False):
        if organization not in ('x','y','uniform'):
            raise ValueError('Unknown material organization')
        # Same velocity, volume-fraction histogram and total amount in x/y arms.
        phase=np.random.default_rng(seed).uniform(0,2*np.pi,2)
        u=.15*np.sin(2*np.pi*self.x+phase[0])*np.cos(2*np.pi*self.y+phase[1])
        v=-.15*np.cos(2*np.pi*self.x+phase[0])*np.sin(2*np.pi*self.y+phase[1])
        coordinate={'x':self.x,'y':self.y}.get(organization,self.x)
        phi=self.mean_fraction+ (0 if organization=='uniform' else self.contrast*np.cos(2*np.pi*(coordinate-.25)))
        phi=np.broadcast_to(phi,(self.n,self.n)).copy()
        h=self.project(np.fft.fft2(np.stack([u,v,phi])))
        if disturbed:
            # One fixed localized divergence-free impulse; no later forcing.
            # Fix the initial Fourier polynomial across all study grids.
            xx,yy=np.meshgrid(np.arange(64)/64,np.arange(64)/64,indexing='ij')
            psi=np.exp(3*(np.cos(2*np.pi*(xx-.25))-1)+1.5*(np.cos(2*np.pi*(yy-.25))-1))
            mode=np.rint(np.fft.fftfreq(self.n)*self.n).astype(int)
            ph=np.fft.fft2(psi)[np.ix_(mode%64,mode%64)]*(self.n/64)**2
            ph*=((abs(mode[:,None])<=8)&(abs(mode[None,:])<=8)&self.keep)
            kick=np.stack([1j*self.ky*ph,-1j*self.kx*ph])*.002
            h[:2]+=kick
        return self.project(h)

    def physical(self,h):return np.fft.ifft2(h).real

    def rhs(self,h):
        f=self.physical(h);dx=self.physical(1j*self.kx*h);dy=self.physical(1j*self.ky*h)
        nu=self.nu0*(1+2.5*(f[2] if self.feedback else self.mean_fraction))
        # div[nu(phi) (grad u + grad u^T)], not nu(phi) Laplacian(u).
        sxx=2*dx[0];syy=2*dy[1];sxy=dy[0]+dx[1]
        txx=np.fft.fft2(nu*sxx)*self.keep;tyy=np.fft.fft2(nu*syy)*self.keep;txy=np.fft.fft2(nu*sxy)*self.keep
        r=-np.fft.fft2(f[0]*dx+f[1]*dy)*self.keep
        r[0]+=1j*self.kx*txx+1j*self.ky*txy
        r[1]+=1j*self.kx*txy+1j*self.ky*tyy
        r[2]-=self.D*self.k2*h[2]
        return self.project(r)

    def diagnostics(self,h):
        f=self.physical(h);dx=self.physical(1j*self.kx*h);dy=self.physical(1j*self.ky*h)
        nu=self.nu0*(1+2.5*(f[2] if self.feedback else self.mean_fraction))
        strain_sq=dx[0]**2+dy[1]**2+.5*(dy[0]+dx[1])**2
        return {'energy':float(.5*np.mean(f[0]**2+f[1]**2)),
                'dissipation':float(np.mean(2*nu*strain_sq)),
                'mass':float(f[2].mean()),'fraction_min':float(f[2].min()),'fraction_max':float(f[2].max()),
                'divergence_max':float(np.max(abs(dx[0]+dy[1]))),
                'speed_max':float(np.max(np.hypot(f[0],f[1]))),
                'mean_u':float(f[0].mean()),'mean_v':float(f[1].mean())}

    def step(self,h,dt):
        a=self.rhs(h);ha=h+dt*a/2
        b=self.rhs(ha);hb=h+dt*b/2
        c=self.rhs(hb);hc=h+dt*c
        d=self.rhs(hc)
        # Same RK4 quadrature for the kinetic-energy dissipation integral.
        loss=dt/6*(self.diagnostics(h)['dissipation']+2*self.diagnostics(ha)['dissipation']+2*self.diagnostics(hb)['dissipation']+self.diagnostics(hc)['dissipation'])
        return self.project(h+dt*(a+2*b+2*c+d)/6),loss


def run_pair(n=32,dt=.004,duration=1.,seed=100,organization='x',feedback=True,contrast=.015):
    steps=round(duration/dt);stride=round(.02/dt)
    if steps<1 or stride<1 or abs(steps*dt-duration)>1e-12 or abs(stride*dt-.02)>1e-12:
        raise ValueError('dt must divide duration and the 0.02 output interval')
    m=CoupledFlow(n,feedback=feedback,contrast=contrast)
    fields=[m.start(seed,organization,disturbed) for disturbed in [False,True]]
    energies=[m.diagnostics(h)['energy'] for h in fields];loss=[0.,0.];rows=[];saved=[]
    initial_difference=m.physical(fields[1]-fields[0])[:2]
    scale=float(np.sqrt(np.mean(np.sum(initial_difference**2,axis=0))))
    for i in range(steps+1):
        if i%stride==0:
            obs=[m.diagnostics(h) for h in fields]
            delta=m.physical(fields[1]-fields[0])[:2]
            response=float(np.sqrt(np.mean(np.sum(delta**2,axis=0))))/scale
            entry={'t':i*dt,'disturbance_velocity_L2_relative_to_initial_impulse':response,'baseline':obs[0],'disturbed':obs[1]}
            for j in [0,1]:
                obs[j]['energy_budget_relative']=(obs[j]['energy']-energies[j]+loss[j])/energies[j]
                obs[j]['cfl']=obs[j]['speed_max']*dt*n
                assert obs[j]['fraction_min']>=-1e-9 and obs[j]['fraction_max']<=.03500001
                assert obs[j]['divergence_max']<1e-10 and abs(obs[j]['mass']-.02)<1e-11
                assert abs(obs[j]['energy_budget_relative'])<1e-5 and obs[j]['cfl']<.5
            rows.append(entry);saved.append(np.stack(fields))
        if i==steps:break
        for j in [0,1]:
            fields[j],increment=m.step(fields[j],dt);loss[j]+=increment
            assert np.isfinite(fields[j]).all()
    return {'n':n,'dt':dt,'duration':duration,'seed':seed,'organization':organization,'feedback':feedback,'contrast':contrast,'initial_impulse_L2':scale,'rows':rows},np.asarray(saved)
