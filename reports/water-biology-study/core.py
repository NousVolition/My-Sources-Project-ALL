"""Companion microphysics pilot. All generated fields are simulations, not samples.

Unit periodic box. Q1 uses prescribed incompressible alternating sine shears,
not a solution of unforced Navier-Stokes. Stratification uses 2D Boussinesq.
"""
import numpy as np

TWOPI = 2*np.pi


def initial_scalar(n, vertical=False):
    a = np.arange(n)/n
    x, z = np.meshgrid(a, a, indexing='ij')
    return .5 + .45*np.cos(TWOPI*(z if vertical else x))


def diffusion(c, diffusivity, dt):
    k = TWOPI*np.fft.fftfreq(c.shape[0], d=1/c.shape[0])
    return np.fft.ifft2(np.fft.fft2(c)*np.exp(-diffusivity*dt*(k[:,None]**2+k[None,:]**2))).real


def shear(c, dt, direction, phase=0., speed=1.):
    """Exact Fourier translation along each line of a single sine shear."""
    n = len(c)
    k = TWOPI*np.fft.fftfreq(n, d=1/n)
    shift = speed*dt*np.sin(TWOPI*np.arange(n)/n+phase)
    exponent = k[:,None]*shift[None,:] if direction == 0 else shift[:,None]*k[None,:]
    return np.fft.ifft(np.fft.fft(c, axis=direction)*np.exp(-1j*exponent), axis=direction).real


def schedule(t, phases, flow):
    if flow == 'still':
        return 0, 0., 0.
    half = int(np.floor(2*t+1e-9))
    return half % 2, phases[half % len(phases)], 1.


def scalar_run(n=64, dt=.01, duration=8., diffusivity=.001, seed=100, flow='stirred', save_every=.1):
    steps = round(duration/dt)
    if abs(steps*dt-duration)>1e-10 or abs(.5/dt-round(.5/dt))>1e-10:
        raise ValueError('dt must divide duration and each half-period')
    phases = np.random.default_rng(seed).uniform(0,TWOPI,16)
    c = initial_scalar(n)
    var0 = c.var()
    rows, maps = [], []
    stride = max(1,round(save_every/dt))
    for it in range(steps+1):
        if it%stride == 0 or it == steps:
            rows.append([it*dt,c.mean(),c.var(),1-c.var()/var0,
                         np.mean(2*np.minimum(c,1-c)),c.min(),c.max()])
            maps.append(c.copy())
        if it == steps:
            break
        direction,phase,speed = schedule((it+.5)*dt,phases,flow)
        c = diffusion(c,diffusivity,dt/2)
        c = shear(c,dt,direction,phase,speed)
        c = diffusion(c,diffusivity,dt/2)
    return np.asarray(rows), np.asarray(maps)


def sample_initial(rng, count):
    chunks = []
    total = 0
    while total < count:
        q = rng.random((max(2*(count-total),64),2))
        keep = rng.random(len(q)) < (1+.9*np.cos(TWOPI*q[:,0]))/1.9
        chunks.append(q[keep]); total += keep.sum()
    return np.concatenate(chunks)[:count]


def particle_run(seed=100, count=8192, dt=.01, duration=8., brownian=0., swimming=0., rotational=1., flow='stirred'):
    """Overdamped one-way particles. ABP swimming is assumed, without chemotaxis.

    Brownian half steps bracket exact shear; independent random streams separate
    initial placement, translation and orientation. No boundaries to settle onto.
    """
    rng = np.random.default_rng(seed+20000)
    orient_rng = np.random.default_rng(seed+30000)
    pos = sample_initial(np.random.default_rng(seed+10000),count)
    theta = orient_rng.uniform(0,TWOPI,count)
    exposure = np.zeros(count)
    phases = np.random.default_rng(seed).uniform(0,TWOPI,16)
    rows, hist, coords = [], [], []
    stride = round(.1/dt)
    for it in range(round(duration/dt)+1):
        if it%stride == 0:
            h = np.histogram2d(*pos.T,bins=16,range=[[0,1],[0,1]])[0]
            p = h/count
            # Multinomial noise floor removed only for CV^2; negatives retained.
            cv2 = 256*np.sum(p*p)-1
            rows.append([it*dt,cv2,cv2-255/count,1-.5*np.abs(p-1/256).sum(),
                         exposure.mean(),np.quantile(exposure,.95),len(pos)])
            hist.append(h); coords.append(pos.copy())
        if it == round(duration/dt): break
        if brownian:
            pos += np.sqrt(brownian*dt)*rng.normal(size=pos.shape)
        if swimming:
            theta += np.sqrt(rotational*dt)*orient_rng.normal(size=count)
            pos += swimming*dt/2*np.column_stack((np.cos(theta),np.sin(theta)))
        direction,phase,speed=schedule((it+.5)*dt,phases,flow)
        other=1-direction
        exposure += dt*abs(speed*TWOPI*np.cos(TWOPI*pos[:,other]+phase))
        pos[:,direction] += speed*dt*np.sin(TWOPI*pos[:,other]+phase)
        if swimming:
            pos += swimming*dt/2*np.column_stack((np.cos(theta),np.sin(theta)))
            theta += np.sqrt(rotational*dt)*orient_rng.normal(size=count)
        if brownian:
            pos += np.sqrt(brownian*dt)*rng.normal(size=pos.shape)
        pos %= 1.
    return np.asarray(rows),np.asarray(hist),np.asarray(coords),exposure


def column_run(seed, count=8192, depth=1e-3, duration=3600., radius=.5e-6, excess_density=50., mu=1e-3, adsorption_rate=0.):
    """Quiescent column: exact Stokes drift, bottom deposition, optional
    first-order bulk removal to a separately counted attached compartment.
    Adsorption is a phenomenological ablation, not resolved wall collision physics.
    """
    rng=np.random.default_rng(seed)
    z0=rng.uniform(0,depth,count)
    velocity=2*excess_density*9.81*radius**2/(9*mu)
    arrival=z0/velocity if velocity>0 else np.full(count,np.inf)
    attach=rng.exponential(1/adsorption_rate,count) if adsorption_rate>0 else np.full(count,np.inf)
    deposited=(arrival<=duration)&(arrival<attach)
    attached=(attach<=duration)&(attach<=arrival)
    suspended=~(deposited|attached)
    return velocity,np.array([suspended.sum(),deposited.sum(),attached.sum()]),z0


class Boussinesq:
    """Periodic perturbations about a linear buoyancy background, 2/3 dealiasing.

    du/dt=P[-u.grad(u)+b e_z]+nu Lap u;
    db/dt=-u.grad(b)-N2*w+kappa Lap b; dc/dt=-u.grad(c)+D Lap c.
    Positive N2 = stable, negative = unstable. Total background is not periodic.
    """
    def __init__(self,n,nu=.01,kappa=.001,diffusivity=.001,N2=0.):
        self.n=n; self.nu=nu; self.kappa=kappa; self.D=diffusivity; self.N2=N2
        m=np.fft.fftfreq(n)*n
        self.kx=TWOPI*m[:,None]; self.kz=TWOPI*m[None,:]
        self.k2=self.kx**2+self.kz**2
        self.inv=np.zeros_like(self.k2); np.divide(1,self.k2,out=self.inv,where=self.k2!=0)
        self.keep=(abs(m[:,None])<n/3)&(abs(m[None,:])<n/3)

    def project(self,h):
        dot=self.kx*h[0]+self.kz*h[1]
        h=h.copy(); h[0]-=self.kx*dot*self.inv; h[1]-=self.kz*dot*self.inv
        return h*self.keep

    def start(self,seed):
        phase=np.random.default_rng(seed).uniform(0,TWOPI,2)
        x,z=np.meshgrid(np.arange(self.n)/self.n,np.arange(self.n)/self.n,indexing='ij')
        u=.15*np.sin(TWOPI*x+phase[0])*np.cos(TWOPI*z+phase[1])
        w=-.15*np.cos(TWOPI*x+phase[0])*np.sin(TWOPI*z+phase[1])
        return self.project(np.fft.fft2(np.stack((u,w,np.zeros_like(u),initial_scalar(self.n,True)))))

    def rhs(self,h,dissipative=True):
        f=np.fft.ifft2(h).real
        dx=np.fft.ifft2(1j*self.kx*h).real
        dz=np.fft.ifft2(1j*self.kz*h).real
        r=-np.fft.fft2(f[0]*dx+f[1]*dz)*self.keep
        r[1]+=h[2]; r[2]-=self.N2*h[1]
        if dissipative:
            r-=np.array([self.nu,self.nu,self.kappa,self.D])[:,None,None]*self.k2*h
        return self.project(r)

    def run(self,seed=200,dt=.01,duration=4.):
        h=self.start(seed); rows=[]; fields=[]
        decay=np.exp(-dt/2*np.array([self.nu,self.nu,self.kappa,self.D])[:,None,None]*self.k2)
        for it in range(round(duration/dt)+1):
            if it%round(.1/dt)==0:
                f=np.fft.ifft2(h).real
                kinetic=.5*np.mean(f[0]**2+f[1]**2)
                potential=np.mean(f[2]**2)/(2*self.N2) if self.N2>0 else 0.
                div=np.max(abs(np.fft.ifft2(1j*(self.kx*h[0]+self.kz*h[1])).real))
                rows.append([it*dt,f[3].mean(),f[3].var(),kinetic,potential,div,f[3].min(),f[3].max()])
                fields.append(f)
            if it==round(duration/dt): break
            # Exact diffusion avoids the thermal explicit-diffusion restriction.
            # Strang split with RK4 for transport/buoyancy, second order overall.
            h=h*decay
            a=self.rhs(h,False); b=self.rhs(h+dt*a/2,False)
            c=self.rhs(h+dt*b/2,False); d=self.rhs(h+dt*c,False)
            h=self.project(h+dt*(a+2*b+2*c+d)/6)*decay
        return np.asarray(rows),np.asarray(fields)


def hazard_model(x, bio=0., mineral=0., surface=1., volume=1., cooling=1.):
    """Established nonhomogeneous Poisson survival, illustrative rate parameters.

    x is cooling below an UNSPECIFIED reference temperature in K. No absolute
    freezing point. All anchors are assumed, NOT measured treatment properties.
    Rate unit min^-1; cooling unit K/min. Bio slope 8.7 K^-1 informed by
    Budke & Koop (2015), local class A result; all other rate parameters assumed.
    """
    if cooling<=0 or min(bio,mineral,surface,volume)<0:
        raise ValueError('nonnegative multipliers and positive cooling required')
    x=np.asarray(x)
    rate=np.zeros_like(x,dtype=float); H=np.zeros_like(x,dtype=float)
    # half background from bulk and half from container: separately ablatable.
    for scale,slope,anchor in [(.5*volume,1.5,8),(.5*surface,1.5,8),
                               (volume*mineral,2.,6),(volume*bio,8.7,3)]:
        rate += scale*np.log(2)*slope*np.exp(slope*(x-anchor))
        H += scale*np.log(2)/cooling*(np.exp(slope*(x-anchor))-np.exp(-slope*anchor))
    return rate,H,np.exp(-H)


def bootstrap_mean(values,seed=990,reps=4000):
    values=np.asarray(values)
    draws=np.random.default_rng(seed).choice(values,(reps,len(values)),replace=True).mean(axis=1)
    return [float(values.mean()),*map(float,np.quantile(draws,[.025,.975]))]
