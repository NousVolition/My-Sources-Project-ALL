"""Batched 3D Fourier Galerkin NS and direct spectral Lagrangian observations.

Strict 2/3 truncation, positive viscosity, pressure projection, rotational
nonlinearity and synchronous RK4 for fluid, tracer positions and tangent maps.
The Fourier polynomial is evaluated directly; no grid-to-marker interpolation.
"""
import numpy as np
from scipy import fft


class Solver:
    def __init__(self, n, nu, backend='cpu'):
        if n % 2 or n < 12 or nu <= 0:
            raise ValueError('Even n >= 12 and positive viscosity required')
        self.n, self.nu, self.dx, self.backend = n, nu, 2*np.pi/n, backend
        if backend == 'gpu':
            import cupy
            self.xp = cupy
        elif backend == 'cpu':
            self.xp = np
        else:
            raise ValueError(backend)
        xp = self.xp
        m = fft.fftfreq(n)*n
        z = fft.rfftfreq(n)*n
        k = np.stack(np.broadcast_arrays(m[:, None, None], m[None, :, None], z[None, None, :]))
        self.k_cpu = np.rint(k).astype(int)
        self.keep_cpu = np.max(abs(k), axis=0) < n/3
        self.k = xp.asarray(k)
        self.k2 = xp.asarray(np.sum(k*k, axis=0))
        inv = np.divide(1., np.sum(k*k,axis=0), out=np.zeros_like(k[0]), where=np.sum(k*k,axis=0)>0)
        self.inv = xp.asarray(inv)
        self.keep = xp.asarray(self.keep_cpu)
        self.high = xp.asarray(np.max(abs(k), axis=0) >= .8*n/3)
        weights = np.full((n, n, n//2+1), 2.)
        weights[:, :, 0] = weights[:, :, -1] = 1.
        self.weights = xp.asarray(weights)
        self.modes_cpu = self.k_cpu[:, self.keep_cpu].T
        self.modes = xp.asarray(self.modes_cpu, dtype=float)
        self.eval_weights = xp.asarray(weights[self.keep_cpu]/n**3)
        self.shape = (3, n, n, n//2+1)

    def cpu(self, a):
        return self.xp.asnumpy(a) if self.backend == 'gpu' else np.asarray(a)

    def hat(self, a):
        if self.backend == 'cpu':
            return fft.rfftn(a, axes=(-3,-2,-1), workers=2)
        return self.xp.fft.rfftn(a, axes=(-3,-2,-1))

    def real(self, h):
        if self.backend == 'cpu':
            return fft.irfftn(h, s=(self.n,)*3, axes=(-3,-2,-1), workers=2)
        return self.xp.fft.irfftn(h, s=(self.n,)*3, axes=(-3,-2,-1))

    def inner(self, a, b):
        return self.xp.sum((a.conj()*b).real*self.weights, axis=(-4,-3,-2,-1))/self.n**6

    def project(self, h):
        xp = self.xp
        dot = xp.sum(h*self.k, axis=-4)*self.inv
        return (h-self.k*dot[..., None, :, :, :])*self.keep

    def curl(self, h):
        xp = self.xp
        return 1j*xp.stack([self.k[1]*h[...,2,:,:,:]-self.k[2]*h[...,1,:,:,:],
                            self.k[2]*h[...,0,:,:,:]-self.k[0]*h[...,2,:,:,:],
                            self.k[0]*h[...,1,:,:,:]-self.k[1]*h[...,0,:,:,:]], axis=-4)

    def rhs(self, h, force):
        xp = self.xp
        u = self.real(h)
        w = self.real(self.curl(h))
        cross = xp.stack([u[:,1]*w[:,2]-u[:,2]*w[:,1], u[:,2]*w[:,0]-u[:,0]*w[:,2],
                          u[:,0]*w[:,1]-u[:,1]*w[:,0]], axis=1)
        rate = self.inner(h, force)-self.nu*self.inner(h, self.k2*h)
        return self.project(self.hat(cross))-self.nu*self.k2*h+force, rate, xp.max(xp.sum(abs(u), axis=1), axis=(-3,-2,-1))

    def sample(self, h, positions, centers):
        """Exact retained Fourier polynomial and its analytic spatial derivative."""
        xp = self.xp
        coeff = xp.moveaxis(h[:, :, self.keep], 1, 2)*self.eval_weights[None,:,None]
        phase = xp.exp(1j*xp.matmul(positions, self.modes.T))
        velocity = xp.matmul(phase, coeff).real
        derivatives = (coeff[:,:,:,None]*(1j*self.modes[None,:,None,:])).reshape(len(h),len(self.modes_cpu),9)
        gradient = xp.matmul(phase[:,:centers], derivatives).real.reshape(len(h),centers,3,3)
        return velocity, gradient

    def step(self, h, x, tangent, dt, force=None):
        xp = self.xp
        if force is None:
            force = xp.zeros_like(h)
        c = tangent.shape[1]
        def stage(hh, xx, ff):
            r, budget, speed = self.rhs(hh, force)
            v, j = self.sample(hh, xx, c)
            return r, v, xp.matmul(j, ff), budget, speed
        a = stage(h,x,tangent)
        b = stage(h+dt*a[0]/2,x+dt*a[1]/2,tangent+dt*a[2]/2)
        cstage = stage(h+dt*b[0]/2,x+dt*b[1]/2,tangent+dt*b[2]/2)
        d = stage(h+dt*cstage[0],x+dt*cstage[1],tangent+dt*cstage[2])
        out = [sum(s[i]*weight for s,weight in zip([a,b,cstage,d],[1,2,2,1]))*dt/6 for i in range(4)]
        return self.project(h+out[0]), x+out[1], tangent+out[2], out[3], xp.maximum(xp.maximum(a[4],b[4]),xp.maximum(cstage[4],d[4]))*dt/self.dx

    def pack(self, h):
        return self.cpu(h[:, :, self.keep])/self.n**3

    def unpack(self, coefficients):
        h = self.xp.zeros((len(coefficients),)+self.shape, dtype=complex)
        h[:,:,self.keep] = self.xp.asarray(coefficients)*self.n**3
        return h

    def diagnostics(self, h, x, tangent, budget, energy0, cfl):
        xp = self.xp
        wh = self.curl(h)
        energy = .5*self.inner(h,h)
        div = 1j*xp.sum(self.k*h, axis=1)
        div2 = xp.sum(abs(div)**2*self.weights, axis=(-3,-2,-1))/self.n**6
        _, j = self.sample(h,x,tangent.shape[1])
        deterror = xp.max(abs(xp.linalg.det(tangent)-1), axis=1)
        return self.cpu(xp.stack([energy, (energy-energy0-budget)/xp.maximum(energy0,1e-30),
                                  xp.sqrt(div2), self.inner(wh*self.high,wh*self.high)/xp.maximum(self.inner(wh,wh),1e-30),
                                  deterror, xp.max(abs(xp.trace(j,axis1=-2,axis2=-1)),axis=1), cfl],axis=1))


def random_field(s, seed, rms=.6):
    """Identical analytic low-mode polynomial on every grid; not grid noise."""
    rng = np.random.default_rng(seed)
    axes = np.arange(s.n)*s.dx
    x,y,z = np.meshgrid(axes,axes,axes,indexing='ij',sparse=True)
    u = np.zeros((3,s.n,s.n,s.n))
    for mode in [(a,b,c) for a in range(3) for b in range(-2,3) for c in range(-2,3)
                 if (a,b,c)!=(0,0,0) and (a>0 or b>0 or (b==0 and c>0))]:
        phase = mode[0]*x+mode[1]*y+mode[2]*z+rng.uniform(0,2*np.pi)
        vector = rng.normal(size=3)/(1+sum(k*k for k in mode))
        u += vector[:,None,None,None]*np.cos(phase)
    h = s.project(s.hat(s.xp.asarray(u))[None])[0]
    h[:,0,0,0] = 0
    return h*rms/float(s.cpu(s.xp.sqrt(s.inner(h,h))))


def curl_gaussian(s, center, width, cutoff=4, high_only=False):
    """Periodized Gaussian vector potential A_z, analytically sampled in Fourier space."""
    xp = s.xp
    psi = xp.exp(-.5*width**2*s.k2)*xp.exp(-1j*sum(s.k[j]*center[j] for j in range(3)))*s.n**3
    psi *= xp.max(abs(s.k),axis=0) <= cutoff
    if high_only:
        psi *= xp.max(abs(s.k),axis=0) > 2
    h = xp.stack([1j*s.k[1]*psi,-1j*s.k[0]*psi,xp.zeros_like(psi)])
    h = s.project(h)
    norm = float(s.cpu(xp.sqrt(s.inner(h,h))))
    return h/norm


def phase_values(modes, seed):
    """Grid-independent hashed phases, odd under k -> -k to preserve reality."""
    k = np.asarray(modes, dtype=np.int64)
    nonzero = k != 0
    first = nonzero.argmax(axis=1)
    sign = np.sign(k[np.arange(len(k)),first]); sign[~nonzero.any(axis=1)] = 1
    canonical = k*sign[:,None]
    with np.errstate(over='ignore'):
        a = np.uint64(seed)+np.uint64(canonical[:,0]+128)*np.uint64(73856093)+np.uint64(canonical[:,1]+128)*np.uint64(19349663)+np.uint64(canonical[:,2]+128)*np.uint64(83492791)
        a = (a+np.uint64(0x9E3779B97F4A7C15))
        a = (a^(a>>np.uint64(30)))*np.uint64(0xBF58476D1CE4E5B9)
        a = (a^(a>>np.uint64(27)))*np.uint64(0x94D049BB133111EB)
        a = a^(a>>np.uint64(31))
    return sign*((a>>np.uint64(11)).astype(float)/2**53*2*np.pi)


def interventions(s, prepared, seed, shift, replicates=4):
    xp = s.xp
    low = xp.max(abs(s.k),axis=0) <= 2
    high = s.keep & ~low
    common_low = prepared.mean(axis=0)*low
    highfields = prepared*high
    power = xp.mean(xp.sum(abs(highfields)**2,axis=1),axis=0)
    norm = xp.sqrt(xp.sum(abs(highfields)**2,axis=1))
    # Identical power per mode and common coarse field; retain each mode's vector direction.
    normalized = highfields*xp.sqrt(power)[None,None]/xp.maximum(norm[:,None],1e-300)
    retained = common_low[None]+normalized
    variants = {'retained':retained}
    for r in range(replicates):
        factors = xp.ones_like(s.k2,dtype=complex)
        factors[s.keep] = xp.asarray(np.exp(1j*phase_values(s.modes_cpu, seed+910000+1009*r)))
        variants['scramble'+str(r)] = common_low[None]+normalized*factors[None,None]
    phase = xp.exp(1j*sum(s.k[j]*shift[j] for j in range(3)))
    variants['sham'] = common_low[None]+normalized*phase[None,None]
    # A neutral, transverse reference. Averaging opposite vectors then normalizing
    # would amplify cancellation roundoff and can destroy conjugate symmetry.
    k0,k1,k2=s.k
    reset=1j*xp.stack([3*k1-2*k2,k2-3*k0,2*k0-k1])
    alternate=1j*xp.stack([k1+k2,2*k2-k0,-k0-2*k1])
    reset=xp.where((xp.sum(abs(reset)**2,axis=0)==0)[None],alternate,reset)
    normreset = xp.sqrt(xp.sum(abs(reset)**2,axis=0))
    reset *= xp.sqrt(power)[None]/xp.maximum(normreset[None],1e-300)
    neutral_phase=xp.ones_like(s.k2,dtype=complex)
    neutral_phase[s.keep]=xp.asarray(np.exp(1j*phase_values(s.modes_cpu,seed+800000)))
    reset*=neutral_phase[None]
    variants['reset'] = xp.stack([common_low+reset,common_low+reset])
    checks = {}
    reference_power = xp.sum(abs(retained)**2,axis=1)
    for key,value in variants.items():
        checks[key] = dict(power_error=float(s.cpu(xp.max(abs(xp.sum(abs(value)**2,axis=1)-reference_power)))/s.n**6),
                           coarse_error=float(s.cpu(xp.max(abs((value-retained)*low)))/s.n**3),
                           divergence_error=float(s.cpu(xp.max(abs(xp.sum(s.k*value,axis=1))))/s.n**3))
    return variants, checks
