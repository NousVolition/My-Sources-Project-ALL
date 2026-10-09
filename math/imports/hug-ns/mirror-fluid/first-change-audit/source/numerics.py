"""Matched-stretch study. Same supplied start and Heun time integrator.

The rotational form P(u cross curl u) equals -P((u.grad)u) for an
incompressible, dealiased Fourier field. Mean velocity is retained.
No forcing, boundary gate, pressure threshold, or feedback is introduced.
"""
import numpy as np
from scipy import fft, ndimage

L = 6.0


class Flow:
    def __init__(self, n, nu=.001, workers=4):
        self.n, self.nu, self.workers = n, nu, workers
        self.dx = L/n
        modes = fft.fftfreq(n)*n
        positive = fft.rfftfreq(n)*n
        self.k = tuple(2*np.pi/L*a for a in
                       (modes[:, None, None], modes[None, :, None], positive[None, None, :]))
        self.k2 = sum(k*k for k in self.k)
        self.inv = np.zeros_like(self.k2)
        np.divide(1, self.k2, out=self.inv, where=self.k2 > 0)
        # Same inclusive cutoff as the supplied program. Study grids are not multiples of 3.
        self.keep = ((abs(modes[:, None, None]) <= n//3) &
                     (abs(modes[None, :, None]) <= n//3) & (positive[None, None, :] <= n//3))
        self.high = ((abs(modes[:, None, None]) >= .8*(n//3)) |
                     (abs(modes[None, :, None]) >= .8*(n//3)) | (positive[None, None, :] >= .8*(n//3))) & self.keep
        self.weights = np.full((1, 1, n//2+1), 2.)
        self.weights[..., 0] = self.weights[..., -1] = 1.
        self.scale = L**3/n**6

    def hat(self, a):
        return fft.rfftn(a, axes=(-3, -2, -1), workers=self.workers)

    def real(self, a):
        return fft.irfftn(a, s=(self.n,)*3, axes=(-3, -2, -1), workers=self.workers)

    def project(self, h):
        dot = sum(k*c for k, c in zip(self.k, h))*self.inv
        return np.stack([(c-k*dot)*self.keep for k, c in zip(self.k, h)])

    def inner(self, a, b):
        return float(np.sum((a.conj()*b).real*self.weights)*self.scale)

    def curl(self, h):
        x, y, z = self.k
        return 1j*np.stack((y*h[2]-z*h[1], z*h[0]-x*h[2], x*h[1]-y*h[0]))

    def initial(self, radius=.2, spin_factor=1., strain_factor=1., angle_degrees=0.):
        axis = np.arange(self.n)*self.dx-L/2
        x, y, z = np.meshgrid(axis, axis, axis, indexing='ij', sparse=True)
        env = np.exp(-.04*(x*x+y*y+z*z))
        angle = np.deg2rad(angle_degrees)
        a = np.array([np.sin(angle), 0., np.cos(angle)])
        strain = 40*strain_factor*(-np.eye(3)+3*np.outer(a, a))
        coords = (x, y, z)
        raw = np.stack([sum(strain[i, j]*coords[j] for j in range(3))*env for i in range(3)])
        # The amplitude scaling keeps the analytic tube peak fixed when its radius changes.
        psi = np.broadcast_to(.4*spin_factor*(radius/.2)**2*np.exp(-(x*x+y*y)/radius**2), (self.n,)*3)
        ph = self.hat(psi)*self.keep
        raw[0] += self.real(1j*self.k[1]*ph)
        raw[1] -= self.real(1j*self.k[0]*ph)
        return self.project(self.hat(raw))

    def perturb(self, h, amplitude, seed):
        # A fixed set of low Fourier modes sampled identically on every grid.
        rng = np.random.default_rng(seed)
        axis = np.arange(self.n)*self.dx-L/2
        xyz = np.meshgrid(axis, axis, axis, indexing='ij', sparse=True)
        raw = np.zeros((3, self.n, self.n, self.n))
        for mode in ((1, 1, 0), (1, 0, 2), (0, 2, 1), (2, -1, 1), (1, 2, -2), (2, 2, 1)):
            phase = sum(2*np.pi/L*m*x for m, x in zip(mode, xyz))+rng.uniform(0, 2*np.pi)
            raw += rng.normal(size=3)[:, None, None, None]*np.cos(phase)
        delta = self.project(self.hat(raw))
        delta[:, 0, 0, 0] = 0
        delta *= amplitude*np.sqrt(self.inner(h, h)/self.inner(delta, delta))
        return h+delta

    def rhs(self, h):
        u = self.real(h)
        wh = self.curl(h)
        w = self.real(wh)
        force = np.stack((u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]))
        rh = self.project(self.hat(force))-self.nu*self.k2*h
        dissE = self.nu*self.inner(wh, wh)
        dE = self.inner(h, rh)
        dZ = self.inner(wh, self.curl(rh))
        dissZ = self.nu*self.inner(wh, self.k2*wh)
        rates = np.array([np.sqrt(np.sum(w*w, axis=0)).max(), dissE, dE, dZ+dissZ, dissZ])
        speed = float(np.sqrt(np.sum(u*u, axis=0)).max())
        return rh, rates, speed

    def step(self, h, dt, cached=None):
        r1, a, speed1 = self.rhs(h) if cached is None else cached
        y = h+dt*r1
        r2, b, speed2 = self.rhs(y)
        return self.project(h+.5*dt*(r1+r2)), .5*dt*(a+b), max(speed1, speed2)

    def width(self, mag, idx, direction):
        if np.linalg.norm(direction) == 0:
            return {'minimum_chord': None, 'minimum_chord_cells': None, 'censored_directions': 12}
        d = direction/np.linalg.norm(direction)
        e = np.cross(d, np.eye(3)[np.argmin(abs(d))]); e /= np.linalg.norm(e)
        f = np.cross(d, e)
        origin = np.array(idx)*self.dx
        radii = np.arange(0, L/2+self.dx/8, self.dx/4)
        level = mag[idx]/2
        widths = []
        for theta in np.linspace(0, np.pi, 12, endpoint=False):
            axis = np.cos(theta)*e+np.sin(theta)*f
            crossings = []
            for sign in (-1, 1):
                points = ((origin+sign*radii[:, None]*axis) % L/self.dx).T
                values = ndimage.map_coordinates(mag, points, order=1, mode='grid-wrap')
                hits = np.flatnonzero(values < level)
                if len(hits):
                    j = int(hits[0])
                    if j == 0: crossings.append(0.); continue
                    crossings.append(float(radii[j-1]+(level-values[j-1])*(radii[j]-radii[j-1])/(values[j]-values[j-1])))
            if len(crossings) == 2: widths.append(sum(crossings))
        return {'minimum_chord': min(widths) if widths else None,
                'maximum_chord': max(widths) if widths else None,
                'minimum_chord_cells': min(widths)/self.dx if widths else None,
                'censored_directions': 12-len(widths)}

    def observe(self, h, t, accum=None, initial=None):
        accum = np.zeros(5) if accum is None else accum
        wh = self.curl(h); w = self.real(wh)
        mag = np.sqrt(np.sum(w*w, axis=0))
        idx = np.unravel_index(np.argmax(mag), mag.shape)
        W = float(mag[idx]); avg = float(mag.mean()); rms = float(np.sqrt(np.mean(mag*mag)))
        E = .5*self.inner(h, h); Z = .5*self.inner(wh, wh)
        stretch = np.zeros_like(w)
        for i in range(3):
            for j in range(3): stretch[i] += w[j]*self.real(1j*self.k[j]*h[i])
        projected_stretch = np.sum(w*stretch, axis=0)/np.maximum(mag, 1e-300)
        visc = self.real(-self.nu*self.k2*wh)
        projected_visc = np.sum(w*visc, axis=0)/np.maximum(mag, 1e-300)
        enstrophy_production = float(np.sum(w*stretch)*self.dx**3)
        u = self.real(h)
        div = float(np.abs(self.real(1j*sum(k*c for k, c in zip(self.k, h)))).max())
        axis = np.arange(self.n)*self.dx-L/2
        roi = np.broadcast_to((axis[:, None, None]**2+axis[None, :, None]**2 <= .4**2), mag.shape)
        ci = np.unravel_index(np.argmax(np.where(roi, mag, -1)), mag.shape)
        row = {'t': float(t), 'Wmax': W, 'mean_spin': avg, 'rms_spin': rms,
               'max_to_mean_spin': W/max(avg, 1e-300), 'max_to_rms_spin': W/max(rms, 1e-300),
               'energy': E, 'enstrophy': Z, 'I': float(accum[0]),
               'peak_speed': float(np.sqrt(np.sum(u*u, axis=0)).max()), 'divergence_max': div,
               'mean_velocity': (h[:, 0, 0, 0].real/self.n**3).tolist(),
               'peak_location': [float(axis[i]) for i in idx],
               'central_roi_Wmax': float(mag[ci]), 'central_roi_peak_location': [float(axis[i]) for i in ci],
               'stretch_parallel_at_Wmax': float(projected_stretch[idx]),
               'viscosity_parallel_at_Wmax': float(projected_visc[idx]),
               'stretch_vector_at_Wmax': stretch[(slice(None),)+idx].tolist(),
               'viscosity_vector_at_Wmax': visc[(slice(None),)+idx].tolist(),
               'maximum_parallel_stretch': float(projected_stretch.max()),
               'maximum_stretch_magnitude': float(np.sqrt(np.sum(stretch*stretch, axis=0)).max()),
               'enstrophy_production': enstrophy_production,
               'viscous_energy_loss_rate': 2*self.nu*Z,
               'viscous_enstrophy_loss_rate': self.nu*self.inner(wh, self.k2*wh),
               'high_band_energy_fraction': self.inner(h*self.high, h*self.high)/max(2*E, 1e-300),
               'high_band_enstrophy_fraction': self.inner(wh*self.high, wh*self.high)/max(2*Z, 1e-300),
               'half_peak_volume': float(np.count_nonzero(mag >= W/2)*self.dx**3),
               'fixed_threshold_volumes': {str(v): float(np.count_nonzero(mag >= v)*self.dx**3) for v in (50, 100, 200, 400)},
               'width_at_global_peak': self.width(mag, idx, w[(slice(None),)+idx]),
               'width_at_central_roi_peak': self.width(mag, ci, w[(slice(None),)+ci]),
               'integrated_viscous_energy_loss': float(accum[1]),
               'integrated_enstrophy_production': float(accum[3]),
               'integrated_viscous_enstrophy_loss': float(accum[4])}
        if initial is not None:
            row['energy_budget_relative_residual'] = (E-initial['energy']+accum[1])/initial['energy']
            row['enstrophy_budget_relative_residual'] = (Z-initial['enstrophy']-accum[3]+accum[4])/max(initial['enstrophy'], Z)
        return row
