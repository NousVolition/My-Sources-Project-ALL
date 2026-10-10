"""Full-inertia 2D Boussinesq pilot with stress-free, fixed-temperature plates.

u_t + u.grad(u) = -grad(p) + Pr*Lap(u) + Pr*Ra*theta*e_z
theta_t + u.grad(theta) = Lap(theta) + w, div(u)=0.
Length H, time H^2/kappa, temperature imposed DeltaT. Synthetic model only.
"""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
PLAN = json.loads((ROOT/'plan.json').read_text())
LX = 2*np.sqrt(2)
KC = np.pi/np.sqrt(2)
Q = KC**2 + np.pi**2
RAC = 27*np.pi**4/4


def modal_matrix(ra, pr=1., k=KC, vertical=1):
    q = k*k+(vertical*np.pi)**2
    return np.array([[-pr*q, pr*ra*k*k/q], [1., -q]])


def leading_rate(ra, pr=1.):
    return float(np.max(np.linalg.eigvals(modal_matrix(ra, pr)).real))


class Layer:
    def __init__(self, n=24, ra=RAC, pr=1., buoyancy=True, nonlinear=True):
        assert n % 2 == 0 and n >= 12
        self.n, self.ra, self.pr = n, ra, pr
        self.buoyancy, self.nonlinear = buoyancy, nonlinear
        m = np.fft.fftfreq(n)*n
        self.kx = (2*np.pi/LX)*m[:, None]
        self.kz = np.pi*m[None, :]
        self.k2 = self.kx**2+self.kz**2
        self.inv = np.divide(1., self.k2, out=np.zeros_like(self.k2), where=self.k2>0)
        self.keep = (abs(m[:, None]) < n/3) & (abs(m[None, :]) < n/3)
        self.reflect_z = (-np.arange(n)) % n
        self.x, self.z = np.meshgrid(np.arange(n)*LX/n, np.arange(n)*2/n, indexing='ij')

    def project(self, h):
        h = h.copy()*self.keep
        # Enforce the plate boundary subspace at EVERY stage. Analytical parity
        # alone is insufficient: roundoff in forbidden kz=0 vertical motions
        # can grow much faster than physical bounded-layer modes.
        mirror = h[:, :, self.reflect_z]
        h[0] = .5*(h[0]+mirror[0])
        h[1:] = .5*(h[1:]-mirror[1:])
        dot = self.kx*h[0]+self.kz*h[1]
        h[0] -= self.kx*dot*self.inv
        h[1] -= self.kz*dot*self.inv
        h[:, 0, 0] = 0.  # zero mean velocity; odd thermal extension has zero mean
        return h

    def initial(self, seed=410, amplitude=1e-8, mixed=True):
        rng = np.random.default_rng(seed)
        f = np.zeros((3, self.n, self.n))
        modes = [(1, 1, 1.)]+([(2, 1, .1), (1, 2, .1), (2, 2, .1), (3, 2, .1)] if mixed else [])
        for mx, mz, weight in modes:
            phase = rng.uniform(0, 2*np.pi)
            theta = amplitude*weight
            w = Q*theta
            a = mx*KC*self.x+phase
            z = mz*np.pi*self.z
            f[0] += -mz*np.pi/(mx*KC)*w*np.sin(a)*np.cos(z)
            f[1] += w*np.cos(a)*np.sin(z)
            f[2] += theta*np.cos(a)*np.sin(z)
        return self.project(np.fft.fft2(f))

    @staticmethod
    def energy(h):
        f = np.fft.ifft2(h).real
        return np.array([.5*np.mean(f[0]**2+f[1]**2), .5*np.mean(f[2]**2)])

    def rhs(self, h, diffusion=False):
        f = np.fft.ifft2(h).real
        if self.nonlinear:
            dx = np.fft.ifft2(1j*self.kx*h).real
            dz = np.fft.ifft2(1j*self.kz*h).real
            r = -np.fft.fft2(f[0]*dx+f[1]*dz)
        else:
            r = np.zeros_like(h)
        coupling = self.pr*self.ra if self.buoyancy else 0.
        r[1] += coupling*h[2]
        r[2] += h[1]
        if diffusion:
            r -= np.array([self.pr, self.pr, 1.])[:, None, None]*self.k2*h
        flux = np.mean(f[1]*f[2])
        return self.project(r), np.array([coupling*flux, flux])

    def step(self, h, dt):
        # Exact diffusion half steps, with RK4 for projected nonlinear transport
        # and thermal/buoyancy coupling. Generally second-order Strang splitting.
        decay = np.exp(-dt/2*np.array([self.pr, self.pr, 1.])[:, None, None]*self.k2)
        e0 = self.energy(h)
        h = h*decay
        dissipated = e0-self.energy(h)
        a, wa = self.rhs(h)
        b, wb = self.rhs(h+dt*a/2)
        c, wc = self.rhs(h+dt*b/2)
        d, wd = self.rhs(h+dt*c)
        h = self.project(h+dt*(a+2*b+2*c+d)/6)
        work = dt*(wa+2*wb+2*wc+wd)/6
        e1 = self.energy(h)
        h = h*decay
        dissipated += e1-self.energy(h)
        return h, work, dissipated

    def diagnostics(self, h):
        f = np.fft.ifft2(h).real
        div = np.fft.ifft2(1j*(self.kx*h[0]+self.kz*h[1])).real
        uz = np.fft.ifft2(1j*self.kz*h[0]).real
        wall = max(np.max(abs(f[1:, :, [0, self.n//2]])), np.max(abs(uz[:, [0, self.n//2]])))
        speed = np.sqrt(np.mean(f[0]**2+f[1]**2))
        # A cos(kx+phase) sin(pi*z) has |Fourier[1,1]| = A*n^2/4.
        mode = 4*abs(h[1, 1, 1])/self.n**2
        return mode, speed, float(np.max(abs(div))), float(wall), f


COLUMNS = ['t', 'critical_w_amplitude', 'rms_velocity', 'kinetic_energy', 'thermal_variance_half',
           'buoyancy_work', 'kinetic_dissipation', 'thermal_production', 'thermal_dissipation',
           'kinetic_budget_residual', 'thermal_budget_residual', 'max_divergence', 'max_wall_residual']


def run(n=24, dt=.01, ratio=1., seed=410, amplitude=1e-8, duration=6., buoyancy=True, nonlinear=True, thermal_initial_multiplier=1.):
    solver = Layer(n, ratio*RAC, buoyancy=buoyancy, nonlinear=nonlinear)
    state = solver.initial(seed, amplitude)
    state[2] *= thermal_initial_multiplier
    count = round(duration/dt)
    assert abs(count*dt-duration) < 1e-10
    stride = max(1, round(PLAN['record_interval']/dt))
    work = np.zeros(2); diss = np.zeros(2)
    initial = solver.energy(state)
    rows, maps, map_times = [], [], []
    for i in range(count+1):
        if i % stride == 0 or i == count:
            mode, speed, div, wall, f = solver.diagnostics(state)
            energy = solver.energy(state)
            residual = energy-initial-work+diss
            rows.append([i*dt, mode, speed, *energy, work[0], diss[0], work[1], diss[1], *residual, div, wall])
        if i in (0, count//2, count):
            maps.append(np.fft.ifft2(state).real)
            map_times.append(i*dt)
        if i == count:
            break
        state, injected, dissipated = solver.step(state, dt)
        if not np.isfinite(state).all():
            raise FloatingPointError(f'Nonfinite PDE state at time {(i+1)*dt}; run rejected.')
        work += injected
        diss += dissipated
    rows = np.asarray(rows)
    fit = (rows[:, 0] >= 1) & (rows[:, 0] <= duration) & (rows[:, 1] > 0)
    rate = float(np.polyfit(rows[fit, 0], np.log(rows[fit, 1]), 1)[0]) if fit.sum() >= 3 else None
    scale = max(float(rows[:, 2].max()), 1e-300)
    expected = leading_rate(ratio*RAC) if buoyancy else -Q
    result = dict(n=n, dt=dt, ratio=ratio, seed=seed, amplitude=amplitude, duration=duration,
                  buoyancy=buoyancy, nonlinear=nonlinear, fitted_rate=rate, analytic_linear_rate=expected,
                  thermal_initial_multiplier=thermal_initial_multiplier,
                  mode_amplification=float(rows[-1, 1]/rows[0, 1]) if amplitude else 0.,
                  rate_error=abs(rate-expected) if rate is not None else None,
                  relative_divergence=float(rows[:, 11].max()/scale), relative_wall=float(rows[:, 12].max()/scale),
                  relative_kinetic_budget=float(np.max(abs(rows[:, 9]))/max(rows[:, 3].max(), 1e-300)),
                  relative_thermal_budget=float(np.max(abs(rows[:, 10]))/max(rows[:, 4].max(), 1e-300)))
    return result, rows, np.asarray(maps), np.asarray(map_times)
