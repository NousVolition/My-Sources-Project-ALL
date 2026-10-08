"""Periodic unforced Navier-Stokes on a cube.

The hug is the starting flow only. A flat gate multiplies the stream
function, the velocity comes from the product rule, then the field is
projected so the divergence is zero. After time zero the envelope does
not act. There is no body force in the default step.

    d_t u + (u·grad) u = -grad p + nu * laplacian u
    div u = 0
"""

import numpy as np

L = 6.0
NU = 0.01
R = 1.5
EPS = 0.45


def bump(t):
    out = np.zeros_like(t, dtype=float)
    mask = t > 0
    out[mask] = np.exp(-1.0 / t[mask])
    return out


def flat(s):
    """0 for s <= 0, 1 for s >= 1, flat at both ends."""
    a = bump(s)
    b = bump(1.0 - s)
    out = np.zeros_like(s, dtype=float)
    live = (s > 0) & (s < 1)
    out[s >= 1] = 1.0
    denom = a[live] + b[live]
    out[live] = a[live] / denom
    return out


def grids(n, length=L):
    dx = length / n
    axis = -length / 2 + dx * np.arange(n)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    idx = np.fft.fftfreq(n) * n
    ix, iy, iz = np.meshgrid(idx, idx, idx, indexing="ij")
    keep = (np.abs(ix) <= n // 3) & (np.abs(iy) <= n // 3) & (np.abs(iz) <= n // 3)
    return x, y, z, kx, ky, kz, k2, keep, dx


def spec(f, keep):
    return np.fft.fftn(f) * keep


def deriv(f, k, keep):
    return np.fft.ifftn(1j * k * spec(f, keep)).real


def project(u, kx, ky, kz, k2, keep):
    uh = [spec(c, keep) for c in u]
    div = kx * uh[0] + ky * uh[1] + kz * uh[2]
    fac = np.zeros_like(k2)
    np.divide(1.0, k2, out=fac, where=k2 > 0)
    lam = div * fac
    return [np.fft.ifftn(uh[i] - (kx, ky, kz)[i] * lam).real for i in range(3)]


def div_max(u, kx, ky, kz, keep):
    d = deriv(u[0], kx, keep) + deriv(u[1], ky, keep) + deriv(u[2], kz, keep)
    return float(np.max(np.abs(d)))


def energy(u, length=L):
    return float(0.5 * np.mean(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) * length**3)


def slope_max(u, kx, ky, kz, keep):
    peak = 0.0
    for c in u:
        for k in (kx, ky, kz):
            peak = max(peak, float(np.max(np.abs(deriv(c, k, keep)))))
    return peak


def peak_speed(u):
    return float(np.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2).max())


def checker(u):
    n = u[0].shape[0]
    s = (
        (np.arange(n)[:, None, None] + np.arange(n)[None, :, None] + np.arange(n)[None, None, :])
        % 2
    ) * 2.0 - 1.0
    num = sum(float(np.mean(c * s) ** 2) for c in u)
    den = sum(float(np.mean(c**2)) for c in u)
    return num / den


def gate(x, y, z, eps=EPS, length=L):
    d = length / 2 - np.maximum(np.maximum(np.abs(x), np.abs(y)), np.abs(z))
    w = flat((d - eps) / eps)
    w[d <= eps] = 0.0
    return w


def ring(x, y, z, radius=R, sharp=1.0):
    q = x * x + y * y
    return np.exp(-sharp * ((q - radius * radius) ** 2 + z * z))


def pulled_tube(x, y, z):
    tube = np.exp(-18.0 * (x * x + y * y))
    strain = 0.35 * x * y * np.exp(-0.15 * (x * x + y * y + z * z))
    return tube + strain


def two_tubes(x, y, z, gap=0.55, sharp=12.0):
    left = np.exp(-sharp * ((x - gap) ** 2 + y * y + 0.04 * z * z))
    right = np.exp(-sharp * ((x + gap) ** 2 + y * y + 0.04 * z * z))
    return left - right


def build(n, psi=ring, nu=NU, eps=EPS):
    """Gated stream function, product-rule velocity, then projection."""
    x, y, z, kx, ky, kz, k2, keep, dx = grids(n)
    w = gate(x, y, z, eps=eps)
    p0 = psi(x, y, z)
    px, py = deriv(p0, kx, keep), deriv(p0, ky, keep)
    wx, wy = deriv(w, kx, keep), deriv(w, ky, keep)
    u = [wy * p0 + w * py, -(wx * p0 + w * px), np.zeros_like(p0)]
    collar = w == 0
    before = float(max(np.max(np.abs(c[collar])) for c in u)) if np.any(collar) else 0.0
    u = project(u, kx, ky, kz, k2, keep)
    after = float(max(np.max(np.abs(c[collar])) for c in u)) if np.any(collar) else 0.0
    ops = (kx, ky, kz, k2, keep, dx)
    report = {
        "n": n,
        "nu": nu,
        "divergence": div_max(u, kx, ky, kz, keep),
        "collarBeforeProjection": before,
        "collarAfterProjection": after,
        "energy": energy(u),
        "slope": slope_max(u, kx, ky, kz, keep),
        "peak": peak_speed(u),
        "checker": checker(u),
    }
    return u, ops, report


def rhs(u, ops, nu):
    kx, ky, kz, k2, keep, dx = ops
    uh = [spec(c, keep) for c in u]
    div = kx * uh[0] + ky * uh[1] + kz * uh[2]
    fac = np.zeros_like(k2)
    np.divide(1.0, k2, out=fac, where=k2 > 0)
    lam = div * fac
    uh = [uh[i] - (kx, ky, kz)[i] * lam for i in range(3)]
    grads = []
    for comp in uh:
        grads.append(tuple(np.fft.ifftn(1j * k * comp).real for k in (kx, ky, kz)))
    real = [np.fft.ifftn(c).real for c in uh]
    adv = [real[0] * grads[i][0] + real[1] * grads[i][1] + real[2] * grads[i][2] for i in range(3)]
    ah = [spec(a, keep) for a in adv]
    out_h = [-nu * k2 * uh[i] - ah[i] for i in range(3)]
    div2 = kx * out_h[0] + ky * out_h[1] + kz * out_h[2]
    lam2 = div2 * fac
    out_h = [out_h[i] - (kx, ky, kz)[i] * lam2 for i in range(3)]
    return [np.fft.ifftn(c).real for c in out_h]


def advance(u, dt, ops, nu):
    r1 = rhs(u, ops, nu)
    y = [u[i] + dt * r1[i] for i in range(3)]
    r2 = rhs(y, ops, nu)
    z = [(u[i] + y[i] + dt * r2[i]) / 2 for i in range(3)]
    kx, ky, kz, k2, keep, dx = ops
    return project(z, kx, ky, kz, k2, keep)


def timestep(u, ops, nu):
    dx = ops[-1]
    speed = max(peak_speed(u), 1e-6)
    return min(0.15 * dx / speed, 0.08 * dx * dx / (6 * nu))


def terms(u, ops, nu):
    """Biggest steepening and biggest smoothing. They may be in different spots."""
    kx, ky, kz, k2, keep, dx = ops
    uh = [spec(c, keep) for c in u]
    grads = []
    for comp in uh:
        grads.append(tuple(np.fft.ifftn(1j * k * comp).real for k in (kx, ky, kz)))
    real = [np.fft.ifftn(c).real for c in uh]
    adv = [real[0] * grads[i][0] + real[1] * grads[i][1] + real[2] * grads[i][2] for i in range(3)]
    visc = [np.fft.ifftn(-nu * k2 * uh[i]).real for i in range(3)]
    steep = float(np.sqrt(adv[0] ** 2 + adv[1] ** 2 + adv[2] ** 2).max())
    smooth = float(np.sqrt(visc[0] ** 2 + visc[1] ** 2 + visc[2] ** 2).max())
    return steep, smooth

def pressure_at_peak(u, ops, nu):
    """At one grid point, signed instantaneous rates of speed-squared/2.

    Use the projected state and filtered advection used by rhs. Positive means
    increasing local kinetic energy; negative means decreasing it. The selected
    fastest grid point can change between samples; it is not a tracked particle.
    The gate only shaped the start. See REVIEW.md for the diagnostic correction.
    """
    kx, ky, kz, k2, keep, dx = ops
    u = project(u, kx, ky, kz, k2, keep)
    uh = [spec(c, keep) for c in u]
    grads = []
    for comp in uh:
        grads.append(tuple(np.fft.ifftn(1j * k * comp).real for k in (kx, ky, kz)))
    real = [np.fft.ifftn(c).real for c in uh]
    adv = [real[0]*grads[i][0] + real[1]*grads[i][1] + real[2]*grads[i][2] for i in range(3)]
    ah = [spec(c, keep) for c in adv]
    adv = [np.fft.ifftn(c).real for c in ah]
    visc = [np.fft.ifftn(-nu * k2 * uh[i]).real for i in range(3)]
    dh = 1j * (kx*ah[0] + ky*ah[1] + kz*ah[2])
    ph = np.zeros_like(dh)
    # Delta p = -div(advection), so -k^2 p_hat = -div_hat.
    np.divide(dh, k2, out=ph, where=k2 > 0)
    gp = [np.fft.ifftn(1j * k * ph).real for k in (kx, ky, kz)]
    speed = np.sqrt(real[0]**2 + real[1]**2 + real[2]**2)
    ip = np.unravel_index(int(np.argmax(speed)), speed.shape)
    uu = [real[i][ip] for i in range(3)]
    def along(term):
        return float(sum(uu[i] * term[i][ip] for i in range(3)))
    carry, pressure, smoothing = -along(adv), -along(gp), along(visc)
    return {
        "peak": float(speed[ip]),
        "index": [int(i) for i in ip],
        "position": [float(-dx*speed.shape[axis]/2 + ip[axis]*dx) for axis in range(3)],
        "carry": carry,
        "pressure": pressure,
        "smoothing": smoothing,
        "total_local_energy_rate": carry + pressure + smoothing,
    }
