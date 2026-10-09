"""The hug as specified.

psi_h = W * psi_0, W a flat C-infinity gate, exactly 0 near the faces.
Velocity from the product rule, then pressure-projected.
Unforced periodic Navier-Stokes, nu = 0.01, no body force.
"""

import json
from pathlib import Path

import numpy as np

L = 6.0
NU = 0.01
R = 1.5
EPS = 0.45
T_END = 0.40


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


def grids(n):
    dx = L / n
    axis = -L / 2 + dx * np.arange(n)
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


def energy(u):
    return float(0.5 * np.mean(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) * L**3)


def slope_max(u, kx, ky, kz, keep):
    peak = 0.0
    for c in u:
        for k in (kx, ky, kz):
            peak = max(peak, float(np.max(np.abs(deriv(c, k, keep)))))
    return peak


def checker(u):
    n = u[0].shape[0]
    s = ((np.arange(n)[:, None, None] + np.arange(n)[None, :, None] + np.arange(n)[None, None, :]) % 2) * 2.0 - 1.0
    num = sum(float(np.mean(c * s) ** 2) for c in u)
    den = sum(float(np.mean(c**2)) for c in u)
    return num / den


def gate(x, y, z):
    d = L / 2 - np.maximum(np.maximum(np.abs(x), np.abs(y)), np.abs(z))
    w = flat((d - EPS) / EPS)
    w[d <= EPS] = 0.0
    return w


def psi0(x, y, z):
    q = x * x + y * y
    return np.exp(-((q - R * R) ** 2) - z * z)


def build(n):
    x, y, z, kx, ky, kz, k2, keep, dx = grids(n)
    w = gate(x, y, z)
    p0 = psi0(x, y, z)
    # Product rule. Do not multiply the velocity by W.
    px, py, pz = deriv(p0, kx, keep), deriv(p0, ky, keep), deriv(p0, kz, keep)
    wx, wy, wz = deriv(w, kx, keep), deriv(w, ky, keep), deriv(w, kz, keep)
    u = [
        wy * p0 + w * py,
        -(wx * p0 + w * px),
        np.zeros_like(p0),
    ]
    collar = w == 0
    before = float(max(np.max(np.abs(c[collar])) for c in u))
    u = project(u, kx, ky, kz, k2, keep)
    after_collar = float(max(np.max(np.abs(c[collar])) for c in u))
    iz = n // 2
    speed = np.sqrt(u[0][:, :, iz] ** 2 + u[1][:, :, iz] ** 2 + u[2][:, :, iz] ** 2)
    report = {
        "n": n,
        "dx": dx,
        "divergence": div_max(u, kx, ky, kz, keep),
        "collarBeforeProjection": before,
        "collarAfterProjection": after_collar,
        "energy": energy(u),
        "slope": slope_max(u, kx, ky, kz, keep),
        "checker": checker(u),
        "gateMin": float(w.min()),
        "gateMax": float(w.max()),
        "support": float(np.mean(w > 0)),
    }
    return u, (kx, ky, kz, k2, keep, dx), report, speed


def rhs(u, kx, ky, kz, k2, keep):
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
    out_h = [-NU * k2 * uh[i] - ah[i] for i in range(3)]
    div2 = kx * out_h[0] + ky * out_h[1] + kz * out_h[2]
    lam2 = div2 * fac
    out_h = [out_h[i] - (kx, ky, kz)[i] * lam2 for i in range(3)]
    return [np.fft.ifftn(c).real for c in out_h]


def advance(u, dt, ops):
    kx, ky, kz, k2, keep, dx = ops
    r1 = rhs(u, kx, ky, kz, k2, keep)
    y = [u[i] + dt * r1[i] for i in range(3)]
    r2 = rhs(y, kx, ky, kz, k2, keep)
    z = [(u[i] + y[i] + dt * r2[i]) / 2 for i in range(3)]
    return project(z, kx, ky, kz, k2, keep)


def run(n):
    u, ops, report, speed0 = build(n)
    kx, ky, kz, k2, keep, dx = ops
    umax = float(np.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2).max())
    dt = min(0.2 * dx / max(umax, 1e-6), 0.1 * dx * dx / (6 * NU))
    steps = int(np.ceil(T_END / dt))
    dt = T_END / steps
    e0 = report["energy"]
    s0 = report["slope"]
    for step in range(1, steps + 1):
        u = advance(u, dt, ops)
        if step % 5 == 0:
            print(f"  n={n} t={step * dt:.3f}", flush=True)
    iz = n // 2
    speed1 = np.sqrt(u[0][:, :, iz] ** 2 + u[1][:, :, iz] ** 2 + u[2][:, :, iz] ** 2)
    report.update(
        {
            "dt": dt,
            "steps": steps,
            "umax0": umax,
            "energy0": e0,
            "energy1": energy(u),
            "slope0": s0,
            "slope1": slope_max(u, kx, ky, kz, keep),
            "divergence1": div_max(u, kx, ky, kz, keep),
            "checker1": checker(u),
            "mid0": speed0.tolist(),
            "mid1": speed1.tolist(),
        }
    )
    return report


def main():
    print("coarse", flush=True)
    a = run(32)
    print("fine", flush=True)
    b = run(48)
    # Compare midplane speeds on the coarse grid by sampling every fine point nearest? 
    # Compare energies and slopes, and midplane L2 after resizing fine by stride.
    fine = np.array(b["mid1"])
    # 48 to 32 is not integer. Compare divergence and energy only, plus a shared statistic.
    out = {
        "equation": "dt u + (u·grad)u = -grad p + 0.01 laplacian u, div u = 0",
        "construction": "psi_h = W * psi_0, product rule, then pressure projection. No force.",
        "nu": NU,
        "endTime": T_END,
        "coarse": {k: a[k] for k in a if k not in ("mid0", "mid1")},
        "fine": {k: b[k] for k in b if k not in ("mid0", "mid1")},
        "mid": b["mid1"],
        "midStart": b["mid0"],
    }
    Path("/workspace/src/lib/clay-hug.json").write_text(json.dumps(out))
    print(json.dumps({"coarse": out["coarse"], "fine": out["fine"]}, indent=2))


if __name__ == "__main__":
    main()
