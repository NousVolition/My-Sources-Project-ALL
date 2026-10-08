"""FLUID. Small break to time 2. This file steps Navier-Stokes.

Produces the series in fluid-break-to-2.json / break-to-2.json.
"""
import json
from pathlib import Path
import numpy as np
import clay_hug as m

m.NU = 0.01


def mirror_of(u):
    n = u[0].shape[0]
    ref = [(n - i) % n for i in range(n)]
    return [-u[0][ref], u[1][ref], u[2][ref]]


def dot(a, b):
    return float(sum(np.mean(x * y) for x, y in zip(a, b)))


def em(u):
    mu = mirror_of(u)
    num = sum(np.mean((a - b) ** 2) for a, b in zip(u, mu))
    den = sum(np.mean(a ** 2) for a in u)
    return float(np.sqrt(num / max(den, 1e-30)))


def main():
    n = 33
    x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
    ops = (kx, ky, kz, k2, keep, dx)
    psi = np.exp(-12 * ((x - 0.9) ** 2 + y * y)) * (x > 0)
    u = [m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)]
    u = [(a + b) / 2 for a, b in zip(u, mirror_of(u))]
    u = m.project(u, kx, ky, kz, k2, keep)
    u = [(a + b) / 2 for a, b in zip(u, mirror_of(u))]
    blob = np.exp(-8 * ((x - 0.9) ** 2 + (y - 0.35) ** 2 + z * z)) * (x > 0)
    phi = [np.zeros_like(blob), blob, np.zeros_like(blob)]
    phi = m.project(phi, kx, ky, kz, k2, keep)
    phi = [(a - b) / 2 for a, b in zip(phi, mirror_of(phi))]
    phi = [c / np.sqrt(max(dot(phi, phi), 1e-30)) for c in phi]
    scale = float(np.sqrt(max(dot(u, u), 1e-30)))
    amp = 0.001
    u = m.project([a + amp * scale * b for a, b in zip(u, phi)], kx, ky, kz, k2, keep)

    def wmax(uu):
        uh = [m.spec(c, keep) for c in uu]
        grads = [tuple(np.fft.ifftn(1j * k * comp).real for k in (kx, ky, kz)) for comp in uh]
        mag = np.sqrt((grads[2][1] - grads[1][2]) ** 2 + (grads[0][2] - grads[2][0]) ** 2 + (grads[1][0] - grads[0][1]) ** 2)
        return float(mag.max())

    p0 = float(np.sqrt(sum(c * c for c in u)).max())
    dt = 0.05 * dx / max(p0, 1e-6)
    T = 2.0
    steps = int(np.ceil(T / dt))
    dt = T / steps
    series, acc, prev = [], 0.0, None
    every = max(1, steps // 8)
    for step in range(steps + 1):
        if step:
            u = m.advance(u, dt, ops)
        w = wmax(u)
        if prev is not None:
            acc += 0.5 * dt * (prev + w)
        prev = w
        if step % every == 0 or step == steps:
            series.append({"t": round(step * dt, 3), "E": em(u), "W": w, "I": acc})
            print(series[-1], flush=True)
    Path("fluid-break-to-2.json").write_text(json.dumps({"kind": "fluid", "series": series}))


if __name__ == "__main__":
    main()
