"""FLUID. Four source cases. This file steps Navier-Stokes.

Produces fluid-sources.json.
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


def wmax(u, ops):
    kx, ky, kz, k2, keep, dx = ops
    uh = [m.spec(c, keep) for c in u]
    grads = [tuple(np.fft.ifftn(1j * k * comp).real for k in (kx, ky, kz)) for comp in uh]
    mag = np.sqrt((grads[2][1] - grads[1][2]) ** 2 + (grads[0][2] - grads[2][0]) ** 2 + (grads[1][0] - grads[0][1]) ** 2)
    return float(mag.max())


def main():
    n = 33
    x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
    ops = (kx, ky, kz, k2, keep, dx)
    psi = np.exp(-12 * ((x - 0.9) ** 2 + y * y)) * (x > 0)
    ab = [m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)]
    ab = [(a + b) / 2 for a, b in zip(ab, mirror_of(ab))]
    ab = m.project(ab, kx, ky, kz, k2, keep)
    ab = [(a + b) / 2 for a, b in zip(ab, mirror_of(ab))]
    blob = np.exp(-8 * ((x - 0.55) ** 2 + (y - 0.4) ** 2 + z * z))
    c = m.project([np.zeros_like(blob), 0.25 * blob, np.zeros_like(blob)], kx, ky, kz, k2, keep)
    c_sym = [(a + b) / 2 for a, b in zip(c, mirror_of(c))]
    c_uneven = [0.7 * a + 0.3 * b for a, b in zip(c, mirror_of(c))]
    phi = [(a - b) / 2 for a, b in zip(c, mirror_of(c))]
    phi = [v / np.sqrt(max(dot(phi, phi), 1e-30)) for v in phi]

    def dab(u):
        return dot([a - b for a, b in zip(u, mirror_of(u))], phi)

    def run(u0):
        u = m.project([v.copy() for v in u0], kx, ky, kz, k2, keep)
        p0 = float(np.sqrt(sum(v * v for v in u)).max())
        dt = 0.05 * dx / max(p0, 1e-6)
        T = 0.5
        steps = int(np.ceil(T / dt))
        dt = T / steps
        series = []
        every = max(1, steps // 5)
        for step in range(steps + 1):
            if step:
                u = m.advance(u, dt, ops)
            if step % every == 0 or step == steps:
                series.append({"t": round(step * dt, 3), "E": em(u), "D": dab(u), "W": wmax(u, ops)})
        return series

    fields = {
        "ab": ab,
        "ab_cd": [a + b for a, b in zip(ab, c_sym)],
        "ab_c": [a + b for a, b in zip(ab, c)],
        "ab_uneven": [a + b for a, b in zip(ab, c_uneven)],
        "ab_c_mirror": [a + b for a, b in zip(ab, mirror_of(c))],
    }
    out = {"kind": "fluid", "n": n, "T": 0.5}
    for name, field in fields.items():
        out[name] = run(field)
        print(name, out[name][-1], flush=True)
    Path("fluid-sources.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
