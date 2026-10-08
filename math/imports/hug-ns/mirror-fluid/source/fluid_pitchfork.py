"""FLUID. Two breaks, then the antisymmetric part is removed.

Produces fluid-pitchfork.json. Steps Navier-Stokes. Not the small model.
D is the projection of u - M[u] onto the antisymmetric template phi.
The break is removed by replacing u with (u + M[u]) / 2, then stepping again.
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
    u_sym = [m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)]
    u_sym = [(a + b) / 2 for a, b in zip(u_sym, mirror_of(u_sym))]
    u_sym = m.project(u_sym, kx, ky, kz, k2, keep)
    u_sym = [(a + b) / 2 for a, b in zip(u_sym, mirror_of(u_sym))]
    blob = np.exp(-8 * ((x - 0.9) ** 2 + (y - 0.35) ** 2 + z * z)) * (x > 0)
    phi = [np.zeros_like(blob), blob, np.zeros_like(blob)]
    phi = m.project(phi, kx, ky, kz, k2, keep)
    phi = [(a - b) / 2 for a, b in zip(phi, mirror_of(phi))]
    phi = [c / np.sqrt(max(dot(phi, phi), 1e-30)) for c in phi]
    scale = float(np.sqrt(max(dot(u_sym, u_sym), 1e-30)))

    def dab(u):
        return dot([a - b for a, b in zip(u, mirror_of(u))], phi)

    def advance(u, T):
        p0 = float(np.sqrt(sum(c * c for c in u)).max())
        dt = 0.04 * dx / max(p0, 1e-6)
        steps = max(1, int(np.ceil(T / dt)))
        dt = T / steps
        for _ in range(steps):
            u = m.advance(u, dt, ops)
        return u

    def strip(u):
        return [(a + b) / 2 for a, b in zip(u, mirror_of(u))]

    runs = {}
    for name, amp in (("plus", 0.005), ("minus", -0.005), ("zero", 0.0)):
        u = m.project([a + amp * scale * b for a, b in zip(u_sym, phi)], kx, ky, kz, k2, keep)
        u = advance(u, 0.12)
        held = {"E": em(u), "D": dab(u)}
        u2 = advance(strip(u), 0.12)
        runs[name] = {"before_strip": held, "after_strip": {"E": em(u2), "D": dab(u2)}}
        print(name, runs[name], flush=True)
    Path("fluid-pitchfork.json").write_text(json.dumps({"kind": "fluid", "runs": runs}))


if __name__ == "__main__":
    main()
