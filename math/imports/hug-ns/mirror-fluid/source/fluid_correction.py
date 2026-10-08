"""FLUID. Selector driven by the field difference. Steps Navier-Stokes.

D = <u - M[u], phi>
q is updated from D after each fluid step: q <- q + dt * (-0.2*q + 0.8*D)
The loop runs to time 0.16. That is the window in fluid-correction.json.
q is not a fluid force. It is a record of the field difference.
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

    def run(amp):
        u = m.project([a + amp * scale * b for a, b in zip(u_sym, phi)], kx, ky, kz, k2, keep)
        p0 = float(np.sqrt(sum(c * c for c in u)).max())
        dt = 0.04 * dx / max(p0, 1e-6)
        steps = max(1, int(np.ceil(0.16 / dt)))
        dt = 0.16 / steps
        q = 0.0
        for _ in range(steps):
            u = m.advance(u, dt, ops)
            q = q + dt * (-0.2 * q + 0.8 * dab(u))
        return {"amp": amp, "E": em(u), "D": dab(u), "q": q}

    rows = [run(a) for a in (0.0, 0.001, -0.001, 0.005, -0.005)]
    print(rows, flush=True)
    Path("fluid-correction.json").write_text(json.dumps({"kind": "fluid", "rows": rows}))


if __name__ == "__main__":
    main()
