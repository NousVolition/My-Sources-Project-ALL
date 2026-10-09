"""Fluid dashboard the bowl does not cover.

Records E, D, W, I, energy, enstrophy, max/average spin,
fixed and half-peak cell counts, terms at the peak-spin cell,
divergence, and spectral energy near the cutoff.
"""
import json
from pathlib import Path
import importlib.util
import numpy as np

spec = importlib.util.spec_from_file_location("clay_hug", Path(__file__).with_name("clay_hug.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.NU = 0.01

def mirror_of(u):
    n = u[0].shape[0]
    ref = [(n - i) % n for i in range(n)]
    return [-u[0][ref], u[1][ref], u[2][ref]]

def curl(u, kx, ky, kz, keep):
    ux, uy, uz = u
    return [
        m.deriv(uy, kz, keep) - m.deriv(uz, ky, keep),
        m.deriv(uz, kx, keep) - m.deriv(ux, kz, keep),
        m.deriv(ux, ky, keep) - m.deriv(uy, kx, keep),
    ]

def snapshot(u, phi, ops, I):
    kx, ky, kz, k2, keep, dx = ops
    w = curl(u, kx, ky, kz, keep)
    mag = np.sqrt(sum(c * c for c in w))
    W = float(mag.max())
    loc = np.unravel_index(int(mag.argmax()), mag.shape)
    adv = [
        u[0] * m.deriv(u[i], kx, keep) + u[1] * m.deriv(u[i], ky, keep) + u[2] * m.deriv(u[i], kz, keep)
        for i in range(3)
    ]
    visc = [m.NU * (m.deriv(m.deriv(u[i], kx, keep), kx, keep) + m.deriv(m.deriv(u[i], ky, keep), ky, keep) + m.deriv(m.deriv(u[i], kz, keep), kz, keep)) for i in range(3)]
    stretch = [
        w[0] * m.deriv(u[i], kx, keep) + w[1] * m.deriv(u[i], ky, keep) + w[2] * m.deriv(u[i], kz, keep)
        for i in range(3)
    ]
    uh = m.spec(u[0], keep)
    shell = float(np.mean(np.abs(uh) ** 2))
    edge = float(np.mean(np.abs(uh * (k2 > 0.6 * k2.max())) ** 2))
    mu = mirror_of(u)
    E = float(np.sqrt(sum(np.mean((a - b) ** 2) for a, b in zip(u, mu)) / max(sum(np.mean(a ** 2) for a in u), 1e-30)))
    D = float(sum(np.mean((a - b) * p) for a, b, p in zip(u, mu, phi)))
    return {
        "E": E, "D": D, "W": W, "I": I,
        "energy": m.energy(u),
        "enstrophy": float(0.5 * np.mean(mag ** 2) * m.L ** 3),
        "ratio": W / max(float(mag.mean()), 1e-30),
        "cells_50": int(np.sum(mag >= 50)),
        "cells_100": int(np.sum(mag >= 100)),
        "half_peak": int(np.sum(mag >= 0.5 * W)),
        "at_peak": {
            "advection": float(np.sqrt(sum(c[loc] ** 2 for c in adv))),
            "viscosity": float(np.sqrt(sum(c[loc] ** 2 for c in visc))),
            "stretch": float(np.sqrt(sum(c[loc] ** 2 for c in stretch))),
        },
        "div": m.div_max(u, kx, ky, kz, keep),
        "cutoff_share": edge / max(shell, 1e-30),
    }

def main():
    n = 33
    x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
    ops = (kx, ky, kz, k2, keep, dx)
    psi = np.exp(-12 * ((x - 0.9) ** 2 + y * y)) * (x > 0)
    base = [m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)]
    base = [(a + b) / 2 for a, b in zip(base, mirror_of(base))]
    base = m.project(base, kx, ky, kz, k2, keep)
    blob = np.exp(-8 * ((x - 0.55) ** 2 + (y - 0.4) ** 2 + z * z))
    c = m.project([np.zeros_like(blob), 0.25 * blob, np.zeros_like(blob)], kx, ky, kz, k2, keep)
    phi = [(a - b) / 2 for a, b in zip(c, mirror_of(c))]
    phi = [v / np.sqrt(max(sum(np.mean(w * w) for w in phi), 1e-30)) for v in phi]
    scale = float(np.sqrt(sum(np.mean(v * v) for v in base)))
    out = {}
    for name, h in (("even", 0.0), ("plus", 0.01)):
        u = m.project([a + h * scale * b for a, b in zip(base, phi)], kx, ky, kz, k2, keep)
        p0 = float(np.sqrt(sum(v * v for v in u)).max())
        T = 0.4
        dt = 0.04 * dx / max(p0, 1e-6)
        steps = int(np.ceil(T / dt))
        dt = T / steps
        I = 0.0
        rows = []
        for step in range(steps + 1):
            if step:
                u = m.advance(u, dt, ops)
                I += dt * snapshot(u, phi, ops, I)["W"]
            if step in (0, steps):
                row = snapshot(u, phi, ops, I)
                row["t"] = round(step * dt, 3)
                rows.append(row)
                print(name, row["t"], row["W"], row["E"], flush=True)
        out[name] = rows
    Path(__file__).with_name("dashboard.json").write_text(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
