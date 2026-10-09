"""Breathing hug rerun. This script wrote results.json.

Six starts. Gap opens and closes inside the fluid loop.
Grid 33 cubed. Viscosity 0.01. Time 0.4.
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

def dot(a, b):
    return float(sum(np.mean(x * y) for x, y in zip(a, b)))

def em(u):
    # E: unsigned size of the broken mirror
    mu = mirror_of(u)
    num = sum(np.mean((a - b) ** 2) for a, b in zip(u, mu))
    den = sum(np.mean(a ** 2) for a in u)
    return float(np.sqrt(num / max(den, 1e-30)))

n = 33
x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
ops = (kx, ky, kz, k2, keep, dx)

# even mirrored hug
psi = np.exp(-12 * ((x - 0.9) ** 2 + y * y)) * (x > 0)
base = [m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)]
base = [(a + b) / 2 for a, b in zip(base, mirror_of(base))]
base = m.project(base, kx, ky, kz, k2, keep)
base = [(a + b) / 2 for a, b in zip(base, mirror_of(base))]

# C blob and the antisymmetric template used for D
blob = np.exp(-8 * ((x - 0.55) ** 2 + (y - 0.4) ** 2 + z * z))
c = m.project([np.zeros_like(blob), 0.25 * blob, np.zeros_like(blob)], kx, ky, kz, k2, keep)
phi = [(a - b) / 2 for a, b in zip(c, mirror_of(c))]
phi = [v / np.sqrt(max(dot(phi, phi), 1e-30)) for v in phi]
scale = float(np.sqrt(max(dot(base, base), 1e-30)))

# opening-and-closing push, projected so it stays divergence-free
breath = np.exp(-8 * (x * x + y * y))
breath_u = m.project([breath * np.sign(x), np.zeros_like(breath), np.zeros_like(breath)], kx, ky, kz, k2, keep)

def dab(u):
    # D: signed side of the broken mirror
    return dot([a - b for a, b in zip(u, mirror_of(u))], phi)

def run(u0, T=0.4):
    u = m.project([v.copy() for v in u0], kx, ky, kz, k2, keep)
    p0 = float(np.sqrt(sum(c * c for c in u)).max())
    dt = 0.04 * dx / max(p0, 1e-6)
    steps = int(np.ceil(T / dt))
    dt = T / steps
    gaps = []
    for step in range(steps + 1):
        if step:
            gap = 0.15 * np.sin(2 * np.pi * step * dt / 0.4)
            u = [a + gap * b for a, b in zip(u, breath_u)]
            u = m.advance(u, dt, ops)
        if step in (0, steps):
            gaps.append({"t": round(step * dt, 3), "E": em(u), "D": dab(u)})
    return gaps

left = np.exp(-10 * ((x + 0.7) ** 2 + (y - 0.35) ** 2 + z * z))
right = np.exp(-10 * ((x - 0.7) ** 2 + (y + 0.35) ** 2 + z * z))
hug = (left - 0.85 * right) * (1 - 0.65 * np.exp(-30 * x * x))
hug_u = m.project(
    [m.deriv(hug, ky, keep), -m.deriv(hug, kx, keep), np.zeros_like(hug)],
    kx, ky, kz, k2, keep,
)
cases = {
    "even": base,
    "small lean": [a + 0.01 * scale * b for a, b in zip(base, phi)],
    "minus lean": [a - 0.01 * scale * b for a, b in zip(base, phi)],
    "extra mirror": [a + b for a, b in zip(base, [(p + q) / 2 for p, q in zip(c, mirror_of(c))])],
    "C no mirror": [a + b for a, b in zip(base, c)],
    "gap hug": hug_u,
}
out = {k: run(v) for k, v in cases.items()}
Path(__file__).with_name("results.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
