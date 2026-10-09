"""One lean, two readings. Steps the ball and the field. Writes one-table.json."""
import json
from pathlib import Path
import importlib.util
import numpy as np

def ball(kind, x0=0.2, T=8.0, dt=0.02):
    x, v = x0, 0.0
    series = [(0.0, x)]
    n = int(T / dt)
    for i in range(1, n + 1):
        pull = -x if kind == "one" else x - x ** 3
        v += dt * (pull - 0.45 * v)
        x += dt * v
        if i % (n // 4) == 0:
            series.append((round(i * dt, 2), x))
    return series

def field(amp, T=0.4):
    spec = importlib.util.spec_from_file_location("clay_hug", Path(__file__).with_name("clay_hug.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.NU = 0.01

    def mirror_of(u):
        n = u[0].shape[0]
        ref = [(n - i) % n for i in range(n)]
        return [-u[0][ref], u[1][ref], u[2][ref]]

    def em(u):
        mu = mirror_of(u)
        num = sum(np.mean((a - b) ** 2) for a, b in zip(u, mu))
        den = sum(np.mean(a ** 2) for a in u)
        return float(np.sqrt(num / max(den, 1e-30)))

    n = 33
    x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
    ops = (kx, ky, kz, k2, keep, dx)
    psi = np.exp(-12 * ((x - 0.9) ** 2 + y * y)) * (x > 0)
    base = [m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)]
    base = [(a + b) / 2 for a, b in zip(base, mirror_of(base))]
    base = m.project(base, kx, ky, kz, k2, keep)
    base = [(a + b) / 2 for a, b in zip(base, mirror_of(base))]
    blob = np.exp(-8 * ((x - 0.55) ** 2 + (y - 0.4) ** 2 + z * z))
    c = m.project([np.zeros_like(blob), 0.25 * blob, np.zeros_like(blob)], kx, ky, kz, k2, keep)
    phi = [(a - b) / 2 for a, b in zip(c, mirror_of(c))]
    phi = [v / np.sqrt(max(sum(np.mean(w * w) for w in phi), 1e-30)) for v in phi]
    scale = float(np.sqrt(sum(np.mean(v * v) for v in base)))
    u = m.project([a + amp * scale * b for a, b in zip(base, phi)], kx, ky, kz, k2, keep)
    p0 = float(np.sqrt(sum(v * v for v in u)).max())
    dt = 0.04 * dx / max(p0, 1e-6)
    steps = int(np.ceil(T / dt))
    dt = T / steps
    out = [(0.0, em(u))]
    breath = np.exp(-8 * (x * x + y * y))
    push = m.project([breath * np.sign(x), np.zeros_like(breath), np.zeros_like(breath)], kx, ky, kz, k2, keep)
    for step in range(1, steps + 1):
        gap = 0.15 * np.sin(2 * np.pi * step * dt / 0.4)
        u = [a + gap * b for a, b in zip(u, push)]
        u = m.advance(u, dt, ops)
        if step in (steps // 2, steps):
            out.append((round(step * dt, 2), em(u)))
    return out

out = {
    "ball_one_low": ball("one"),
    "ball_two_low": ball("two"),
    "field_even": field(0.0),
    "field_lean": field(0.01),
}
Path(__file__).with_name("one-table.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
