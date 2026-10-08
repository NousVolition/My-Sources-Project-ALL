"""Matched-stretch grid comparison.

This is the run behind grids-compared.json.
It did not save the full fields. Those were discarded after the counts.
Pass --save-fields to keep the velocity arrays.

The start is a tight tube in a stretch, projected so it stays divergence-free.
The gate in clay_hug is not used here. This start is built directly.
"""

import argparse
import json
from pathlib import Path

import numpy as np

import clay_hug as m

m.NU = 0.001


def build(n):
    x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
    envelope = np.exp(-0.04 * (x * x + y * y + z * z))
    strain = 40.0
    u = [-strain * x * envelope, -strain * y * envelope, 2 * strain * z * envelope]
    psi = 0.4 * np.exp(-25.0 * (x * x + y * y))
    u[0] = u[0] + m.deriv(psi, ky, keep)
    u[1] = u[1] - m.deriv(psi, kx, keep)
    return m.project(u, kx, ky, kz, k2, keep), (kx, ky, kz, k2, keep, dx)


def spin(u, ops):
    kx, ky, kz, k2, keep, dx = ops
    uh = [m.spec(c, keep) for c in u]
    grads = [tuple(np.fft.ifftn(1j * k * comp).real for k in (kx, ky, kz)) for comp in uh]
    wx = grads[2][1] - grads[1][2]
    wy = grads[0][2] - grads[2][0]
    wz = grads[1][0] - grads[0][1]
    return np.sqrt(wx * wx + wy * wy + wz * wz)


def measure(mag, dx, fixed=50.0):
    big = float(mag.max())
    half = mag >= 0.5 * big
    fixed_mask = mag >= fixed
    half_width = (int(half.sum()) * dx ** 3) ** (1 / 3)
    return {
        "biggest": big,
        "halfCells": int(half.sum()),
        "halfVolume": float(half.sum() * dx ** 3),
        "width": half_width,
        "dx": dx,
        "widthsPerCell": half_width / dx,
        "fixedCutoff": fixed,
        "fixedCells": int(fixed_mask.sum()),
        "fixedVolume": float(fixed_mask.sum() * dx ** 3),
    }


def run(n, t_end, samples, save_fields, outdir):
    u, ops = build(n)
    dx = ops[-1]
    p0 = float(np.sqrt(sum(c * c for c in u)).max())
    dt = 0.04 * dx / max(p0, 1e-6)
    steps = int(np.ceil(t_end / dt))
    dt = t_end / steps
    series = []
    acc = 0.0
    prev = None
    every = max(1, steps // samples)
    for step in range(steps + 1):
        if step:
            u = m.advance(u, dt, ops)
        mag = spin(u, ops)
        big = float(mag.max())
        if prev is not None:
            acc += 0.5 * dt * (prev + big)
        prev = big
        if step % every == 0 or step == steps:
            row = measure(mag, dx)
            row["t"] = round(step * dt, 4)
            row["integral"] = acc
            series.append(row)
            print(
                f"n={n} t={row['t']:.3f} biggest={big:.1f} "
                f"halfCells={row['halfCells']} fixedCells={row['fixedCells']}",
                flush=True,
            )
            if save_fields:
                path = outdir / f"field-n{n}-t{row['t']:.3f}.npz"
                np.savez_compressed(path, u0=u[0], u1=u[1], u2=u[2], spin=mag, dx=dx)
    return {"n": n, "dx": dx, "series": series}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-fields", action="store_true")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = {
        "what": "matched stretch, biggest spin and half-peak width",
        "note": "The published grids-compared.json used the half-peak cutoff only. A fixed cutoff of 50 is added here.",
        "n48": run(48, 0.20, 8, args.save_fields, outdir),
        "n64": run(64, 0.16, 8, args.save_fields, outdir),
        "n80": run(80, 0.12, 6, args.save_fields, outdir),
    }
    (outdir / "grids-compared.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
