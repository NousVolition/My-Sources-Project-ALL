"""Compare the original and smooth starts on periodic cubes.

Default: original Gaussian, L=6, N=8,16,32, matching the submitted experiment.
--suite adds L=18 cases, a fixed-spacing control, the smooth field from a
common projected start, and one half-time-step check. NumPy is required for
the FFT projection and diagnostics; the original pure-Python step is used.
"""
import argparse
import copy
import csv
from functools import lru_cache
import math
from pathlib import Path

import numpy as np

from navier import PurePythonNavierStokes3D
from initial_field import SmoothRingField
from hugged_ring import HuggedRingField
from stokes import meridional_velocity, swirl

R0, ALPHA = 1.5, 1.0
# Match the old local radial width at the peak; this changes the full shape.
SMOOTH = SmoothRingField(radial_scale=math.sqrt(2 * R0 * ALPHA), axial_scale=ALPHA)


@lru_cache(maxsize=16)
def hug_for_box(box):
    return HuggedRingField(base=SMOOTH, box_length=box)


def field_velocity(profile, x, y, z, box=6.0):
    if profile == "smooth":
        return SMOOTH.velocity(x, y, z)
    if profile == "hug":
        return hug_for_box(box).velocity(x, y, z)
    if profile != "legacy":
        raise ValueError("profile must be legacy, smooth, or hug")
    r = math.hypot(x, y)
    ur, uz = meridional_velocity(r, z, R0, ALPHA)
    if r == 0.0:
        return 0.0, 0.0, uz
    uth = swirl(r, z, R0, ALPHA)
    return ur * x / r - uth * y / r, ur * y / r + uth * x / r, uz


def sample_stokes(sim, profile="legacy"):
    """Use sim.dx consistently; no radius clamp is needed."""
    mid = (sim.N - 1) / 2.0
    for i in range(sim.N):
        x = (i - mid) * sim.dx
        for j in range(sim.N):
            y = (j - mid) * sim.dx
            for k in range(sim.N):
                z = (k - mid) * sim.dx
                sim.u[i][j][k], sim.v[i][j][k], sim.w[i][j][k] = field_velocity(
                    profile, x, y, z, box=sim.N*sim.dx)
                sim.S[i][j][k] = 0.0


def velocity_array(sim):
    return np.array((sim.u, sim.v, sim.w))


def derivative(a, axis, dx):
    return (np.roll(a, -1, axis=axis) - np.roll(a, 1, axis=axis)) / (2 * dx)


def diagnostics(sim):
    velocity = velocity_array(sim)
    div = sum(derivative(velocity[c], c, sim.dx) for c in range(3))
    norm_grad2 = np.zeros_like(div)
    for c in range(3):
        for d in range(3):
            norm_grad2 += derivative(velocity[c], d, sim.dx)**2
    max_index = np.unravel_index(np.argmax(np.abs(div)), div.shape)
    return dict(
        max_div=float(np.max(np.abs(div))),
        interior_div=float(np.max(np.abs(div[1:-1, 1:-1, 1:-1]))),
        max_div_on_edge=int(any(i in (0, sim.N - 1) for i in max_index)),
        max_grad=float(np.sqrt(norm_grad2.max())),
        energy=float(0.5 * sim.dx**3 * np.sum(velocity**2)),
    )


def boundary_jump(profile, box, N):
    """Max sampled velocity mismatch at opposite physical faces."""
    coords = (np.arange(N) - (N - 1) / 2) * box / N
    worst = 0.0
    for axis in range(3):
        tangents = [d for d in range(3) if d != axis]
        for a in coords:
            for b in coords:
                plus, minus = [0.0] * 3, [0.0] * 3
                plus[axis], minus[axis] = box / 2, -box / 2
                for d, v in zip(tangents, (a, b)):
                    plus[d] = minus[d] = float(v)
                up, um = field_velocity(profile, *plus, box=box), field_velocity(profile, *minus, box=box)
                worst = max(worst, math.sqrt(sum((p - m)**2 for p, m in zip(up, um))))
    return worst


def run_case(profile, start, box, N, dt=0.02, steps=4):
    sim = PurePythonNavierStokes3D(N=N, dx=box / N)
    # Isolate fluid dynamics: tanh's tail otherwise leaves a tiny constant force.
    sim.sigma = 0.0
    sample_stokes(sim, profile)
    metadata = dict(profile=profile, start=start, box=box, N=N, dx=sim.dx, dt=dt,
                    face_jump=boundary_jump(profile, box, N))
    rows = [dict(metadata, phase="sampled", step=0, time=0.0, branch="both",
                 correction_fraction=0.0, **diagnostics(sim))]
    raw = velocity_array(sim)
    if start == "projected":
        sim.project(backend="fft")
    common = velocity_array(sim)
    correction = float(np.linalg.norm(common - raw) / np.linalg.norm(raw))
    rows.append(dict(metadata, phase="start", step=0, time=0.0, branch="both",
                     correction_fraction=correction, **diagnostics(sim)))
    bare, held = sim, copy.deepcopy(sim)
    print(f"{profile}/{start}  L={box:g} N={N} dx={sim.dx:g}"
          f"  sampled div={rows[0]['max_div']:.5g}"
          f"  face jump={metadata['face_jump']:.5g}", flush=True)
    for n in range(1, steps + 1):
        bare.step(dt=dt, P_U=0.0)
        held.step_with_pressure(dt=dt, P_U=0.0, backend="fft")
        for label, cube in (("no_projection", bare), ("projected", held)):
            rows.append(dict(metadata, phase="step", step=n, time=n*dt, branch=label,
                             correction_fraction=correction, **diagnostics(cube)))
        d0, dp = rows[-2], rows[-1]
        print(f"  step {n}: div {d0['max_div']:.5g} / {dp['max_div']:.3g}"
              f"  grad {d0['max_grad']:.5g} / {dp['max_grad']:.5g}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("legacy", "smooth", "hug"), default="legacy")
    parser.add_argument("--start", choices=("raw", "projected"), default="raw")
    parser.add_argument("--box", type=float, default=6.0)
    parser.add_argument("--points", nargs="+", type=int, default=[8, 16, 32])
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--suite", action="store_true")
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    if args.box <= 0 or args.dt <= 0 or args.steps < 1 or any(n < 3 for n in args.points):
        parser.error("Require box, dt > 0, steps >= 1 and points >= 3.")
    rows = []
    if args.suite:
        cases = [(6.0, 8), (6.0, 16), (6.0, 32), (18.0, 16), (18.0, 32), (18.0, 48)]
        for profile, start in (("legacy", "raw"), ("smooth", "projected")):
            for box, N in cases:
                rows.extend(run_case(profile, start, box, N))
        rows.extend(run_case("smooth", "projected", 6.0, 32, dt=0.01, steps=8))
    else:
        for N in args.points:
            rows.extend(run_case(args.profile, args.start, args.box, N, args.dt, args.steps))
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved {len(rows)} diagnostic rows to {args.csv}", flush=True)


if __name__ == "__main__":
    main()
