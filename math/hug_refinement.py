"""Fixed-time refinement of the existing hugged ring, using the cube stencils.

The NumPy backend changes evaluation speed, not the centered spatial stencil,
forward-Euler update, pressure projection, initial field, viscosity, or force.
Final arrays are cached under scratch/ so comparisons can be reproduced.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.ndimage import map_coordinates

from box_experiment import derivative, diagnostics, hug_for_box, sample_stokes, velocity_array
from navier import PurePythonNavierStokes3D


def sample_on_grid(values, points):
    """Periodic cubic interpolation onto the SAME physical cell centers.

    Different even grids are not nested. Merely taking every other point
    introduces a half-cell displacement. Order-3 interpolation is O(h^4) for
    sufficiently resolved smooth data; it is not an exact field comparison.
    """
    fine = values.shape[0]
    index = (np.arange(points)+.5)*fine/points-.5
    coords = np.asarray(np.meshgrid(index, index, index, indexing="ij"))
    return map_coordinates(values, coords, order=3, mode="grid-wrap", prefilter=True)


def compare_fields(coarse, fine, box=6.):
    nc, nf = coarse.shape[1], fine.shape[1]
    def transfer(a):
        return a if nc == nf else sample_on_grid(a, nc)
    sampled = np.array([transfer(a) for a in fine])
    velocity_error = float(np.linalg.norm(coarse-sampled)/np.linalg.norm(sampled))
    numerator, denominator = 0., 0.
    for c in range(3):
        for d in range(3):
            a = derivative(coarse[c], d, box/nc)
            b = transfer(derivative(fine[c], d, box/nf))
            numerator += float(np.sum((a-b)**2))
            denominator += float(np.sum(b*b))
    return dict(velocity_relative_l2=velocity_error,
                gradient_relative_l2=float(np.sqrt(numerator/denominator)))


def run(points, dt, end, cache, box=6., nu=.01):
    steps = round(end/dt)
    if steps < 1 or abs(steps*dt-end) > 1e-12:
        raise ValueError("End time must be a positive multiple of dt")
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    key = f"n{points}-dt{dt:g}-t{end:g}"
    summary_path = cache/(key+".json")
    arrays_path = cache/(key+".npz")
    digest = hashlib.sha256()
    for name in ("navier.py", "box_experiment.py", "hugged_ring.py", "hug_envelope.py",
                 "initial_field.py", "hug_refinement.py"):
        digest.update(Path(__file__).with_name(name).read_bytes())
    fingerprint = digest.hexdigest()
    if summary_path.exists() and arrays_path.exists():
        saved = json.loads(summary_path.read_text())
        if (saved.get("source_sha256") == fingerprint and saved["box"] == box
                and saved["nu"] == nu):
            print(f"Reusing completed {key}", flush=True)
            return saved, arrays_path
    began = time.perf_counter()
    sim = PurePythonNavierStokes3D(points, box/points)
    sim.nu, sim.sigma = nu, 0.
    sample_stokes(sim, "hug")
    raw = velocity_array(sim)
    sim.u, sim.v, sim.w = raw.copy()
    sim.S = np.asarray(sim.S)
    sampled = diagnostics(sim)
    sim.project(backend="fft")
    initial = velocity_array(sim)
    correction = float(np.linalg.norm(initial-raw)/np.linalg.norm(raw))
    del raw
    def measure(step):
        row = dict(N=points, dt=dt, step=step, time=step*dt, **diagnostics(sim))
        speed2 = sim.u**2+sim.v**2+sim.w**2
        row["max_speed"] = float(np.sqrt(speed2.max()))
        row["advective_cfl"] = float(dt/sim.dx*np.max(np.abs(sim.u)+np.abs(sim.v)+np.abs(sim.w)))
        row["diffusion_number"] = 6*nu*dt/sim.dx**2
        row["linear_advection_diffusion_ratio"] = float(dt*speed2.max()/(2*nu))
        return row
    rows = [measure(0)]
    print(f"N={points} dt={dt:g}: initial gradient {rows[0]['max_grad']:.6f}", flush=True)
    for step in range(1, steps+1):
        sim.step_with_pressure(dt, 0., backend="fft", step_backend="numpy")
        if step % max(1, steps//4) == 0 or step == steps:
            row = measure(step)
            rows.append(row)
            print(f"  N={points} dt={dt:g} t={row['time']:.3f}: "
                  f"gradient={row['max_grad']:.6f} energy={row['energy']:.6f} "
                  f"div={row['max_div']:.2g}", flush=True)
    final = velocity_array(sim)
    if not np.all(np.isfinite(final)) or np.any(sim.S != 0):
        raise RuntimeError("Non-finite velocity or nonzero scalar in unforced run")
    np.savez_compressed(arrays_path, initial=initial, final=final)
    result = dict(N=points, dt=dt, steps=steps, end_time=end, box=box, nu=nu,
                  sigma=0., P_U=0., source_sha256=fingerprint,
                  seconds=time.perf_counter()-began,
                  sampled=sampled, correction_fraction=correction, rows=rows)
    summary_path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result, arrays_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, nargs="+", default=[32, 64, 128])
    parser.add_argument("--dt", type=float, default=.001)
    parser.add_argument("--end", type=float, default=.08)
    parser.add_argument("--cache", type=Path, default=Path("scratch/hug-refinement"))
    parser.add_argument("--out", type=Path, default=Path("math/results/hug-refinement.json"))
    args = parser.parse_args()
    if args.dt <= 0 or args.end <= 0 or any(n < 4 for n in args.points):
        parser.error("Require positive dt/end and grids >= 4")
    points = sorted(set(args.points))
    runs, paths = [], []
    for n in points:
        saved, path = run(n, args.dt, args.end, args.cache)
        runs.append(saved)
        paths.append(path)
    half, half_path = run(points[-1], args.dt/2, args.end, args.cache)
    comparisons = []
    for left, right, n0, n1 in zip(paths[:-1], paths[1:], points[:-1], points[1:]):
        with np.load(left) as a, np.load(right) as b:
            comparisons.append(dict(kind="grid", coarse=n0, fine=n1,
                                    initial=compare_fields(a["initial"], b["initial"]),
                                    final=compare_fields(a["final"], b["final"])))
    with np.load(paths[-1]) as a, np.load(half_path) as b:
        time_error = compare_fields(a["final"], b["final"])
    result = dict(field=hug_for_box(6.).metadata(), domain="periodic cube [-3,3]^3",
                  external_force=0., method="original centered stencils, Euler, centered FFT projection",
                  verification="NumPy backend compared to original loop updates and Jacobi projection",
                  comparison="periodic cubic interpolation at matching physical cell centers",
                  runs=runs+[half], spatial_comparisons=comparisons, time_comparison=time_error,
                  limits="Finite-time numerical evidence; no error bound or convergence theorem. "
                         "Centered derivatives have checkerboard null modes. Interpolation is approximate.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    rows = [row for r in runs+[half] for row in r["rows"]]
    with args.out.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(dict(spatial=comparisons, time=time_error), indent=2), flush=True)


if __name__ == "__main__":
    main()
