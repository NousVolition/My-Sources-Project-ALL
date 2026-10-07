"""Extend the saved hug grid study without repeating completed runs.

--half-step adds the 256^3 time-step control while retaining previous results.
--finer-grid 384 adds a finer comparison at the saved longer time, t=0.16.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from hug_refinement import compare_fields, run
from continue_hug_refinement import continue_run
from navier import PurePythonNavierStokes3D


def extend_to_finer_grid(points):
    """Use the same dt on a finer grid and preserve intermediate checkpoints."""
    previous = json.loads(Path('math/results/hug-longer-time.json').read_text())
    end, dt, cache = previous['end_time'], .001, Path('scratch/hug-refinement')
    if points <= 256 or end != .16:
        raise ValueError('Require a grid finer than 256 and the saved time-0.16 study')
    low, low_path = run(256, dt, end, cache)
    expected = next(r for r in previous['runs'] if r['N'] == 256 and r['dt'] == dt)
    if low['source_sha256'] != expected['source_sha256'] or low['rows'] != expected['rows']:
        raise ValueError('Saved comparison does not match its cached 256-grid run')

    original = PurePythonNavierStokes3D.step_with_pressure
    progress_count = 0
    phase_start = 0
    began = time.perf_counter()
    def progress(self, *args, **kwargs):
        nonlocal progress_count
        result = original(self, *args, **kwargs)
        progress_count += 1
        print(f'N={self.N}: global step {phase_start+progress_count}/160; '
              f'phase step {progress_count}/40; elapsed {(time.perf_counter()-began)/60:.1f} min',
              flush=True)
        return result

    PurePythonNavierStokes3D.step_with_pressure = progress
    try:
        print(f'Preparing N={points}, dt={dt:g}, unforced ring; checkpoints every 0.04', flush=True)
        high, high_path = run(points, dt, .04, cache)
        for stop in (.08, .12, .16):
            phase_start, progress_count = high['steps'], 0
            print(f'Checkpoint saved through t={high["end_time"]:g}; continuing to {stop:g}', flush=True)
            high, high_path = continue_run(high_path.with_suffix('.json'), stop, cache)
    finally:
        PurePythonNavierStokes3D.step_with_pressure = original

    if high['source_sha256'] != low['source_sha256']:
        raise ValueError('Grid comparison uses different numerical source versions')
    print('Comparing complete fields at matching physical positions', flush=True)
    with np.load(low_path) as a, np.load(high_path) as b:
        comparison = dict(coarse=256, fine=points, dt=dt, time=end,
                          initial=compare_fields(a['initial'], b['initial']),
                          final=compare_fields(a['final'], b['final']))
    result = dict(field=previous['field'], domain=previous['domain'],
                  external_force=0., method=previous['method'], end_time=end,
                  runs=[low, high], spatial_comparison=comparison,
                  earlier_spatial_comparison=previous['spatial_comparison'],
                  earlier_time_comparison=previous['time_comparison'],
                  checkpoint_times=[.04,.08,.12,.16],
                  comparison='Periodic cubic interpolation at matching physical cell centers',
                  limits=previous['limits']+' The grid refinement ratios differ: '
                         f'128 to 256 is 2; 256 to {points} is {points/256:g}. Smaller differences alone '
                         'do not establish an observed convergence order.')
    path = Path(f'math/results/hug-grid-{points}.json')
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with path.with_suffix('.csv').open('w', newline='') as handle:
        rows = [row for saved in result['runs'] for row in saved['rows']]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(comparison, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--half-step", action="store_true")
    parser.add_argument("--finer-grid", type=int)
    args = parser.parse_args()
    if args.finer_grid is not None:
        if args.half_step:
            parser.error('Choose either the earlier half-step control or a finer-grid extension')
        extend_to_finer_grid(args.finer_grid)
        return
    path = Path("math/results/hug-refinement.json")
    study = json.loads(path.read_text())
    # Progress reporting does not modify the update or diagnostic calculation.
    original = PurePythonNavierStokes3D.step_with_pressure
    counters = {}
    def progress(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        key = id(self)
        counters[key] = counters.get(key, 0)+1
        if counters[key] % 5 == 0:
            dt = args[0] if args else kwargs['dt']
            print(f"N={self.N} dt={dt:g}: finished step {counters[key]}/{round(.08/dt)}", flush=True)
        return result
    PurePythonNavierStokes3D.step_with_pressure = progress
    try:
        low, low_path = run(128, .001, .08, "scratch/hug-refinement")
        high, high_path = run(256, .001, .08, "scratch/hug-refinement")
        if args.half_step:
            half, half_path = run(256, .0005, .08, "scratch/hug-refinement")
    finally:
        PurePythonNavierStokes3D.step_with_pressure = original
    # Retain the original reviewed results and time-step control. Matching
    # source fingerprints are checked by run() before any cache is reused.
    existing = next(r for r in study["runs"] if r["N"] == 128 and r["dt"] == .001)
    if existing["source_sha256"] != low["source_sha256"]:
        raise RuntimeError("Saved study and cached comparison use different source versions")
    print("Comparing complete fields at matching physical positions", flush=True)
    with np.load(low_path) as a, np.load(high_path) as b:
        comparison = dict(kind="grid", coarse=128, fine=256,
                          initial=compare_fields(a["initial"], b["initial"]),
                          final=compare_fields(a["final"], b["final"]))
    study["runs"] = [r for r in study["runs"] if (r["N"], r["dt"]) != (256, .001)]+[high]
    study["runs"].sort(key=lambda r: (r["N"], -r["dt"]))
    study["spatial_comparisons"] = [p for p in study["spatial_comparisons"]
                                    if (p["coarse"], p["fine"]) != (128, 256)]+[comparison]
    study["spatial_comparisons"].sort(key=lambda p: p["coarse"])
    time_checks = study.get("time_comparisons", [dict(
        grid=study.get("time_comparison_grid",128), steps=[.001,.0005],
        **study["time_comparison"])])
    if args.half_step:
        if half['source_sha256'] != high['source_sha256']:
            raise RuntimeError('Time comparison uses different source versions')
        with np.load(high_path) as a, np.load(half_path) as b:
            time_error = compare_fields(a['final'], b['final'])
            if not np.array_equal(a['initial'], b['initial']):
                raise RuntimeError('Time-step comparison must have identical initial arrays')
        study['runs'] = [r for r in study['runs'] if (r['N'],r['dt']) != (256,.0005)]+[half]
        study['runs'].sort(key=lambda r: (r['N'],-r['dt']))
        time_checks = [c for c in time_checks if c['grid'] != 256]+[
            dict(grid=256,steps=[.001,.0005],**time_error)]
        study['time_comparison'] = time_error
        study['time_comparison_grid'] = 256
        study['extension'] = '256^3 half-step completed; all earlier completed runs reused'
        print(json.dumps(dict(time_comparison_256=time_error),indent=2),flush=True)
    else:
        study.setdefault('time_comparison_grid',128)
        study.setdefault('extension','256^3 spatial run added; completed runs reused')
    study['time_comparisons'] = time_checks
    study['time_comparison_steps'] = [.001,.0005]
    path.write_text(json.dumps(study, indent=2, allow_nan=False)+"\n")
    rows = [row for r in study["runs"] for row in r["rows"]]
    with path.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(comparison, indent=2), flush=True)


if __name__ == "__main__":
    main()
