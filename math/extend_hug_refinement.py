"""Extend the saved hug grid study without repeating completed runs.

--half-step adds the 256^3 time-step control while retaining previous results.
--finer-grid 384 adds a finer comparison at the saved longer time, t=0.16.
--time-control-grid 384 halves the time step at that same grid and final time.
--continue-to 0.24 resumes the saved 256/384 grid and 384 time-step controls.
--from-study selects a later published continuation as the exact restart baseline.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import time

import numpy as np

from hug_refinement import compare_fields, run
from continue_hug_refinement import continue_run, numerical_source_hash
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


def run_half_step_sequence(full, full_path, cache, stops):
    """Check the initial state and run half-sized steps through saved checkpoints."""
    points, dt = full['N'], full['dt']/2
    if (full['source_sha256'] != numerical_source_hash()
            or full['sigma'] != 0 or full['P_U'] != 0):
        raise ValueError('Require the unchanged unforced numerical source')
    if not stops or stops[-1] != full['end_time']:
        raise ValueError('Checkpoint sequence must end at the baseline final time')
    if any(b <= a for a, b in zip((0.,)+tuple(stops[:-1]), stops)):
        raise ValueError('Checkpoint times must increase from zero')
    total = round(stops[-1]/dt)
    original = PurePythonNavierStokes3D.step_with_pressure
    progress_count, phase_start = 0, 0
    began = time.perf_counter()

    def progress(self, *args, **kwargs):
        nonlocal progress_count
        if phase_start == 0 and progress_count == 0:
            with np.load(full_path) as reference:
                initial = reference['initial']
            if not all(np.array_equal(a, b) for a, b in
                       zip((self.u, self.v, self.w), initial)):
                raise ValueError('Time-step control initial arrays differ from the baseline')
            print('Initial arrays exactly match the saved full-step field', flush=True)
        result = original(self, *args, **kwargs)
        progress_count += 1
        print(f'N={points} dt={dt:g}: global step {phase_start+progress_count}/{total}; '
              f'elapsed {(time.perf_counter()-began)/60:.1f} min', flush=True)
        return result

    PurePythonNavierStokes3D.step_with_pressure = progress
    try:
        half, half_path = run(points, dt, stops[0], cache, box=full['box'], nu=full['nu'])
        for stop in stops[1:]:
            phase_start, progress_count = half['steps'], 0
            print(f'Checkpoint saved through t={half["end_time"]:g}; continuing to {stop:g}',
                  flush=True)
            half, half_path = continue_run(half_path.with_suffix('.json'), stop, cache)
    finally:
        PurePythonNavierStokes3D.step_with_pressure = original
    if half['source_sha256'] != full['source_sha256']:
        raise ValueError('Time-step control numerical sources differ')
    with np.load(full_path) as a, np.load(half_path) as b:
        if not np.array_equal(a['initial'], b['initial']):
            raise ValueError('Time-step control initial arrays differ from the baseline')
    return half, half_path


def extend_time_control(points):
    """Reuse a completed finer-grid baseline and measure time-step sensitivity."""
    previous = json.loads(Path(f'math/results/hug-grid-{points}.json').read_text())
    cache, end, dt = Path('scratch/hug-refinement'), .16, .001
    full_path = cache/f'n{points}-dt{dt:g}-t{end:g}.npz'
    full = json.loads(full_path.with_suffix('.json').read_text())
    expected = next(r for r in previous['runs'] if r['N'] == points and r['dt'] == dt)
    if full != expected or full['end_time'] != end or not full_path.exists():
        raise ValueError('Saved baseline must match the published finer-grid study')
    stops = (.04, .08, .12, .16)
    print(f'Reusing completed N={points}, dt={dt:g}, t={end:g} baseline', flush=True)
    half, half_path = run_half_step_sequence(full, full_path, cache, stops)
    print('Comparing complete fields on the identical grid', flush=True)
    with np.load(full_path) as a, np.load(half_path) as b:
        comparison = compare_fields(a['final'], b['final'], box=full['box'])
    result = dict(field=previous['field'], domain=previous['domain'], external_force=0.,
                  method=previous['method'], end_time=end, runs=[full, half],
                  time_comparison=dict(grid=points, steps=[dt, dt/2], **comparison),
                  initial_arrays_identical=True, checkpoint_times=list(stops),
                  earlier_time_comparison=previous['earlier_time_comparison'],
                  spatial_comparison=previous['spatial_comparison'],
                  comparison='Same physical grid; no interpolation in the time-step comparison',
                  limits=previous['limits']+' Two time steps measure sensitivity; '
                         'they do not establish a temporal convergence order.')
    path = Path(f'math/results/hug-time-{points}.json')
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with path.with_suffix('.csv').open('w', newline='') as handle:
        rows = [row for saved in result['runs'] for row in saved['rows']]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result['time_comparison'], indent=2), flush=True)


def resume_through_checkpoints(summary_path, stops, cache):
    """Continue only the new interval, retaining exact restart and cache checks."""
    summary_path = Path(summary_path)
    saved = json.loads(summary_path.read_text())
    dt, start = saved['dt'], saved['end_time']
    if (not stops or not math.isfinite(dt) or dt <= 0
            or any(not math.isfinite(t) for t in stops)
            or any(b <= a for a, b in zip((start,)+tuple(stops[:-1]), stops))
            or any(abs(round(t/dt)*dt-t) > 1e-12 for t in stops)):
        raise ValueError('Checkpoint times must increase after the saved time and align with dt')
    original = PurePythonNavierStokes3D.step_with_pressure
    current_step = saved['steps']
    total = round(stops[-1]/dt)
    began = time.perf_counter()

    def progress(self, *args, **kwargs):
        nonlocal current_step
        result = original(self, *args, **kwargs)
        current_step += 1
        print(f'N={self.N} dt={dt:g}: global step {current_step}/{total}; '
              f'elapsed {(time.perf_counter()-began)/60:.1f} min', flush=True)
        return result

    PurePythonNavierStokes3D.step_with_pressure = progress
    try:
        for stop in stops:
            current_step = saved['steps']
            saved, arrays_path = continue_run(summary_path, stop, cache)
            summary_path = arrays_path.with_suffix('.json')
            print(f'Saved N={saved["N"]} dt={dt:g} through t={stop:g}', flush=True)
    finally:
        PurePythonNavierStokes3D.step_with_pressure = original
    return saved, arrays_path


def continuation_baseline(baseline_path=None):
    """Select a recorded three-run continuation, or the original 0.16 controls."""
    if baseline_path is not None:
        previous = json.loads(Path(baseline_path).read_text())
        expected = previous['runs']
        if [(r['N'], r['dt']) for r in expected] != [(256, .001), (384, .001), (384, .0005)]:
            raise ValueError('Require the published 256/384 grid and 384 time-step controls')
        start = previous['end_time']
        if any(r['end_time'] != start for r in expected) or previous['external_force'] != 0:
            raise ValueError('Baseline runs must share the unforced endpoint')
        return previous, expected, previous['spatial_comparison']
    previous = json.loads(Path('math/results/hug-time-384.json').read_text())
    grid = json.loads(Path('math/results/hug-grid-384.json').read_text())
    if previous['end_time'] != grid['end_time']:
        raise ValueError('Baseline studies must share the endpoint')
    expected = [next(r for r in grid['runs'] if r['N'] == 256), *previous['runs']]
    return previous, expected, previous['spatial_comparison']['final']


def continue_grid_time_study(end, baseline_path=None):
    """Extend the selected published controls with unchanged grids and time steps."""
    previous, expected, earlier_spatial = continuation_baseline(baseline_path)
    start, cache = previous['end_time'], Path('scratch/hug-refinement')
    if not math.isfinite(end) or end <= start:
        raise ValueError('End must exceed the matching published checkpoint times')
    stops, t = [], start
    while t < end:
        t = min(round(t+.04, 12), end)
        stops.append(t)
    if [(r['N'], r['dt']) for r in expected] != [(256, .001), (384, .001), (384, .0005)]:
        raise ValueError('Require the published 256/384 grid and 384 time-step controls')
    source, paths = numerical_source_hash(), []
    for saved in expected:
        path = cache/f'n{saved["N"]}-dt{saved["dt"]:g}-t{start:g}.json'
        if (json.loads(path.read_text()) != saved or not path.with_suffix('.npz').exists()
                or saved['source_sha256'] != source or saved['sigma'] != 0
                or saved['P_U'] != 0 or saved['end_time'] != start
                or saved['box'] != expected[0]['box'] or saved['nu'] != expected[0]['nu']
                or any(abs(round(t/saved['dt'])*saved['dt']-t) > 1e-12 for t in stops)):
            raise ValueError('Checkpoint must exactly match the published unforced baseline')
        paths.append(path)
    print('All three checkpoints match the published results and unchanged numerical source',
          flush=True)
    with np.load(paths[1].with_suffix('.npz')) as a, np.load(paths[2].with_suffix('.npz')) as b:
        if not np.array_equal(a['initial'], b['initial']):
            raise ValueError('Time-step controls have different initial arrays')
    print('384-grid time-step controls have exactly identical original starting arrays', flush=True)
    runs, arrays = [], []
    for path in paths:
        saved, output = resume_through_checkpoints(path, stops, cache)
        runs.append(saved)
        arrays.append(output)
    print('Comparing complete fields at the new endpoint', flush=True)
    with np.load(arrays[0]) as low, np.load(arrays[1]) as full:
        spatial = compare_fields(low['final'], full['final'], box=expected[0]['box'])
    with np.load(arrays[1]) as full, np.load(arrays[2]) as half:
        if not np.array_equal(full['initial'], half['initial']):
            raise ValueError('Time-step controls have different initial arrays')
        temporal = compare_fields(full['final'], half['final'], box=expected[0]['box'])
    result = dict(field=previous['field'], domain=previous['domain'], external_force=0.,
                  method=previous['method'], start_time=start, end_time=end, runs=runs,
                  checkpoint_times=stops, initial_arrays_identical=True,
                  spatial_comparison=dict(coarse=256, fine=384, dt=.001, **spatial),
                  time_comparison=dict(grid=384, steps=[.001, .0005], **temporal),
                  earlier_comparisons=dict(time=start,
                      spatial=earlier_spatial,
                      temporal=previous['time_comparison']),
                  comparison='Grid comparison uses periodic cubic interpolation at matching '
                             'physical cell centers; time-step comparison uses the identical grid.',
                  limits=previous['limits'])
    path = Path(f'math/results/hug-continued-{end:g}.json')
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with path.with_suffix('.csv').open('w', newline='') as handle:
        rows = [row for saved in runs for row in saved['rows']]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(dict(spatial=spatial, time=temporal), indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--half-step", action="store_true")
    mode.add_argument("--finer-grid", type=int)
    mode.add_argument("--time-control-grid", type=int)
    mode.add_argument("--continue-to", type=float)
    parser.add_argument("--from-study", type=Path,
                        help="Published continuation JSON to resume instead of the original 0.16 controls")
    args = parser.parse_args()
    if args.from_study is not None and args.continue_to is None:
        parser.error('--from-study requires --continue-to')
    if args.continue_to is not None:
        continue_grid_time_study(args.continue_to, args.from_study)
        return
    if args.time_control_grid is not None:
        extend_time_control(args.time_control_grid)
        return
    if args.finer_grid is not None:
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
