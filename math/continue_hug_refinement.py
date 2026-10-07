"""Continue saved unforced hug runs, without repeating the earlier interval.

The Euler update is one-step: velocity and the zero scalar are sufficient to
resume. The pressure potential is recomputed by the existing projection at
each step; the saved velocity must not be projected a second time on restart.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from box_experiment import diagnostics, velocity_array
from hug_refinement import compare_fields
from navier import PurePythonNavierStokes3D


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def numerical_source_hash():
    digest = hashlib.sha256()
    for name in ('navier.py', 'box_experiment.py', 'hugged_ring.py',
                 'hug_envelope.py', 'initial_field.py', 'hug_refinement.py'):
        digest.update(Path(__file__).with_name(name).read_bytes())
    return digest.hexdigest()


def continue_run(summary_path, end, output_cache):
    """Resume one compatible checkpoint and retain its complete time history."""
    summary_path = Path(summary_path)
    arrays_path = summary_path.with_suffix('.npz')
    saved = json.loads(summary_path.read_text())
    if saved['source_sha256'] != numerical_source_hash():
        raise ValueError('Checkpoint numerical source differs from current source')
    if saved['sigma'] != 0 or saved['P_U'] != 0:
        raise ValueError('Continuation requires an unforced checkpoint with zero scalar')
    n, dt, start = saved['N'], saved['dt'], saved['end_time']
    steps = round(end / dt)
    if dt <= 0 or end <= start or abs(steps * dt - end) > 1e-12:
        raise ValueError('End time must exceed checkpoint time and be a multiple of dt')
    first_step = saved['steps']
    if abs(first_step * dt - start) > 1e-12 or saved['rows'][-1]['step'] != first_step:
        raise ValueError('Checkpoint time history is inconsistent')
    provenance = dict(summary_sha256=file_hash(summary_path),
                      arrays_sha256=file_hash(arrays_path),
                      continuation_source_sha256=file_hash(__file__),
                      resumed_from_time=start,
                      resumed_from_file=summary_path.name,
                      state='exact saved velocity; scalar remains zero; no restart projection')
    output_cache = Path(output_cache)
    output_cache.mkdir(parents=True, exist_ok=True)
    key = f'n{n}-dt{dt:g}-t{end:g}'
    result_path = output_cache / (key + '.json')
    result_arrays = output_cache / (key + '.npz')
    if result_path.exists() and result_arrays.exists():
        previous = json.loads(result_path.read_text())
        if previous.get('continuation') == provenance and previous['end_time'] == end:
            print(f'Reusing completed continuation {key}', flush=True)
            return previous, result_arrays

    began = time.perf_counter()
    with np.load(arrays_path) as arrays:
        initial, checkpoint = arrays['initial'], arrays['final']
    if checkpoint.shape != (3, n, n, n) or initial.shape != checkpoint.shape:
        raise ValueError('Checkpoint array shape does not match its grid')
    if not np.isfinite(checkpoint).all():
        raise ValueError('Checkpoint contains non-finite velocity')
    sim = PurePythonNavierStokes3D(n, saved['box'] / n)
    sim.nu, sim.sigma = saved['nu'], 0.
    sim.u, sim.v, sim.w = checkpoint
    sim.S = np.zeros((n, n, n), dtype=float)
    resumed = diagnostics(sim)
    for name in ('energy', 'max_grad', 'max_div'):
        if not np.isclose(resumed[name], saved['rows'][-1][name], rtol=1e-12, atol=1e-14):
            raise ValueError(f'Checkpoint diagnostic mismatch: {name}')
    rows = list(saved['rows'])
    interval = max(1, (steps - first_step) // 4)
    print(f'Resuming N={n} dt={dt:g} from t={start:g} to {end:g}; '
          f'{steps-first_step} new steps', flush=True)
    for step in range(first_step + 1, steps + 1):
        sim.step_with_pressure(dt, 0., backend='fft', step_backend='numpy')
        added = step - first_step
        if added % 5 == 0 or step == steps:
            print(f'N={n} dt={dt:g}: added step {added}/{steps-first_step}', flush=True)
        if added % interval == 0 or step == steps:
            row = dict(N=n, dt=dt, step=step, time=step*dt, **diagnostics(sim))
            speed2 = sim.u**2 + sim.v**2 + sim.w**2
            row['max_speed'] = float(np.sqrt(speed2.max()))
            row['advective_cfl'] = float(dt/sim.dx*np.max(np.abs(sim.u)+np.abs(sim.v)+np.abs(sim.w)))
            row['diffusion_number'] = 6*sim.nu*dt/sim.dx**2
            row['linear_advection_diffusion_ratio'] = float(dt*speed2.max()/(2*sim.nu))
            rows.append(row)
            print(f"  t={row['time']:.3f} gradient={row['max_grad']:.6f} "
                  f"energy={row['energy']:.6f} div={row['max_div']:.2g}", flush=True)
    final = velocity_array(sim)
    if not np.isfinite(final).all() or np.any(sim.S != 0):
        raise RuntimeError('Non-finite velocity or nonzero scalar in unforced continuation')
    np.savez_compressed(result_arrays, initial=initial, final=final)
    result = dict(saved, steps=steps, end_time=end, rows=rows,
                  continuation=provenance, previous_seconds=saved['seconds'],
                  seconds=time.perf_counter()-began,
                  seconds_scope='new continuation only, including checkpoint restore and save')
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result, result_arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-time', type=float, default=.08)
    parser.add_argument('--end', type=float, default=.16)
    parser.add_argument('--cache', type=Path, default=Path('scratch/hug-refinement'))
    parser.add_argument('--out', type=Path, default=Path('math/results/hug-longer-time.json'))
    args = parser.parse_args()
    if args.end <= args.from_time:
        parser.error('End time must exceed checkpoint time')
    earlier = json.loads(Path('math/results/hug-refinement.json').read_text())
    if any(run['end_time'] != args.from_time for run in earlier['runs']):
        parser.error('Checkpoint time must match the saved refinement study')
    runs, paths = [], []
    for n, dt in ((128, .001), (256, .001), (256, .0005)):
        path = args.cache / f'n{n}-dt{dt:g}-t{args.from_time:g}.json'
        run, arrays = continue_run(path, args.end, args.cache)
        runs.append(run)
        paths.append(arrays)
    print('Comparing the complete longer-time fields', flush=True)
    with np.load(paths[0]) as low, np.load(paths[1]) as high:
        spatial = compare_fields(low['final'], high['final'])
    with np.load(paths[1]) as full, np.load(paths[2]) as half:
        if not np.array_equal(full['initial'], half['initial']):
            raise ValueError('Time-step comparison did not start from identical arrays')
        temporal = compare_fields(full['final'], half['final'])
    result = dict(field=earlier['field'], domain=earlier['domain'],
                  external_force=0., method=earlier['method'],
                  start_time=args.from_time, end_time=args.end, runs=runs,
                  spatial_comparison=dict(coarse=128, fine=256, dt=.001, **spatial),
                  time_comparison=dict(grid=256, steps=[.001, .0005], **temporal),
                  earlier_comparisons=dict(time=args.from_time,
                      spatial=earlier['spatial_comparisons'][-1]['final'],
                      temporal=earlier['time_comparison']),
                  limits=earlier['limits'])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with args.out.with_suffix('.csv').open('w', newline='') as handle:
        rows = [row for run in runs for row in run['rows']]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(dict(spatial=spatial, time=temporal), indent=2), flush=True)


if __name__ == '__main__':
    main()
