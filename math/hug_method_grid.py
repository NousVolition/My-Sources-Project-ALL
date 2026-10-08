"""Repeat the longer method comparison on 256^3; reuse the saved 384^3 result.

Each grid's four controls share its own saved input. Prior histories differ
between grids, so these are within-grid sensitivity checks, not a pure
spatial-convergence comparison. All numerical stepping code is reused.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from continue_hug_refinement import file_hash, numerical_source_hash
from hug_method_interval import observe, run_branch, write_json
from hug_refinement import compare_fields
from hug_time_method import difference


def study(input_path, reference_path, checkpoint_cache, cache, out):
    input_path, reference_path, checkpoint_cache, cache, out = map(
        Path, (input_path, reference_path, checkpoint_cache, cache, out))
    source = json.loads(input_path.read_text())
    reference = json.loads(reference_path.read_text())
    fingerprint = numerical_source_hash()
    module = Path(__file__).parent
    runner_hash = file_hash(module / 'hug_method_interval.py')
    integrator_hash = file_hash(module / 'hug_time_method.py')
    configurations = (('euler', .001), ('euler', .0005),
                      ('heun', .001), ('heun', .0005))
    if (source['external_force'], source['end_time'], source['numerical_source_sha256']) != (
            0., .4, fingerprint):
        raise ValueError('Require the verified unforced saved-flow record at 0.40')
    if (reference['grid'], reference['start_time'], reference['end_time'],
            reference['box'], reference['nu'], reference['external_force'],
            reference['input_history_dt'], reference['numerical_source_sha256'],
            reference['analysis_source_sha256'], reference['integrator_sha256']) != (
            384, .4, .41, 6., .01, 0., .0005, fingerprint, runner_hash, integrator_hash):
        raise ValueError('Saved 384^3 comparison has different settings or source')
    if tuple((r['method'], r['dt']) for r in reference['runs']) != configurations:
        raise ValueError('Saved 384^3 method controls are incomplete')
    origin = next(r for r in source['row_sources']
                  if (r['N'], r['dt'], r['time']) == (256, .001, .4))
    path = checkpoint_cache / origin['file']
    if (file_hash(path) != origin['summary_sha256']
            or file_hash(path.with_suffix('.npz')) != origin['arrays_sha256']):
        raise ValueError('Original 256^3 checkpoint hash changed')
    saved = json.loads(path.read_text())
    if (saved['N'], saved['box'], saved['nu'], saved['sigma'],
            saved['P_U'], saved['source_sha256']) != (256, 6., .01, 0., 0., fingerprint):
        raise ValueError('Original equation settings changed')
    with np.load(path.with_suffix('.npz'), allow_pickle=False) as arrays:
        initial = arrays['final']
    initial_observations = observe(initial, 6., .01, .001)
    previous = next(r for r in source['rows']
                    if (r['N'], r['dt'], r['time']) == (256, .001, .4))
    if not np.isclose(initial_observations['energy'], previous['energy'],
                      rtol=2e-12, atol=1e-13):
        raise ValueError('Starting energy differs from prior record')
    common = dict(checkpoint=origin, input_file=input_path.name,
                  input_sha256=file_hash(input_path), numerical_source_sha256=fingerprint)
    runs = []
    for method, dt in configurations:
        points = run_branch(initial, .4, (.405, .410), dt, method, cache, common)
        runs.append(dict(method=method, dt=dt, checkpoints=points,
                         **{k: v for k, v in points[-1].items()
                            if k not in ('provenance', 'file')}))
    del initial
    comparisons = []
    for method in ('euler', 'heun'):
        full, half = [r for r in runs if r['method'] == method]
        for j, stop in enumerate((.405, .410)):
            with np.load((cache / full['checkpoints'][j]['file']).with_suffix('.npz'),
                         allow_pickle=False) as arrays:
                coarse = arrays['final']
            with np.load((cache / half['checkpoints'][j]['file']).with_suffix('.npz'),
                         allow_pickle=False) as arrays:
                fine = arrays['final']
            measured = dict(method=method, time=stop,
                            **difference(coarse, fine, 6.),
                            **compare_fields(coarse, fine, 6.))
            comparisons.append(measured)
            del coarse, fine
            print(f'Compared 256^3 {method} at {stop:g}: '
                  f'energy difference={measured["energy_difference"]:.6e}', flush=True)
    result = dict(
        start_time=.4, end_time=.41, grid=256, input_history_dt=.001,
        external_force=0., box=6., nu=.01, common_source=common,
        numerical_source_sha256=fingerprint,
        analysis_source_sha256=file_hash(__file__),
        runner_sha256=runner_hash, integrator_sha256=integrator_hash,
        initial_observations=initial_observations, runs=runs, comparisons=comparisons,
        reused_reference=dict(
            file=reference_path.name, sha256=file_hash(reference_path),
            grid=384, input_history_dt=.0005, recomputed=False,
            comparisons=reference['comparisons']),
        limits=[
            'Each grid has four controls from identical arrays within that grid.',
            'The 256^3 input has Euler dt 0.001 history; the reused 384^3 input has dt 0.0005 history.',
            'These compare time-step sensitivity within each grid; they do not isolate spatial convergence between grids.',
            'The completed 384^3 controls and original trajectory through 0.40 are reused without evolution.',
            'Full-step/half-step differences are sensitivities, not exact errors or rigorous error bounds.',
            'No long-time stability or Clay proof follows from these finite numerical observations.'])
    out.parent.mkdir(parents=True, exist_ok=True)
    write_json(out, result)
    rows = [dict(method=r['method'], dt=r['dt'], **row) for r in runs for row in r['rows']]
    with out.with_suffix('.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print('256^3 comparison complete; 384^3 results reused without evolution.', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('math/results/hug-resolution.json'))
    parser.add_argument('--reference', type=Path, default=Path('math/results/hug-method-interval.json'))
    parser.add_argument('--checkpoint-cache', type=Path, default=Path('scratch/hug-refinement'))
    parser.add_argument('--cache', type=Path, default=Path('scratch/hug-method-grid/n256'))
    parser.add_argument('--out', type=Path, default=Path('math/results/hug-method-grid.json'))
    args = parser.parse_args()
    study(args.input, args.reference, args.checkpoint_cache, args.cache, args.out)
