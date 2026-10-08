"""Add one 256^3 Heun dt=0.00025 control; reuse both completed larger steps."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from continue_hug_refinement import file_hash, numerical_source_hash
from hug_method_interval import observe, read_point, run_branch, write_json
from hug_refinement import compare_fields
from hug_time_method import difference


def study(reference_path, checkpoint_cache, cache, out):
    reference_path, checkpoint_cache, cache, out = map(
        Path, (reference_path, checkpoint_cache, cache, out))
    reference = json.loads(reference_path.read_text())
    module = Path(__file__).parent
    fingerprint = numerical_source_hash()
    runner_hash = file_hash(module / 'hug_method_interval.py')
    integrator_hash = file_hash(module / 'hug_time_method.py')
    if (reference['grid'], reference['start_time'], reference['end_time'],
            reference['input_history_dt'], reference['box'], reference['nu'],
            reference['external_force'], reference['numerical_source_sha256'],
            reference['runner_sha256'], reference['integrator_sha256'],
            reference['analysis_source_sha256']) != (
            256, .4, .41, .001, 6., .01, 0., fingerprint, runner_hash,
            integrator_hash, file_hash(module / 'hug_method_grid.py')):
        raise ValueError('Require the verified 256^3 method comparison through 0.41')
    full, half = [next(r for r in reference['runs']
                       if r['method'] == 'heun' and r['dt'] == dt)
                  for dt in (.001, .0005)]
    common = reference['common_source']
    if file_hash(reference_path.parent / common['input_file']) != common['input_sha256']:
        raise ValueError('The saved-flow source record changed')
    origin = common['checkpoint']
    path = checkpoint_cache / origin['file']
    if (origin['N'], origin['dt'], origin['time']) != (256, .001, .4):
        raise ValueError('Wrong starting checkpoint')
    if (file_hash(path), file_hash(path.with_suffix('.npz'))) != (
            origin['summary_sha256'], origin['arrays_sha256']):
        raise ValueError('Starting checkpoint hashes changed')
    retained = []
    for run in (full, half):
        provenance = dict(source=common, N=256, start_time=.4, dt=run['dt'],
                          method='heun', box=6., nu=.01, external_force=0.,
                          sample_interval=.001, runner_sha256=runner_hash,
                          integrator_sha256=integrator_hash)
        for point, stop in zip(run['checkpoints'], (.405, .41)):
            saved, _ = read_point(cache / point['file'], provenance)
            if (saved != {k: v for k, v in point.items() if k != 'file'}
                    or saved['end_time'] != stop):
                raise ValueError('Retained Heun checkpoint differs from published data')
            retained.append(dict(file=point['file'], summary_sha256=file_hash(cache / point['file']),
                                 arrays_sha256=point['arrays_sha256']))
    with np.load(path.with_suffix('.npz'), allow_pickle=False) as arrays:
        initial = arrays['final']
    initial_observations = observe(initial, 6., .01, .00025)
    if not np.isclose(initial_observations['energy'],
                      reference['initial_observations']['energy'], rtol=2e-12, atol=1e-13):
        raise ValueError('Initial energy differs from saved reference')
    # Only this third time step is evolved. Larger-step branches stay untouched.
    points = run_branch(initial, .4, (.405, .410), .00025, 'heun', cache, common)
    del initial
    new_run = dict(method='heun', dt=.00025, checkpoints=points,
                   **{k: v for k, v in points[-1].items() if k not in ('file', 'provenance')})
    comparisons = []
    reductions = []
    for j, stop in enumerate((.405, .41)):
        previous = next(c for c in reference['comparisons']
                        if c['method'] == 'heun' and c['time'] == stop)
        comparisons.append(dict(pair='full_half', coarse_dt=.001, fine_dt=.0005,
                                reused=True, **previous))
        with np.load((cache / half['checkpoints'][j]['file']).with_suffix('.npz'),
                     allow_pickle=False) as arrays:
            middle = arrays['final']
        with np.load((cache / points[j]['file']).with_suffix('.npz'),
                     allow_pickle=False) as arrays:
            fine = arrays['final']
        measured = dict(pair='half_quarter', coarse_dt=.0005, fine_dt=.00025,
                        reused=False, method='heun', time=stop,
                        **difference(middle, fine, 6.), **compare_fields(middle, fine, 6.))
        del middle, fine
        comparisons.append(measured)
        reductions.append(dict(
            time=stop,
            energy_factor=abs(previous['energy_difference'] / measured['energy_difference']),
            velocity_l2_factor=previous['velocity_l2_difference'] / measured['velocity_l2_difference'],
            gradient_relative_factor=previous['gradient_relative_l2'] / measured['gradient_relative_l2']))
        print(f'Compared Heun steps at {stop:g}: {reductions[-1]}', flush=True)
    energy_series = []
    for a, b, c in zip(full['rows'], half['rows'], new_run['rows']):
        if a['time'] != b['time'] or b['time'] != c['time']:
            raise ValueError('Comparison times differ')
        energy_series.append(dict(time=a['time'], full_energy=a['energy'], half_energy=b['energy'],
                                  quarter_energy=c['energy'],
                                  full_half_gap=abs(a['energy'] - b['energy']),
                                  half_quarter_gap=abs(b['energy'] - c['energy'])))
    result = dict(
        grid=256, start_time=.4, end_time=.41, input_history_dt=.001,
        box=6., nu=.01, external_force=0., common_source=common,
        numerical_source_sha256=fingerprint, analysis_source_sha256=file_hash(__file__),
        runner_sha256=runner_hash, integrator_sha256=integrator_hash,
        initial_observations=initial_observations, new_run=new_run,
        reused_reference=dict(file=reference_path.name, sha256=file_hash(reference_path),
                              methods=['heun'], steps=[.001, .0005], evolution_repeated=False,
                              checkpoints=retained),
        comparisons=comparisons, reductions=reductions, energy_series=energy_series,
        limits=[
            'All three Heun branches start from exactly the same 256^3 state at 0.40.',
            'Only dt 0.00025 is new; dt 0.001 and 0.0005 are reused without evolution.',
            'This tests time-step sensitivity over 0.40 to 0.41 at fixed spatial resolution.',
            'The starting field retains its previous Euler history; this is not Heun from zero.',
            'Shrinking numerical differences do not supply exact continuum errors, a rigorous convergence bound, or a Clay proof.'])
    out.parent.mkdir(parents=True, exist_ok=True)
    write_json(out, result)
    with out.with_suffix('.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(energy_series[0]))
        writer.writeheader()
        writer.writerows(energy_series)
    print('Quarter-step Heun comparison complete; both prior branches reused.', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, default=Path('math/results/hug-method-grid.json'))
    parser.add_argument('--checkpoint-cache', type=Path, default=Path('scratch/hug-refinement'))
    parser.add_argument('--cache', type=Path, default=Path('scratch/hug-method-grid/n256'))
    parser.add_argument('--out', type=Path, default=Path('math/results/hug-method-refinement.json'))
    args = parser.parse_args()
    study(args.reference, args.checkpoint_cache, args.cache, args.out)
