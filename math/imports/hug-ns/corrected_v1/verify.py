"""Save bounded repair checks; no historical trajectory or live study is resumed."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform

import numpy as np

if __package__:
    from . import run, solver as s
else:
    import run
    import solver as s


def collect():
    folder = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('historical_hug_solver', folder.parent / 'solver.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    joins = []
    a = s.L/2 - 1.5*s.EPS
    for h in [1e-4, 1e-5, 1e-6]:
        row = {'h': h}
        for label, module in [('historical', old), ('corrected_v1', s)]:
            x, y, z = np.array([a]), np.array([a]), np.array([0.0])
            left = float(((module.gate(x, y, z)-module.gate(x-h, y, z))/h)[0])
            right = float(((module.gate(x+h, y, z)-module.gate(x, y, z))/h)[0])
            row[label] = {'left_derivative': left, 'right_derivative': right,
                          'absolute_jump': abs(left-right)}
        joins.append(row)
    cutoff = []
    for n in [12, 48]:
        # Historical endpoint probe; corrected v1 excludes its input mode.
        row = {'n': n}
        for label, module in [('historical', old), ('corrected_v1', s)]:
            x, y, z, *ops = module.grids(n)
            keep = ops[-2]
            f = np.cos(2*np.pi*(n//3)*x/s.L)
            row[label] = {
                'endpoint_retained': bool(keep[n//3, 0, 0]),
                'endpoint_square_alias_after_filter': float(np.max(abs(np.fft.ifftn(module.spec(f*f, keep)).real-.5))),
            }
        cutoff.append(row)
    comparisons, runs = [], []
    for name, (psi, nu) in run.STARTS.items():
        old_u, _, old_start = old.build(16, psi=psi, nu=nu)
        new_u, _, new_start = s.build(16, psi=psi, nu=nu)
        comparisons.append({
            'case': name, 'n': 16,
            'initial_velocity_relative_l2_change': float(np.sqrt(sum(np.sum((a-b)**2) for a, b in zip(new_u, old_u))/sum(np.sum(b*b) for b in old_u))),
            'historical_start': old_start, 'corrected_v1_start': new_start,
            'historical_mean': [float(c.mean()) for c in old_u],
            'corrected_v1_mean': [float(c.mean()) for c in new_u],
        })
        runs.append(run.evolve(name, 16, .01, samples=10))
    x, y, z, *ops = s.grids(12)
    ops = tuple(ops)
    shear = [np.sin(2*np.pi*y/s.L), x*0, x*0]
    estimate = s.timestep(shear, ops, .01)
    _, rounding_rows, rounding_steps = run.integrate(shear, ops, .01, 1.49*estimate, samples=1)
    errors = []
    for steps in [4, 8, 16]:
        u = [c.copy() for c in shear]
        for _ in range(steps):
            u = s.advance(u, .2/steps, ops, .5)
        exact = shear[0]*np.exp(-.5*(2*np.pi/s.L)**2*.2)
        errors.append({'steps': steps, 'dt': .2/steps, 'max_exact_shear_error': float(np.max(abs(u[0]-exact)))})
    max_ratio = max(row['dt']/min(row[key] for key in ['start_limit', 'predictor_limit', 'end_limit'])
                    for result in runs for row in result['step_log'])
    return {
        'version': s.VERSION, 'python': platform.python_version(), 'numpy': np.__version__,
        'scope': 'Four 16^3 checks to 0.01, exact-shear time refinement, and direct gate/cutoff/scheduler probes; no physical convergence claim.',
        'parent_commit': '8e1525f64b27909afb6dcf704848b8666ad2260a',
        'source_sha256': {name: hashlib.sha256((folder/name).read_bytes()).hexdigest()
                          for name in ['__init__.py', 'solver.py', 'run.py', 'verify.py']},
        'historical_source_sha256': {name: hashlib.sha256((folder.parent/name).read_bytes()).hexdigest()
                                     for name in ['solver.py', 'run.py', 'provenance.json', 'review-results.json']},
        'gate_join': joins, 'cutoff_endpoint_probe': cutoff,
        'initial_field_changes': comparisons,
        'rounding_reproducer': {'requested_time': 1.49*estimate, 'initial_limit': estimate,
                               'actual_final_time': rounding_rows[-1]['t'], 'step_log': rounding_steps},
        'exact_shear_refinement': errors, 'short_runs': runs,
        'short_run_summary': {
            'count': len(runs), 'grid': 16, 'end_time': .01,
            'max_recorded_divergence': max(row['divergence'] for r in runs for row in r['series']),
            'largest_sample_energy_change': max(float(np.max(np.diff([row['energy'] for row in r['series']]))) for r in runs),
            'max_accepted_dt_over_limit': max_ratio,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output already exists; choose a new path')
    result = collect()
    with args.out.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result['short_run_summary']))


if __name__ == '__main__':
    main()
