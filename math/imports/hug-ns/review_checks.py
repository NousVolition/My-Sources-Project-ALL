"""Small, reproducible checks of the supplied files; no changes to their solver.

Run from any directory: python path/to/review_checks.py --out review-results.json
The reported limitations are measurements, not passing scientific claims.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import platform
import sys

import numpy as np

HERE = Path(__file__).resolve().parent


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


solver = load_module('hug_ns_supplied_solver', 'solver.py')
# run.py uses a local absolute import. Restore sys.modules to avoid changing
# imports for the rest of the repository's tests.
previous = sys.modules.get('solver')
sys.modules['solver'] = solver
try:
    runner = load_module('hug_ns_supplied_runner', 'run.py')
finally:
    if previous is None:
        del sys.modules['solver']
    else:
        sys.modules['solver'] = previous


def operators(n=12):
    x, y, z, *ops = solver.grids(n)
    return x, y, z, tuple(ops)


def alias_measurement(n=12):
    """At N=3m, cos(m*x)^2 aliases its 2m mode onto retained -m."""
    _, _, _, ops = operators(n)
    keep = ops[4]
    m = n // 3
    phase = 2 * np.pi * m * np.arange(n) / n
    wave = np.broadcast_to(np.cos(phase)[:, None, None], (n, n, n))
    product = wave * wave
    filtered = np.fft.ifftn(solver.spec(product, keep)).real
    return {
        'n': n, 'mode': m, 'boundary_mode_retained': bool(keep[m, 0, 0]),
        'max_alias_amplitude_after_filter': float(np.max(np.abs(filtered - 0.5))),
        'expected_for_alias_free_projection': 0.0,
    }


def gate_join_measurement():
    """A tie in max(abs(x),abs(y),abs(z)) inside the transition band."""
    a = solver.L / 2 - 1.5 * solver.EPS
    measurements = []
    for h in (1e-4, 1e-5, 1e-6):
        at = lambda x: float(solver.gate(np.array(x), np.array(a), np.array(0.0)))
        left = (at(a) - at(a-h)) / h
        right = (at(a+h) - at(a)) / h
        measurements.append({'h': h, 'left_derivative': left, 'right_derivative': right,
                             'derivative_jump': abs(left-right)})
    return {'x': a, 'y': a, 'z': 0.0, 'measurements': measurements}


def projection_measurement():
    x, y, z, ops = operators()
    kx, ky, kz, k2, keep, _ = ops
    u = [np.sin(2*np.pi*x/solver.L), np.cos(2*np.pi*y/solver.L),
         np.sin(2*np.pi*(x+z)/solver.L)]
    p = solver.project(u, kx, ky, kz, k2, keep)
    pp = solver.project(p, kx, ky, kz, k2, keep)
    return {'divergence': solver.div_max(p, kx, ky, kz, keep),
            'idempotence_error': max(float(np.max(np.abs(a-b))) for a,b in zip(p,pp)),
            'energy_before': solver.energy(u), 'energy_after': solver.energy(p)}


def shear_measurement():
    x, y, z, ops = operators()
    nu, dt = 0.01, 0.001
    u = [np.sin(2*np.pi*y/solver.L), np.zeros_like(x), np.zeros_like(x)]
    got = solver.advance(u, dt, ops, nu)
    lam = -nu * (2*np.pi/solver.L)**2
    # This flow has zero advection; exact PDE decay and RK2 polynomial are known.
    rk_factor = 1 + lam*dt + 0.5*(lam*dt)**2
    return {'rk2_error': max(float(np.max(np.abs(a-rk_factor*b))) for a,b in zip(got,u)),
            'exact_decay_error': max(float(np.max(np.abs(a-np.exp(lam*dt)*b))) for a,b in zip(got,u)),
            'energy_before': solver.energy(u), 'energy_after': solver.energy(got)}


def energy_rate_measurement():
    """Independent discrete L2 identity on low Fourier modes, away from cutoff."""
    x, y, z, ops = operators(16)
    kx, ky, kz, k2, keep, _ = ops
    phase = 2*np.pi/solver.L
    u = [np.sin(phase*y) + .3*np.cos(phase*z),
         np.sin(phase*z) + .4*np.cos(phase*x),
         np.sin(phase*x) + .2*np.cos(phase*y)]
    nu = .01
    rate = solver.rhs(u, ops, nu)
    observed = solver.L**3 * sum(float(np.mean(a*b)) for a,b in zip(u,rate))
    expected = -nu*solver.L**3 * sum(float(np.mean(solver.deriv(c,k,keep)**2))
                                     for c in u for k in (kx,ky,kz))
    return {'observed_energy_rate': observed, 'viscous_energy_rate': expected,
            'absolute_error': abs(observed-expected)}


def short_runs():
    cases = []
    for name, (psi, nu) in runner.STARTS.items():
        u, ops, initial = solver.build(16, psi=psi, nu=nu)
        kx, ky, kz, _, keep, _ = ops
        # Fixed, deliberately small review step; not the runner's step scheduler.
        dt = min(0.001, solver.timestep(u, ops, nu)/2)
        worst_div = initial['divergence']
        all_finite = all(bool(np.isfinite(c).all()) for c in u)
        largest_energy_increase = 0.0
        last_energy = initial['energy']
        for _ in range(10):
            u = solver.advance(u, dt, ops, nu)
            all_finite &= all(bool(np.isfinite(c).all()) for c in u)
            worst_div = max(worst_div, solver.div_max(u,kx,ky,kz,keep))
            new_energy = solver.energy(u)
            largest_energy_increase = max(largest_energy_increase,new_energy-last_energy)
            last_energy = new_energy
        cases.append({'case': name, 'n':16, 'dt':dt, 'steps':10, 't_end':10*dt,
                      'initial':initial, 'final_peak':solver.peak_speed(u),
                      'final_energy':last_energy, 'finite_at_every_step':all_finite,
                      'max_divergence':worst_div,
                      'largest_step_energy_increase':largest_energy_increase})
    return cases


def runner_measurement():
    u, ops, _ = solver.build(12)
    recommended = solver.timestep(u,ops,.01)
    end = 1.49*recommended
    with contextlib.redirect_stdout(io.StringIO()) as stdout:
        result = runner.evolve('hug',12,end,samples=2)
    # round(1.49) yields one step, so the actual step is 49% above its estimate.
    return {'recommended_initial_dt':recommended, 'requested_end':end,
            'actual_steps':len(result['series'])-1,
            'actual_dt_over_recommended':end/recommended,
            'comparison_terms_in_json': all('steepening' in s and 'smoothing' in s for s in result['series']),
            'comparison_terms_in_stdout': 'steepening=' in stdout.getvalue()}


def pressure_measurement():
    x, y, z, ops = operators(16)
    kx, ky, kz, k2, keep, _ = ops
    phase = 2*np.pi/solver.L
    low = [np.sin(phase*y)+.3*np.cos(phase*z)+.4,
           np.sin(phase*z)+.4*np.cos(phase*x)-.2,
           np.sin(phase*x)+.2*np.cos(phase*y)+.1]
    rng = np.random.default_rng(731)
    high = solver.project([rng.normal(size=x.shape) for _ in range(3)],kx,ky,kz,k2,keep)
    shear = [np.sin(phase*y), np.zeros_like(x), np.zeros_like(x)]
    records = []
    for name, u in [('low_modes',low),('retained_high_modes',high),('analytic_shear',shear)]:
        p = solver.project(u,kx,ky,kz,k2,keep)
        report = solver.pressure_at_peak(u,ops,.01)
        idx = tuple(report['index'])
        rate = solver.rhs(u,ops,.01)
        actual = sum(float(p[i][idx]*rate[i][idx]) for i in range(3))
        records.append({'case':name, 'split':report, 'rhs_local_energy_rate':actual,
                        'absolute_sum_error':abs(report['total_local_energy_rate']-actual)})
    return records


def collect():
    return {
        'scope':'Review of the supplied implementation, not a continuation or a verification of the README long-run claims.',
        'python':platform.python_version(), 'numpy':np.__version__,
        'source_sha256': {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in ('solver.py','run.py','requirements.txt')},
        'projection':projection_measurement(), 'analytic_shear':shear_measurement(),
        'energy_rate':energy_rate_measurement(), 'short_runs':short_runs(),
        'gate_join':gate_join_measurement(),
        'cutoff_aliasing':[alias_measurement(n) for n in (12,48)],
        'runner':runner_measurement(),
        'pressure_split':pressure_measurement(),
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args=parser.parse_args()
    result=collect()
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
