"""Small, independently versioned runs with audited adaptive Heun steps."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

if __package__:
    from . import solver as s
else:
    import solver as s

STARTS = {
    'hug': (s.ring, 0.01),
    'sharper': (lambda x, y, z: s.ring(x, y, z, sharp=4.0), 0.01),
    'tube': (s.pulled_tube, 0.002),
    'tubes': (s.two_tubes, 0.002),
}


def snapshot(u, ops, nu, t):
    kx, ky, kz, k2, keep, dx = ops
    steep, smooth = s.terms(u, ops, nu)
    return {
        't': float(t), 'peak': s.peak_speed(u), 'energy': s.energy(u),
        'slope': s.slope_max(u, kx, ky, kz, keep),
        'divergence': s.div_max(u, kx, ky, kz, keep),
        'steepening': steep, 'smoothing': smooth,
        'peak_split': s.pressure_at_peak(u, ops, nu),
    }


def integrate(u, ops, nu, t_end, samples=12, max_steps=1_000_000):
    """Recompute the estimate at every step; shorten to each output time.

    A trial is halved if either the predictor or result requires a smaller
    step, or becomes nonfinite. This enforces this code's estimate only;
    it is not an error estimator or a proof of Heun stability/accuracy.
    """
    if not np.isfinite(t_end) or t_end <= 0:
        raise ValueError('t_end must be finite and positive')
    for name, value in [('samples', samples), ('max_steps', max_steps)]:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f'{name} must be a positive integer')
    s.timestep(u, ops, nu)  # Validate before diagnostics or stepping.
    t = 0.0
    series = [snapshot(u, ops, nu, t)]
    step_log = []
    for target in np.linspace(0.0, t_end, samples + 1)[1:]:
        target = float(target)
        while t < target:
            if len(step_log) >= max_steps:
                raise RuntimeError(f'step budget exhausted at t={t!r}')
            start_limit = s.timestep(u, ops, nu)
            remaining = target - t
            dt = min(start_limit, remaining)
            for retries in range(33):
                if t + dt <= t or dt <= 0:
                    raise RuntimeError(f'time step underflow at t={t!r}')
                with np.errstate(over='ignore', invalid='ignore'):
                    candidate, predictor = s.heun_trial(u, dt, ops, nu)
                if all(np.isfinite(c).all() for c in [*predictor, *candidate]):
                    stage_limit = s.timestep(predictor, ops, nu)
                    end_limit = s.timestep(candidate, ops, nu)
                    if dt <= min(start_limit, stage_limit, end_limit):
                        break
                dt *= 0.5
            else:
                raise RuntimeError(f'no admissible finite step at t={t!r}')
            next_t = target if dt == remaining else t + dt
            step_log.append({
                't_start': t, 't_end': next_t, 'dt': dt,
                'start_limit': start_limit, 'predictor_limit': stage_limit,
                'end_limit': end_limit, 'rejections': retries,
            })
            u, t = candidate, next_t
        series.append(snapshot(u, ops, nu, target))
    return u, series, step_log


def evolve(name, n, t_end, samples=12):
    psi, nu = STARTS[name]
    u, ops, report = s.build(n, psi=psi, nu=nu)
    _, series, step_log = integrate(u, ops, nu, t_end, samples=samples)
    folder = Path(__file__).resolve().parent
    return {
        'version': s.VERSION, 'case': name, 'n': n, 'nu': nu,
        'domain_length': s.L, 'force': 0, 'method': 'Heun',
        'gate': 'product of three flat coordinate gates; new starting field',
        'cutoff': '3*abs(integer_mode)<N on each axis',
        'source_sha256': {name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
                          for name in ['solver.py', 'run.py']},
        't_end': t_end, 'samples': samples, 'steps': len(step_log),
        'start': report, 'series': series, 'step_log': step_log,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=STARTS)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--time', type=float, default=0.01)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output already exists; choose a new version-specific path')
    if not args.out.parent.is_dir():
        parser.error('output parent directory must already exist')
    result = evolve(args.case, args.n, args.time, args.samples)
    with args.out.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(f"{s.VERSION}: {result['steps']} accepted steps, t={result['series'][-1]['t']}; {args.out}")


if __name__ == '__main__':
    main()
