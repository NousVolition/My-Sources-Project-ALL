"""Centered-drag staggered leapfrog; a dimensionless spring validation.

Equation: x'' = force(x)/mass - gamma*x'. Velocities live at half steps.
The starting half-step is constructed from the supplied on-step velocity.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm


def leapfrog(force, x0, v0, gamma, dt, steps, mass=1.0):
    if gamma < 0 or dt <= 0 or mass <= 0:
        raise ValueError('gamma >= 0, dt > 0 and mass > 0 required')
    # Center drag across the two half-step velocities.
    q = gamma*dt/2
    a, b = (1-q)/(1+q), dt/(1+q)
    x = float(x0)
    half_v = (1+q)*v0 - dt*force(x)/(2*mass)
    xs, vs = [x], [float(v0)]
    for _ in range(steps):
        half_v = a*half_v + b*force(x)/mass
        x += dt*half_v
        next_half = a*half_v + b*force(x)/mass
        xs.append(x)
        vs.append((half_v+next_half)/2)
    return np.arange(steps+1)*dt, np.array(xs), np.array(vs)


def validate(out):
    out.mkdir(parents=True, exist_ok=True)
    rows, arrays = [], {}
    for gamma in [0.0, 0.5, 2.0, 5.0]:
        errors = []
        for dt in [0.02, 0.01]:
            t, x, v = leapfrog(lambda x: -x, 1, 0, gamma, dt, round(12/dt))
            matrix = np.array([[0, 1], [-1, -gamma]])
            exact = np.array([expm(matrix*s) @ [1., 0.] for s in t])
            errors.append(float(np.max(abs(x-exact[:, 0]))))
        assert errors[1] < errors[0]/3.5
        assert np.isfinite(x).all()
        crossings = int(np.count_nonzero(x[1:]*x[:-1] < 0))
        if gamma >= 2:
            assert crossings == 0 and np.all(np.diff(x) <= 1e-12)
        arrays[f'g{gamma:g}_time'] = t
        arrays[f'g{gamma:g}_position'] = x
        arrays[f'g{gamma:g}_velocity'] = v
        rows.append(dict(gamma=gamma, regime='undamped' if gamma == 0 else 'underdamped' if gamma < 2 else 'critical' if gamma == 2 else 'overdamped',
                         dt=0.01, max_position_error=errors[1], error_reduction_on_halving_dt=errors[0]/errors[1], zero_crossings=crossings))
    # Free drag must converge to v=v0*exp(-gamma*t), x=x0+v0*(1-exp(-gamma*t))/gamma.
    t, x, v = leapfrog(lambda x: 0, 0, 1, 5, .001, 1000)
    free_error = float(np.max(abs(x-(-np.expm1(-5*t)/5))))
    assert free_error < 2e-6
    result = {'spring_mass': 1, 'spring_k': 1, 'x0': 1, 'v0': 0, 'time_end': 12,
              'rows': rows, 'free_drag_max_position_error': free_error,
              'scope': 'Dimensionless analytic oscillator validation, not water molecular dynamics.'}
    (out/'oscillator.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    np.savez_compressed(out/'oscillator.npz', **arrays)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent/'data')
    args = parser.parse_args()
    validate(args.out)
