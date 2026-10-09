"""Independent regressions for the separate repair; historical xfails stay."""
import importlib.util
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

FOLDER = Path(__file__).resolve().parents[1] / 'imports/hug-ns/corrected_v1'
spec = importlib.util.spec_from_file_location('hug_ns_corrected_v1', FOLDER / '__init__.py',
                                             submodule_search_locations=[str(FOLDER)])
package = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = package
spec.loader.exec_module(package)
from hug_ns_corrected_v1 import run, solver as s


def setup(n=12):
    x, y, z, *ops = s.grids(n)
    return (x, y, z), tuple(ops)


def test_gate_has_no_coordinate_tie_corner():
    a = s.L / 2 - 1.5 * s.EPS
    for signs in [(1, 1, 0), (-1, 1, 0), (1, 1, 1)]:
        x, y, z = [np.array([a * v]) for v in signs]
        jumps = []
        for h in [1e-4, 1e-5, 1e-6]:
            left = (s.gate(x, y, z) - s.gate(x-h, y, z)) / h
            right = (s.gate(x+h, y, z) - s.gate(x, y, z)) / h
            jumps.append(float(abs(left-right)[0]))
        assert jumps[-1] < 1e-4


def test_gate_plateau_collar_and_reflection():
    a = np.linspace(-3, 3, 601)
    zero = np.zeros_like(a)
    w = s.gate(a, a / 2, zero)
    assert np.all(w[np.abs(a) <= 2.1] == 1)
    assert np.all(w[np.abs(a) >= 2.55] == 0)
    assert np.all((w >= 0) & (w <= 1))
    np.testing.assert_allclose(w, s.gate(-a, -a / 2, zero), atol=0, rtol=0)
    # Flat endpoint probes, including the abs-coordinate origin.
    for center, expected in [(0, 1), (2.1, 1), (2.55, 0), (3, 0)]:
        probes = np.array([center-1e-4, center, center+1e-4])
        np.testing.assert_allclose(s.gate(probes, zero[:3], zero[:3]), expected, atol=1e-14)


@pytest.mark.parametrize('eps', [0, -1, 1.5, np.inf, np.nan])
def test_gate_rejects_invalid_collar(eps):
    with pytest.raises(ValueError):
        s.gate(np.zeros(1), np.zeros(1), np.zeros(1), eps=eps)


@pytest.mark.parametrize('n', [12, 16, 17, 18, 48])
def test_retained_quadratic_products_match_exact_convolution(n):
    _, ops = setup(n)
    keep = ops[-2]
    k = (n-1) // 3
    # Real trigonometric polynomials containing near-cutoff modes on all axes.
    a = {(k, k, 0): .3+.2j, (0, 1, -k): -.1j, (1, 0, 0): .4}
    b = {(k, 0, k): .2, (0, k, -1): .25j, (1, 0, 0): .5}
    for coeff in [a, b]:
        coeff.update({tuple(-v for v in key): complex(value).conjugate()
                      for key, value in list(coeff.items())})
    def field(coeff):
        arr = np.zeros((n, n, n), dtype=complex)
        for key, value in coeff.items():
            arr[tuple(v % n for v in key)] = value * n**3
        return np.fft.ifftn(arr).real
    product = s.spec(field(a) * field(b), keep) / n**3
    exact = np.zeros_like(product)
    # No modular wrapping in the reference convolution; discard nonretained sums.
    for ka, va in a.items():
        for kb, vb in b.items():
            kc = tuple(x+y for x, y in zip(ka, kb))
            if all(3 * abs(v) < n for v in kc):
                exact[tuple(v % n for v in kc)] += va * vb
    np.testing.assert_allclose(product, exact, atol=2e-15, rtol=0)
    if n % 3 == 0:
        assert not keep[n // 3, 0, 0]
        assert not keep[0, n // 3, 0]
        assert not keep[0, 0, n // 3]


def test_old_rounding_reproducer_respects_actual_step_limits():
    (x, y, z), ops = setup()
    u = [np.sin(2*np.pi*y/s.L), np.zeros_like(x), np.zeros_like(x)]
    estimate = s.timestep(u, ops, .01)
    final, rows, steps = run.integrate(u, ops, .01, 1.49 * estimate, samples=1)
    assert len(steps) >= 2
    assert max(row['dt'] for row in steps) <= estimate
    assert rows[-1]['t'] == 1.49 * estimate
    assert steps[-1]['t_end'] == rows[-1]['t']
    assert s.energy(final) < s.energy(u)


def test_scheduler_recomputes_for_accelerating_state(monkeypatch):
    (x, _, _), ops = setup()
    u = [np.ones_like(x), np.zeros_like(x), np.zeros_like(x)]
    def manufactured_growth(state, dt, ops, nu):
        result = [state[0] + 10*dt, state[1].copy(), state[2].copy()]
        return result, result
    monkeypatch.setattr(s, 'heun_trial', manufactured_growth)
    final, rows, steps = run.integrate(u, ops, 0, .15, samples=1)
    assert sum(row['rejections'] for row in steps) > 0
    assert steps[-1]['start_limit'] < steps[0]['start_limit']
    for row in steps:
        assert row['dt'] <= min(row[key] for key in ['start_limit', 'predictor_limit', 'end_limit'])
        assert row['t_end'] > row['t_start']
    assert rows[-1]['t'] == .15
    np.testing.assert_allclose(final[0], 2.5, atol=2e-15)


def test_nonfinite_trials_are_rejected(monkeypatch):
    (x, _, _), ops = setup()
    u = [np.ones_like(x), np.zeros_like(x), np.zeros_like(x)]
    def trial(state, dt, ops, nu):
        output = [np.full_like(x, np.nan)] * 3 if dt > .01 else [c.copy() for c in state]
        return output, output
    monkeypatch.setattr(s, 'heun_trial', trial)
    final, _, steps = run.integrate(u, ops, .01, .02, samples=1)
    assert steps[0]['rejections'] == 1
    assert all(np.isfinite(c).all() for c in final)
    assert all(row['dt'] <= .01 for row in steps)


@pytest.mark.parametrize('end', [0, -1, np.inf, np.nan])
def test_invalid_horizon_is_rejected(end):
    (x, _, _), ops = setup()
    with pytest.raises(ValueError):
        run.integrate([x*0]*3, ops, .01, end)


def test_zero_flow_zero_viscosity_and_invalid_viscosity():
    (x, _, _), ops = setup()
    u = [x*0]*3
    assert s.checker(u) == 0
    final, rows, steps = run.integrate(u, ops, 0, .01, samples=3)
    assert rows[-1]['t'] == .01
    assert s.energy(final) == 0
    for nu in [-.1, np.nan, np.inf]:
        with pytest.raises(ValueError):
            run.integrate(u, ops, nu, .01)


def test_heun_has_second_order_convergence_on_exact_shear():
    (x, y, _), ops = setup()
    u = [np.sin(2*np.pi*y/s.L), x*0, x*0]
    errors = []
    for steps in [4, 8, 16]:
        v = [c.copy() for c in u]
        dt = .2/steps
        for _ in range(steps):
            v = s.advance(v, dt, ops, .5)
        exact = np.exp(-.5*(2*np.pi/s.L)**2 * .2) * u[0]
        errors.append(float(np.max(abs(v[0]-exact))))
    assert 3.8 < errors[0]/errors[1] < 4.3
    assert 3.8 < errors[1]/errors[2] < 4.3


def test_energy_identity_and_pressure_split_on_retained_high_modes():
    (x, _, _), ops = setup(18)
    rng = np.random.default_rng(921)
    u = s.project([rng.normal(size=x.shape) for _ in range(3)], *ops[:-1])
    rate = s.rhs(u, ops, .01)
    actual = s.L**3 * sum(np.mean(a*b) for a, b in zip(u, rate))
    expected = -.01*s.L**3*sum(np.mean(s.deriv(c, k, ops[-2])**2)
                              for c in u for k in ops[:3])
    assert abs(actual-expected) < 1e-11
    split = s.pressure_at_peak(u, ops, .01)
    ip = tuple(split['index'])
    assert abs(split['total_local_energy_rate']-sum(a[ip]*b[ip] for a, b in zip(u, rate))) < 1e-12
    assert s.div_max(u, *ops[:3], ops[-2]) < 1e-12


@pytest.mark.parametrize('case', list(run.STARTS))
def test_four_short_starts_record_limits_and_diagnostics(case):
    result = run.evolve(case, 16, .01, samples=10)
    assert result['version'] == 'hug-ns-corrected-v1'
    assert result['series'][-1]['t'] == .01
    assert len(result['series']) == 11
    assert max(row['divergence'] for row in result['series']) < 1e-11
    assert np.max(np.diff([row['energy'] for row in result['series']])) < 1e-11
    for row in result['step_log']:
        assert row['dt'] <= min(row[key] for key in ['start_limit', 'predictor_limit', 'end_limit'])
    json.dumps(result, allow_nan=False)


def test_cli_preserves_existing_output(tmp_path):
    output = tmp_path / 'existing.json'
    output.write_text('original result', encoding='utf-8')
    command = [sys.executable, '-B', str(FOLDER/'run.py'), 'hug', '--n', '8',
               '--time', '.001', '--samples', '1', '--out', str(output)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert 'output already exists' in result.stderr
    assert output.read_text(encoding='utf-8') == 'original result'
    fresh = tmp_path / 'new.json'
    command[-1] = str(fresh)
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(fresh.read_text(encoding='utf-8'))['series'][-1]['t'] == .001


def test_saved_verification_matches_versioned_and_historical_sources():
    record = json.loads((FOLDER/'verification.json').read_text(encoding='utf-8'))
    assert record['version'] == s.VERSION
    for key, folder in [('source_sha256', FOLDER), ('historical_source_sha256', FOLDER.parent)]:
        for name, expected in record[key].items():
            assert hashlib.sha256((folder/name).read_bytes()).hexdigest() == expected
    assert record['short_run_summary']['max_recorded_divergence'] < 1e-11
    assert record['short_run_summary']['largest_sample_energy_change'] < 0
    assert record['short_run_summary']['max_accepted_dt_over_limit'] <= 1
