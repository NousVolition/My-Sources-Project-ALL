"""Check the supplied ZIP and expose unresolved assumptions as strict xfails."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

FOLDER = Path(__file__).resolve().parents[1] / 'imports' / 'hug-ns'
spec = importlib.util.spec_from_file_location('hug_ns_review', FOLDER / 'review_checks.py')
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


@pytest.fixture(scope='module')
def measured():
    return review.collect()


@pytest.mark.parametrize('filename', ['README.md', 'BOUND.md', 'requirements.txt'])
def test_reviewed_content_matches_recorded_version(filename):
    manifest = json.loads((FOLDER / 'provenance.json').read_text(encoding='utf-8'))
    source = next(row for row in manifest['source_files'] if row['published_path'] == filename)
    content = (FOLDER / filename).read_bytes()
    # Notes are editable; their original and reviewed hashes are recorded separately.
    assert hashlib.sha256(content).hexdigest() == source['published_sha256']
    if filename == 'requirements.txt':
        assert source['published_sha256'] == source['sha256']


def test_projection_removes_divergence_without_adding_energy(measured):
    r = measured['projection']
    assert r['divergence'] < 1e-12
    assert r['idempotence_error'] < 1e-12
    assert r['energy_after'] <= r['energy_before']


def test_known_shear_decay(measured):
    r = measured['analytic_shear']
    assert r['rk2_error'] < 1e-12
    assert r['exact_decay_error'] < 1e-12
    assert r['energy_after'] < r['energy_before']


def test_low_mode_energy_identity(measured):
    assert measured['energy_rate']['absolute_error'] < 1e-10


@pytest.mark.parametrize('case', ['hug', 'sharper', 'tube', 'tubes'])
def test_short_run_stays_finite_and_divergence_small(measured, case):
    r = next(row for row in measured['short_runs'] if row['case'] == case)
    assert r['finite_at_every_step']
    assert r['max_divergence'] < 1e-11
    assert r['largest_step_energy_increase'] < 1e-11


@pytest.mark.xfail(strict=True, reason='Supplied <= N//3 filter retains the aliased cutoff mode; see REVIEW.md')
def test_cutoff_product_has_no_alias(measured):
    assert max(r['max_alias_amplitude_after_filter'] for r in measured['cutoff_aliasing']) < 1e-12


@pytest.mark.xfail(strict=True, reason='Supplied max-distance gate has a derivative jump at coordinate ties; see REVIEW.md')
def test_gate_join_is_differentiable(measured):
    assert measured['gate_join']['measurements'][-1]['derivative_jump'] < 1e-3


@pytest.mark.xfail(strict=True, reason='Supplied round() scheduler can exceed its own initial step estimate; see REVIEW.md')
def test_runner_respects_initial_step_estimate(measured):
    assert measured['runner']['actual_dt_over_recommended'] <= 1.0 + 1e-12


def test_runner_prints_comparison_terms(measured):
    assert measured['runner']['comparison_terms_in_stdout']


@pytest.mark.parametrize('case', ['low_modes', 'retained_high_modes', 'analytic_shear'])
def test_pressure_split_agrees_with_solver_at_same_point(measured, case):
    row = next(r for r in measured['pressure_split'] if r['case'] == case)
    assert row['absolute_sum_error'] < 1e-10
    if case == 'analytic_shear':
        assert abs(row['split']['pressure']) < 1e-12
        assert abs(row['split']['carry']) < 1e-12
        assert row['split']['smoothing'] < 0


def test_gate_initial_fields_and_evolution_preserved():
    manifest = json.loads((FOLDER / 'provenance.json').read_text(encoding='utf-8'))
    core = (FOLDER/'solver.py').read_bytes().split(b'def pressure_at_peak(',1)[0]
    assert hashlib.sha256(core).hexdigest() == manifest['received_solver_core_sha256']
