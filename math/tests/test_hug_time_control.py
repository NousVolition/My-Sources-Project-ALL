"""Time-step controls must share the same initial field and reuse saved work."""
import numpy as np
import pytest

from extend_hug_refinement import run_half_step_sequence
from hug_refinement import run
from navier import PurePythonNavierStokes3D


def test_half_step_sequence_matches_direct_run_and_reuses_cache(tmp_path, monkeypatch):
    full, full_path = run(8, .001, .008, tmp_path)
    half, half_path = run_half_step_sequence(full, full_path, tmp_path, (.004, .008))
    direct, direct_path = run(8, .0005, .008, tmp_path/'direct')
    with np.load(half_path) as a, np.load(direct_path) as b:
        assert np.array_equal(a['initial'], b['initial'])
        assert np.array_equal(a['final'], b['final'])
    assert half['rows'][-1] == direct['rows'][-1]

    def unexpected_step(*args, **kwargs):
        pytest.fail('A completed time-step control must be reused')
    monkeypatch.setattr(PurePythonNavierStokes3D, 'step_with_pressure', unexpected_step)
    again, again_path = run_half_step_sequence(full, full_path, tmp_path, (.004, .008))
    assert again == half and again_path == half_path


def test_different_initial_field_is_rejected_before_first_step(tmp_path, monkeypatch):
    full, full_path = run(8, .001, .008, tmp_path)
    with np.load(full_path) as saved:
        initial, final = saved['initial'], saved['final']
    initial[0, 0, 0, 0] += 1.
    np.savez_compressed(full_path, initial=initial, final=final)

    def unexpected_step(*args, **kwargs):
        pytest.fail('An incompatible initial field must be rejected before evolution')
    monkeypatch.setattr(PurePythonNavierStokes3D, 'step_with_pressure', unexpected_step)
    with pytest.raises(ValueError, match='initial arrays differ'):
        run_half_step_sequence(full, full_path, tmp_path, (.004, .008))
    assert PurePythonNavierStokes3D.step_with_pressure is unexpected_step
