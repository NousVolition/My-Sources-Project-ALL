"""A staged continuation must preserve the field and reuse completed work."""
import numpy as np
import pytest

from extend_hug_refinement import resume_through_checkpoints
from hug_refinement import run
from navier import PurePythonNavierStokes3D


def test_staged_continuation_matches_direct_run_and_reuses_cache(tmp_path, monkeypatch):
    baseline, baseline_path = run(8, .001, .004, tmp_path)
    direct, direct_path = run(8, .001, .008, tmp_path/'direct')
    original = PurePythonNavierStokes3D.step_with_pressure
    calls = []

    def counted(self, *args, **kwargs):
        calls.append(args[0])
        return original(self, *args, **kwargs)

    monkeypatch.setattr(PurePythonNavierStokes3D, 'step_with_pressure', counted)
    resumed, resumed_path = resume_through_checkpoints(
        baseline_path.with_suffix('.json'), (.006, .008), tmp_path)
    assert calls == [.001]*4
    assert PurePythonNavierStokes3D.step_with_pressure is counted
    assert resumed['rows'][:len(baseline['rows'])] == baseline['rows']
    assert resumed['rows'][-1] == direct['rows'][-1]
    with np.load(resumed_path) as a, np.load(direct_path) as b:
        assert np.array_equal(a['initial'], b['initial'])
        assert np.array_equal(a['final'], b['final'])
    again, again_path = resume_through_checkpoints(
        baseline_path.with_suffix('.json'), (.006, .008), tmp_path)
    assert again == resumed and again_path == resumed_path
    assert len(calls) == 4


def test_invalid_schedule_is_rejected_before_evolution(tmp_path, monkeypatch):
    _, path = run(8, .001, .004, tmp_path)

    def unexpected_step(*args, **kwargs):
        pytest.fail('Invalid checkpoint schedule must not evolve the field')

    monkeypatch.setattr(PurePythonNavierStokes3D, 'step_with_pressure', unexpected_step)
    for stops in ((), (.004, .008), (.006, .005), (.0065,), (float('nan'),)):
        with pytest.raises(ValueError, match='Checkpoint times'):
            resume_through_checkpoints(path.with_suffix('.json'), stops, tmp_path)
    assert PurePythonNavierStokes3D.step_with_pressure is unexpected_step
