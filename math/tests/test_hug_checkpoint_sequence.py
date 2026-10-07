"""A staged continuation must preserve the field and reuse completed work."""
import json

import numpy as np
import pytest

from extend_hug_refinement import continuation_baseline, resume_through_checkpoints
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


def test_later_published_baseline_supplies_its_endpoint_and_comparisons(tmp_path):
    saved = dict(end_time=.24, external_force=0.,
        runs=[dict(N=n, dt=dt, end_time=.24) for n,dt in
              ((256,.001),(384,.001),(384,.0005))],
        spatial_comparison=dict(coarse=256, fine=384, velocity_relative_l2=.001),
        time_comparison=dict(grid=384, velocity_relative_l2=.0003))
    path = tmp_path/'later.json'
    path.write_text(json.dumps(saved))
    previous, runs, spatial = continuation_baseline(path)
    assert previous == saved
    assert runs == saved['runs']
    assert spatial == saved['spatial_comparison']
    saved['runs'][0]['end_time'] = .16
    path.write_text(json.dumps(saved))
    with pytest.raises(ValueError, match='share the unforced endpoint'):
        continuation_baseline(path)


def test_original_two_study_baseline_is_preserved(tmp_path, monkeypatch):
    folder = tmp_path/'math/results'
    folder.mkdir(parents=True)
    low = dict(N=256, dt=.001, end_time=.16)
    full, half = [dict(N=384, dt=dt, end_time=.16) for dt in (.001,.0005)]
    spatial = dict(velocity_relative_l2=.0006)
    time_study = dict(end_time=.16, runs=[full,half], spatial_comparison=dict(final=spatial))
    (folder/'hug-time-384.json').write_text(json.dumps(time_study))
    (folder/'hug-grid-384.json').write_text(json.dumps(dict(end_time=.16,runs=[low,full])))
    monkeypatch.chdir(tmp_path)
    previous, runs, comparison = continuation_baseline()
    assert previous == time_study
    assert runs == [low,full,half]
    assert comparison == spatial
