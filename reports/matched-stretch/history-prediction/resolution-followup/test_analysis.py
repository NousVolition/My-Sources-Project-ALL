"""Analytic deformation and material-label alignment checks for comparisons."""
import numpy as np
import pytest
from analyze import material_comparison, source_hash_matches


def flow(rate):
    time = np.array([0., .1, .2, .3])
    tangent = np.array([np.diag([np.exp(rate*t), np.exp(-rate*t), 1.]) for t in time])[:, None]
    return {"time": time, "positions": np.zeros((4, 1, 3)), "tangent": tangent}


def test_known_incompressible_strain_difference():
    p = {"anchors": [.1], "gap": .1, "primary_horizon": .1}
    result = material_comparison(flow(.5), flow(.6), p, .55)
    assert result["target_RMSE"] == pytest.approx(.1)
    assert result["relative_target_RMSE"] == pytest.approx(1/6)
    assert result["target_event_Jaccard"] == 0.
    assert result["target_event_flip_fraction"] == 1.
    assert result["fine_target_volume_error_max"] < 1e-14


def test_different_material_labels_are_rejected():
    a, b = flow(.5), flow(.6)
    b["positions"][0, 0, 0] = .01
    with pytest.raises(AssertionError, match="material labels"):
        material_comparison(a, b, {"anchors": [.1], "gap": .1, "primary_horizon": .1}, .55)


def test_different_times_are_rejected():
    a, b = flow(.5), flow(.6)
    b["time"][1] += .01
    with pytest.raises(AssertionError, match="observation times"):
        material_comparison(a, b, {"anchors": [.1], "gap": .1, "primary_horizon": .1}, .55)


def test_source_check_accepts_only_newline_equivalence():
    import hashlib
    recorded = hashlib.sha256(b"x = 1\r\n").hexdigest()
    assert source_hash_matches(b"x = 1\n", recorded)
    assert not source_hash_matches(b"x = 2\n", recorded)
