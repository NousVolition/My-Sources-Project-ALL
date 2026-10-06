"""A saved restart must give the same field as uninterrupted evolution."""
import json

import numpy as np
import pytest

from continue_hug_refinement import continue_run
from hug_refinement import run


def test_checkpoint_restart_matches_uninterrupted_evolution(tmp_path):
    saved, arrays = run(8, .001, .004, tmp_path)
    resumed, resumed_path = continue_run(arrays.with_suffix('.json'), .008, tmp_path)
    direct, direct_path = run(8, .001, .008, tmp_path / 'direct')
    with np.load(resumed_path) as a, np.load(direct_path) as b:
        assert np.array_equal(a['initial'], b['initial'])
        assert np.array_equal(a['final'], b['final'])
    assert resumed['rows'][:len(saved['rows'])] == saved['rows']
    assert resumed['rows'][-1] == direct['rows'][-1]
    assert resumed['continuation']['resumed_from_time'] == .004
    # A matching completed continuation is reused without any new steps.
    again, again_path = continue_run(arrays.with_suffix('.json'), .008, tmp_path)
    assert again == resumed and again_path == resumed_path


def test_checkpoint_from_changed_numerical_source_is_rejected(tmp_path):
    _, arrays = run(8, .001, .004, tmp_path)
    path = arrays.with_suffix('.json')
    saved = json.loads(path.read_text())
    saved['source_sha256'] = 'different-source'
    path.write_text(json.dumps(saved))
    with pytest.raises(ValueError, match='numerical source'):
        continue_run(path, .008, tmp_path)
