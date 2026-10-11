"""The recording layer must preserve parent arithmetic and matched initial fields."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_history import run, initial, Flow
from simulate import run as parent_run
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'spatial-convergence'))
from run_refinement import transfer


def test_complete_recording_matches_parent_bit_for_bit(tmp_path):
    c={'n':12,'nu':.02,'dt':.005,'save_dt':.025,'end':.1,'markers':12}
    a=parent_run(617,c,tmp_path/'parent'/'data')
    b=run(617,c,tmp_path/'recorder'/'data')
    with np.load(a,allow_pickle=False) as x,np.load(b,allow_pickle=False) as y:
        assert set(x.files)==set(y.files)
        for key in x.files: np.testing.assert_array_equal(x[key],y[key])


def test_initial_velocity_is_same_polynomial_on_non_nested_grids():
    small=Flow(14,.02,workers=1); large=Flow(22,.02,workers=1)
    a=initial(small,1024); b=initial(large,1024)
    np.testing.assert_allclose(large.real(transfer(a,22)),large.real(b),rtol=0,atol=8e-15)
