"""Resume, provenance and known-solution checks for the longer control."""
import json
import numpy as np
import pytest

from hug_method_interval import run_branch
from hug_time_method import evolve


def field():
    u=np.zeros((3,6,6,6))
    u[0]=np.sin(np.arange(6)*2*np.pi/6)[None,:,None]
    return u


@pytest.mark.parametrize('method',['euler','heun'])
def test_resume_matches_uninterrupted_solution_and_preserves_first_checkpoint(tmp_path,method):
    u=field();before=u.copy()
    first=run_branch(u,0.,(.002,),.001,method,tmp_path,{'input':'known shear'})
    retained=(tmp_path/first[0]['file']).read_bytes()
    points=run_branch(u,0.,(.002,.004),.001,method,tmp_path,{'input':'known shear'})
    assert (tmp_path/first[0]['file']).read_bytes()==retained
    expected,ledger=evolve(u,.004,4,method)
    with np.load((tmp_path/points[-1]['file']).with_suffix('.npz')) as arrays:
        np.testing.assert_array_equal(arrays['final'],expected)
    np.testing.assert_array_equal(u,before)
    assert points[-1]['energy_change']==pytest.approx(ledger['energy_change'],abs=1e-13)
    assert points[-1]['work_quadrature']==pytest.approx(ledger['work_quadrature'],abs=1e-13)
    assert len(points[-1]['rows'])==5
    # Exact modal amplification independently checks the continued solver.
    decay=.01*(2*np.sin(np.pi/6))**2
    amplification=1-.001*decay+(.5*(.001*decay)**2 if method=='heun' else 0.)
    assert points[-1]['observations']['energy']==pytest.approx(54*amplification**8,abs=1e-12)


def test_completed_controls_reuse_without_evolution(tmp_path,monkeypatch):
    import hug_method_interval as module
    first=run_branch(field(),0.,(.002,.004),.001,'heun',tmp_path,{'input':'same'})
    def fail(*args,**kwargs):raise AssertionError('Completed calculation was repeated')
    monkeypatch.setattr(module,'evolve',fail)
    assert run_branch(field(),0.,(.002,.004),.001,'heun',tmp_path,{'input':'same'})==first


def test_changed_input_provenance_and_corrupt_array_are_rejected(tmp_path):
    points=run_branch(field(),0.,(.002,),.001,'euler',tmp_path,{'input':'original'})
    with pytest.raises(ValueError,match='different inputs'):
        run_branch(field(),0.,(.002,),.001,'euler',tmp_path,{'input':'changed'})
    path=tmp_path/points[0]['file']
    saved=json.loads(path.read_text());saved['arrays_sha256']='0'*64
    path.write_text(json.dumps(saved))
    with pytest.raises(ValueError,match='hash changed'):
        run_branch(field(),0.,(.002,),.001,'euler',tmp_path,{'input':'original'})


@pytest.mark.parametrize('stops,dt',[((.003,.002),.001),((0.,),.001),((.0015,),.001),((.002,),0.)])
def test_invalid_schedule_fails_before_evolution(tmp_path,stops,dt):
    with pytest.raises(ValueError):
        run_branch(field(),0.,stops,dt,'heun',tmp_path,{})
