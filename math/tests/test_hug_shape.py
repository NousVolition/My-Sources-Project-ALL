"""Independent small distributions check the spatial-moment diagnostic."""
import json
from pathlib import Path

import numpy as np
import pytest

import hug_shape
from hug_shape import energy_shape, extend_study


def test_two_point_distribution_has_known_center_and_widths():
    u = np.zeros((3, 4, 4, 4))
    u[0, 0, 1, 2] = 1.
    u[2, 3, 2, 0] = 2.
    result = energy_shape(u, box=4.)
    # Positions (-1.5,-.5,.5) and (1.5,.5,-1.5), weights 1 and 4.
    assert result['energy'] == 2.5
    np.testing.assert_allclose([result['center_x'], result['center_y'], result['center_z']], [.9, .3, -1.1])
    assert result['radial_width'] == pytest.approx(np.sqrt(1.6))
    assert result['axial_width'] == pytest.approx(.8)
    assert result['rms_radius'] == pytest.approx(np.sqrt(2.24))


def test_speed_rescaling_preserves_width_and_counts_edge_union_once():
    u = np.zeros((3, 20, 20, 20))
    u[0, 0, 0, 0] = 1.
    u[1, 10, 10, 10] = 2.
    before, after = energy_shape(u), energy_shape(-3*u)
    assert before['edge_energy_fraction'] == pytest.approx(.2)
    assert after['energy'] == pytest.approx(9*before['energy'])
    for key in before.keys() - {'energy'}:
        assert after[key] == pytest.approx(before[key], abs=1e-14)


def test_undefined_or_invalid_inputs_are_rejected():
    for u, box in [(np.zeros((3,4,4,4)),6), (np.ones((3,4,4,4)),0),
                   (np.ones((3,4,4,3)),6), (np.full((3,4,4,4),np.nan),6)]:
        with pytest.raises(ValueError):
            energy_shape(u, box)


def test_extension_reads_only_new_arrays_and_preserves_prior_observations(tmp_path, monkeypatch):
    configurations=((256,.001),(384,.001),(384,.0005))
    previous=dict(end_time=.24,external_force=0.,
        numerical_source_sha256=hug_shape.numerical_source_hash(),
        analysis_source_sha256='earlier-analysis',sources=[dict(file='earlier-checkpoint')],
        rows=[dict(N=n,dt=dt,time=t,energy=2.,radial_width=w,axial_width=w,rms_radius=w*2**.5)
              for n,dt in configurations for t,w in ((0.,2.),(.24,1.8))])
    prior=tmp_path/'prior.json'
    prior.write_text(json.dumps(previous))
    old_bytes=prior.read_bytes()
    opened=[]

    class SavedArrays:
        def __init__(self,path):
            self.path=Path(path)
            self.n=int(self.path.name.split('-')[0][1:])
        def __enter__(self):
            opened.append(self.path.name)
            return self
        def __exit__(self,*args):
            return False
        def __getitem__(self,key):
            assert key=='final'
            return type('SavedVelocity',(),dict(shape=(3,self.n,self.n,self.n)))()

    for n,dt in configurations:
        for t in (.28,.32):
            path=tmp_path/f'n{n}-dt{dt:g}-t{t:g}.json'
            path.write_text(json.dumps(dict(N=n,dt=dt,end_time=t,box=6.,nu=.01,
                sigma=0.,P_U=0.,source_sha256=previous['numerical_source_sha256'],
                rows=[dict(energy=2.)])))
            path.with_suffix('.npz').write_bytes(b'array fixture')
    monkeypatch.setattr(hug_shape.np,'load',SavedArrays)
    monkeypatch.setattr(hug_shape,'energy_shape',lambda *args: dict(
        energy=2.,radial_width=1.7,axial_width=1.7,rms_radius=1.7*2**.5))
    result=extend_study(prior,tmp_path,(.28,.32),tmp_path/'new.json')
    assert len(opened)==6 and all('t0.28' in p or 't0.32' in p for p in opened)
    assert prior.read_bytes()==old_bytes
    assert [r for r in result['rows'] if r['time']<=.24]==previous['rows']
    assert result['sources'][:1]==previous['sources']
    assert result['continuation']['reused_rows']==6
    assert result['continuation']['new_measurements']==6
    assert result['summary'][0]['sample_times']==[0.,.24,.28,.32]
    assert result['summary'][0]['radial_width_change_percent']==pytest.approx(-15.)
    # Equal last two samples must not be labeled strictly decreasing.
    assert result['summary'][0]['radial_width_sample_trend']=='mixed or constant'
    for stops,out in [((.24,),tmp_path/'bad.json'),((.32,.28),tmp_path/'bad.json'),
                      ((.2801,),tmp_path/'bad.json'),((.28,),prior)]:
        with pytest.raises(ValueError):
            extend_study(prior,tmp_path,stops,out)
    assert len(opened)==6
