"""Independent checks for the alternate time integrator and local comparison."""
import numpy as np
import pytest

from hug_time_method import compare_methods, energy, evolve
from navier import PurePythonNavierStokes3D


def shear(n,box):
    u=np.zeros((3,n,n,n))
    u[0]=np.sin(np.arange(n)*2*np.pi/n)[None,:,None]
    return u


@pytest.mark.parametrize('method',['euler','heun'])
@pytest.mark.parametrize('amplitude',[0.,2.])
def test_constant_is_unchanged(method,amplitude):
    u=np.full((3,6,6,6),amplitude)
    v,m=evolve(u,.01,2,method)
    np.testing.assert_array_equal(v,u)
    assert m['temporal_energy_defect']==m['energy_change']==0.


@pytest.mark.parametrize('n',[8,10])
@pytest.mark.parametrize('method',['euler','heun'])
def test_shear_has_known_amplification_and_energy(n,method):
    box,nu,dt=4.,.07,.02
    u=shear(n,box);before=u.copy()
    decay=nu*(2*np.sin(np.pi/n)/(box/n))**2
    factor=1-dt*decay+(dt**2*decay**2/2 if method=='heun' else 0.)
    v,m=evolve(u,2*dt,2,method,box,nu)
    np.testing.assert_allclose(v,u*factor**2,atol=2e-15)
    np.testing.assert_array_equal(u,before)
    assert m['energy_after']==pytest.approx(energy(u,box)*factor**4,abs=2e-13)
    assert abs(m['total_identity_residual'])<1e-12


def test_heun_matches_two_original_python_maps_and_average():
    n,box,dt,nu=6,3.,.002,.02
    sim=PurePythonNavierStokes3D(n,box/n)
    sim.u,sim.v,sim.w=np.random.default_rng(73).normal(size=(3,n,n,n))
    sim.project(backend='fft')
    u=np.asarray((sim.u,sim.v,sim.w));before=u.copy()
    sim.nu,sim.sigma=nu,0.
    for _ in range(2):
        sim.step_with_pressure(dt,0.,backend='fft',step_backend='python')
    expected=.5*u+.5*np.asarray((sim.u,sim.v,sim.w))
    actual,m=evolve(u,dt,1,'heun',box,nu)
    np.testing.assert_allclose(actual,expected,rtol=3e-14,atol=3e-14)
    np.testing.assert_array_equal(u,before)
    assert m['max_stage_divergence']<1e-12


def test_refinement_against_exact_semidiscrete_shear_solution():
    u=shear(12,4.)
    data=compare_methods(u,.1,4.,.1)
    exact=u*np.exp(-.1*.1*(2*np.sin(np.pi/12)/(4/12))**2)
    errors={}
    for method in ('euler','heun'):
        errors[method]=[]
        for count in (1,2,4):
            got,m=evolve(u,.1,count,method,4.,.1)
            errors[method].append(np.linalg.norm(got-exact))
        ratios=np.array(errors[method][:-1])/errors[method][1:]
        target=2 if method=='euler' else 4
        np.testing.assert_allclose(ratios,target,rtol=.03)
    assert errors['heun'][0]<errors['euler'][0]/50
    assert len(data['rows'])==6 and len(data['comparisons'])==4


@pytest.mark.parametrize('changes',[{'interval':0.},{'substeps':0},{'substeps':1.5},
    {'method':'unknown'},{'box':0.},{'nu':-1.}])
def test_invalid_settings(changes):
    settings=dict(interval=.01,substeps=1,method='heun',box=6.,nu=.01)
    settings.update(changes)
    with pytest.raises(ValueError):evolve(np.zeros((3,4,4,4)),**settings)


def test_invalid_and_unprojected_fields():
    bad=np.zeros((3,8,8,8));bad[0]=np.sin(np.arange(8)*2*np.pi/8)[:,None,None]
    for u in (np.zeros((3,4,4,5)),np.ones((3,4,4,4),complex),
              np.full((3,4,4,4),np.nan),bad):
        with pytest.raises(ValueError):evolve(u,.01,1,'heun')
