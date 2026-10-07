"""Known energy identities and independent-loop checks for one-step probes."""
import numpy as np
import pytest

from hug_step_energy import step_energy
from navier import PurePythonNavierStokes3D


def taylor_green(n=12):
    x,y,z=np.meshgrid(*([np.arange(n)*2*np.pi/n]*3),indexing='ij')
    return np.array([np.sin(x)*np.cos(y),-np.cos(x)*np.sin(y),np.zeros_like(z)])


@pytest.mark.parametrize('amplitude',[0.,2.])
def test_constant_field_has_no_energy_change(amplitude):
    u=np.full((3,6,6,6),amplitude)
    before=u.copy()
    m,v=step_energy(u,.01,return_velocity=True)
    np.testing.assert_array_equal(u,before)
    np.testing.assert_array_equal(v,u)
    for name in ('energy_change','transport_change','viscous_change','quadratic_euler_change','projection_change'):
        assert m[name]==0.


@pytest.mark.parametrize('n',[8,10])
def test_shear_matches_exact_forward_euler_decay(n):
    box,dt,nu=4.,.003,.07
    y=np.arange(n)*2*np.pi/n
    u=np.zeros((3,n,n,n));u[0]=np.sin(y)[None,:,None]
    eigenvalue=(2*np.sin(np.pi/n)/(box/n))**2
    initial=box**3/4
    m,v=step_energy(u,dt,box,nu,True)
    np.testing.assert_allclose(v,(1-dt*nu*eigenvalue)*u,atol=2e-15)
    assert m['energy_after']==pytest.approx(initial*(1-dt*nu*eigenvalue)**2,abs=2e-13)
    assert m['viscous_change']==pytest.approx(-2*dt*nu*eigenvalue*initial)
    assert m['quadratic_euler_change']==pytest.approx(dt**2*(nu*eigenvalue)**2*initial)
    assert abs(m['transport_change'])<1e-15
    assert abs(m['projection_change'])<1e-15


def test_pressure_removes_the_entire_quadratic_term_for_steady_inviscid_flow():
    u=taylor_green()
    m,v=step_energy(u,.01,box=2*np.pi,nu=0.,return_velocity=True)
    np.testing.assert_allclose(v,u,atol=2e-15)
    assert m['quadratic_euler_change']>1e-4
    assert m['projection_change']==pytest.approx(-m['quadratic_euler_change'],abs=1e-13)
    assert abs(m['transport_change'])<1e-13
    assert abs(m['energy_change'])<1e-12


def test_random_field_matches_original_loop_step_and_does_not_mutate_input():
    n,dt,box,nu=6,.002,3.,.02
    sim=PurePythonNavierStokes3D(n,box/n)
    rng=np.random.default_rng(102)
    sim.u,sim.v,sim.w=rng.normal(size=(3,n,n,n))
    sim.project(backend='fft')
    u=np.array([sim.u,sim.v,sim.w]);before=u.copy()
    sim.nu,sim.sigma=nu,0.
    sim.step_with_pressure(dt,0.,backend='fft',step_backend='python')
    metrics,got=step_energy(u,dt,box,nu,True)
    np.testing.assert_allclose(got,np.array([sim.u,sim.v,sim.w]),rtol=2e-14,atol=2e-14)
    np.testing.assert_array_equal(u,before)
    assert abs(metrics['ledger_residual'])<1e-11
    assert metrics['after_max_div']<1e-12


def test_linear_and_quadratic_terms_scale_with_dt_on_identical_input():
    u=taylor_green()
    full=step_energy(u,.01,2*np.pi,.02)
    half=step_energy(u,.005,2*np.pi,.02)
    assert full['viscous_change']==pytest.approx(2*half['viscous_change'])
    assert full['quadratic_euler_change']==pytest.approx(4*half['quadratic_euler_change'])
    assert full['projection_change']==pytest.approx(4*half['projection_change'],abs=1e-13)


def test_projection_accounts_for_initial_gradient_energy():
    n=12
    u=np.zeros((3,n,n,n));u[0]=np.cos(np.arange(n)*2*np.pi/n)[:,None,None]
    m,v=step_energy(u,.01,2*np.pi,0.,True)
    assert m['initial_max_div']>.5
    assert np.max(np.abs(v))<1e-14
    assert m['projection_change']==pytest.approx(-m['energy_predictor'],abs=2e-12)


@pytest.mark.parametrize('dt,box,nu',[(0,6,.01),(-.1,6,.01),(.1,0,.01),(.1,6,-.01),(.1,6,np.nan)])
def test_invalid_settings_are_rejected(dt,box,nu):
    with pytest.raises(ValueError):
        step_energy(np.zeros((3,4,4,4)),dt,box,nu)


def test_invalid_arrays_are_rejected():
    for u in (np.zeros((4,4,4)),np.zeros((3,4,4,5)),np.ones((3,4,4,4),dtype=complex),np.full((3,4,4,4),np.nan)):
        with pytest.raises(ValueError):
            step_energy(u,.01)
