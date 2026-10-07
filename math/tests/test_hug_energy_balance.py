"""Analytic controls for the saved-flow energy budget, without long evolution."""
import numpy as np
import pytest

from hug_energy_balance import dissipation_rate, energy_budget
from hug_resolution import grid_spectrum


def test_viscous_rate_matches_independent_laplacian_work():
    rng=np.random.default_rng(57)
    velocity=rng.normal(size=(3,8,8,8))
    box,nu=3.5,.07
    dx=box/8
    lap=sum((np.roll(velocity,-1,axis)+np.roll(velocity,1,axis)-2*velocity)/dx**2
            for axis in (1,2,3))
    work=-nu*dx**3*np.sum(velocity*lap)
    metrics,_=grid_spectrum(velocity,box)
    rate=dissipation_rate(metrics['neighbor_gradient_rms'],box,nu)
    assert rate == pytest.approx(work,rel=2e-14)


def test_linear_loss_has_exact_budget_on_uneven_times():
    t=np.array([0.,.03,.1,.27,.4])
    before=t.copy()
    e=8.-2.*t
    g=np.full_like(t,np.sqrt(2./(.1*3.**3)))
    result=energy_budget(t,e,g,box=3.,nu=.1)
    for key in ('viscous_trapezoid','viscous_simpson','viscous_coarsened_trapezoid'):
        assert result['summary'][key] == pytest.approx(.8,abs=2e-15)
    assert abs(result['summary']['residual_simpson']) < 2e-15
    np.testing.assert_array_equal(t,before)


def test_quadratic_loss_rate_uses_actual_uneven_times():
    t=np.array([0.,.02,.08,.16,.4])
    q=1.+2*t+3*t*t
    e=5.-(t+t*t+t**3)
    result=energy_budget(t,e,np.sqrt(q),box=1.,nu=1.)['summary']
    assert result['viscous_simpson'] == pytest.approx(.624,abs=2e-15)
    assert abs(result['residual_simpson']) < 2e-15
    assert result['quadrature_method_difference'] > 0
    assert result['coarsening_difference'] > 0


def test_exact_decaying_shear_exposes_only_quadrature_error():
    # Spatially discretized sine shear decays at the Laplacian eigenvalue.
    n,box,nu,mode=16,4.,.03,2
    eigenvalue=(2*np.sin(np.pi*mode/n)/(box/n))**2
    exact_loss=box**3/4*(1-np.exp(-2*nu*eigenvalue*.4))
    errors=[]
    for count in (5,9):
        t=np.linspace(0,.4,count)
        amplitude=np.exp(-nu*eigenvalue*t)
        e=box**3/4*amplitude**2
        g=np.sqrt(eigenvalue/2)*amplitude
        s=energy_budget(t,e,g,box,nu)['summary']
        assert s['observed_energy_loss'] == pytest.approx(exact_loss)
        errors.append(abs(s['viscous_simpson']-exact_loss))
    assert errors[1] < errors[0]/10
    # Composite Simpson error <= T*h^4*max|q''''|/180 for this exponential.
    bound=.4/180*(.4/8)**4*(2*nu*eigenvalue)**5*(box**3/4)
    assert errors[1] <= bound+1e-13


@pytest.mark.parametrize('observed_loss,sign',[(.2,1),(.6,-1)])
def test_residual_sign_is_not_relabelled_as_loss(observed_loss,sign):
    t=np.array([0.,.2,.4])
    s=energy_budget(t,2-observed_loss*t/.4,np.ones(3),box=1.,nu=1.)['summary']
    assert np.sign(s['residual_simpson']) == sign
    assert s['simpson_gap_percent_loss'] == pytest.approx(100*(.4-observed_loss)/observed_loss)


def test_constant_velocity_has_no_viscous_loss_or_percentage():
    s=energy_budget([0,.1,.2],[2,2,2],[0,0,0])['summary']
    assert s['viscous_simpson'] == s['residual_simpson'] == 0
    assert s['simpson_gap_percent_loss'] is None


@pytest.mark.parametrize('times,energies,gradients',[
    ([0,.1],[1,1],[1,1]),
    ([0,.1,.1],[1,1,1],[1,1,1]),
    ([0,.2,.1],[1,1,1],[1,1,1]),
    ([0,.1,.2],[1,-1,1],[1,1,1]),
    ([0,.1,.2],[1,1,1],[1,-1,1]),
    ([0,.1,.2],[1,1,1],[1,np.nan,1]),
    ([0,.1,.2],[1,np.inf,1],[1,1,1]),
])
def test_invalid_observations_are_rejected(times,energies,gradients):
    with pytest.raises(ValueError):
        energy_budget(times,energies,gradients)


def test_invalid_physical_scales_are_rejected():
    for box,nu in [(0,.01),(6,-.01),(np.inf,.01),(6,np.nan)]:
        with pytest.raises(ValueError):
            dissipation_rate(1.,box,nu)
