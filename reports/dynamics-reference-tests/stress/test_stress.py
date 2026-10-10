"""Independent checks of the new benchmark and batch-scan implementation."""
from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.integrate import solve_ivp

sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_stress import PLAN, lorenz_batch, transport_exact
from core import integrate, lorenz


def test_batched_equations_equal_individual_models():
    states=np.random.default_rng(71).normal(size=(20,3))*10
    rhos=np.linspace(.9,50,20)
    expected=np.array([lorenz(y,rho) for y,rho in zip(states,rhos)])
    np.testing.assert_allclose(lorenz_batch(states,rhos),expected,rtol=0,atol=0)


@pytest.mark.parametrize("rho", [0.9,10.,24.746842105263154,50.])
def test_batch_trajectory_against_independent_adaptive_solver(rho):
    initial=np.array([1.,2.,3.])
    # integrate() expects a flat state; reshape only for the batch RHS.
    t,y=integrate(lambda a:lorenz_batch(a.reshape(1,3),np.array([rho])).ravel(),initial,2.,.00125)
    ref=solve_ivp(lambda t,y:lorenz(y,rho),(0,2),initial,method="DOP853",t_eval=t,rtol=2e-13,atol=2e-14)
    assert ref.success
    assert np.max(np.linalg.norm(y-ref.y.T,axis=1))<1e-5


def test_transport_benchmark_pde_with_independent_finite_differences():
    # Avoid reusing the FFT derivative or the benchmark's closed-form time derivative.
    theta=np.random.default_rng(44).uniform(0,2*np.pi,32)
    h=1e-5
    t=.37
    angular=(transport_exact(theta-2*h,t)-8*transport_exact(theta-h,t)
             +8*transport_exact(theta+h,t)-transport_exact(theta+2*h,t))/(12*h)
    settings=PLAN["transport"]
    source=np.full(len(theta),settings["forcing_mean"])
    for k,amplitude in settings["forcing_cosine_amplitudes"].items():
        source+=amplitude*np.cos(int(k)*theta)
    rhs=source-settings["leak"]*transport_exact(theta,t)-settings["omega"]*angular
    errors=[]
    for time_step in (1e-5,5e-6,2.5e-6):
        temporal=(transport_exact(theta,t+time_step)-transport_exact(theta,t-time_step))/(2*time_step)
        errors.append(np.max(abs(temporal-rhs)))
    # The original 1e-5 difference had 1.24e-7 truncation error. Refine the
    # independent derivative, retaining the original acceptance tolerance.
    assert errors[-1]<1e-7
    assert errors[0]/errors[-1]>8


def test_transport_initial_condition_and_positive_source():
    theta=np.linspace(0,2*np.pi,8192,endpoint=False)
    np.testing.assert_allclose(transport_exact(theta,0.),60.,rtol=0,atol=0)
    settings=PLAN["transport"]
    # A lower bound by the triangle inequality is valid for every angle.
    lower_bound=settings["forcing_mean"]-sum(settings["forcing_cosine_amplitudes"].values())
    assert lower_bound==14.


@pytest.mark.parametrize("grid,aliased_mode", [(16,1),(32,15)])
def test_intended_alias_is_present(grid,aliased_mode):
    theta=2*np.pi*np.arange(grid)/grid
    np.testing.assert_allclose(np.cos(17*theta),np.cos(aliased_mode*theta),rtol=0,atol=2e-14)


def test_long_plan_has_exact_integer_step_and_block_alignment():
    long=PLAN["long_run"]
    for h in long["steps"]:
        for duration in [long["transient"],long["duration"],long["qr_interval"],long["block"]]:
            assert abs(duration/h-round(duration/h))<1e-8
    for duration in long["prefix_durations"]:
        assert abs(duration/long["block"]-round(duration/long["block"]))<1e-8
