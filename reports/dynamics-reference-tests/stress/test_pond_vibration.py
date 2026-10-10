"""Independent model identities, reference checks, and pulse-scaling tests."""
from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pond_vibration as p


@pytest.mark.parametrize("width", [.2, 2.0])
def test_bed_velocity_is_independent_derivative(width):
    t = p.START+width*np.array([.1, .3, .5, .8])
    h = width*1e-6
    b_plus, _ = p.bed(t+h, 1e-4, width)
    b_minus, _ = p.bed(t-h, 1e-4, width)
    _, db = p.bed(t, 1e-4, width)
    np.testing.assert_allclose((b_plus-b_minus)/(2*h), db, rtol=1e-8, atol=1e-12)


def test_modal_equations_satisfy_shallow_water_pde():
    x = np.linspace(0, p.L, 41)
    a, u, t = 2e-5, -3e-5, p.START+.073
    deriv = p.rhs(t, [a, u, 0, 0], 1e-4, .2)
    _, db = p.bed(t, 1e-4, .2)
    np.testing.assert_allclose(deriv[0]*np.cos(p.K*x)+p.H*p.K*u*np.cos(p.K*x), db*np.cos(p.K*x), atol=1e-18)
    np.testing.assert_allclose(deriv[1]*np.sin(p.K*x)-p.G*p.K*a*np.sin(p.K*x), -p.DRAG*u*np.sin(p.K*x), atol=1e-18)
    assert abs(u*np.sin(p.K*p.L)) < 1e-18


def test_energy_work_and_dissipation_identity():
    state = np.array([.0003, -.0002, 0, 0])
    deriv = p.rhs(p.START+.31, state, .0001, 2)
    dE = p.RHO*p.L/2*(p.G*state[0]*deriv[0]+p.H*state[1]*deriv[1])
    assert abs(dE-(deriv[2]-deriv[3])) < 1e-15
    assert deriv[3] >= 0


@pytest.mark.parametrize("width", [.2, 2.0])
def test_convolution_against_independent_adaptive_ode(width):
    # Bound the maximum interval so the solver cannot step over an unseen short pulse.
    t = np.linspace(0, 8, 1001)
    solved = solve_ivp(lambda t, y: p.rhs(t, y, 1e-4, width), (0, 8), np.zeros(4),
        method="DOP853", rtol=1e-11, atol=1e-14, max_step=width/20, t_eval=t)
    assert solved.success
    np.testing.assert_allclose(solved.y[:2].T, p.reference(t, 1e-4, width), rtol=0, atol=2e-11)


def test_strength_scaling_and_zero_input():
    t = np.linspace(0, 10, 701)
    np.testing.assert_allclose(p.reference(t, 1e-4, .2), 10*p.reference(t, 1e-5, .2), atol=1e-18)
    _, state = p.integrate(0, .2, .05)
    np.testing.assert_array_equal(state, np.zeros_like(state))


def test_undisturbed_depth_and_volume_have_no_added_water():
    x = np.linspace(0, p.L, 65)
    t = np.linspace(0, 20, 401)
    state = p.reference(t, 1e-4, .2)
    b, _ = p.bed(t, 1e-4, .2)
    depth_change = (state[:, 0]-b)[:, None]*np.cos(p.K*x)[None, :]
    assert np.abs(np.trapezoid(depth_change, x, axis=1)).max() < 1e-15
    assert np.min(p.H+depth_change) > 0
