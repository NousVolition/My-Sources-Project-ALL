"""Independent checks of the screenshot's formula and Cartesian dynamics."""
from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from limit_cycle import exact_radius, exact_state, rhs, integrate, metrics


@pytest.mark.parametrize("r0", [0, .2, .8, 1, 1.5, 2])
def test_formula_initial_value_and_independent_numerical_solution(r0):
    assert np.isclose(exact_radius(0, r0), r0, atol=1e-14)
    t = np.linspace(0, 20, 501)
    initial = r0*np.array([np.cos(.3), np.sin(.3)])
    ref = solve_ivp(lambda t, q: rhs(q), (0, 20), initial, method="DOP853", rtol=1e-12, atol=1e-13, t_eval=t)
    assert ref.success
    np.testing.assert_allclose(ref.y.T, exact_state(t, r0, .3), atol=3e-10, rtol=0)


def test_cartesian_polar_identities_independent():
    for r in (.1, .7, 1, 2.3):
        for theta in (-2.2, .37, 2):
            x, y = r*np.array([np.cos(theta), np.sin(theta)])
            dx, dy = rhs([x, y])
            assert np.isclose((x*dx+y*dy)/r, (1-r*r)*r, atol=1e-13)
            assert np.isclose((x*dy-y*dx)/(r*r), 1, atol=1e-13)


def test_radial_formula_residual_and_semigroup():
    for r0 in (.2, .8, 1.5, 2):
        t, h = .71, 1e-5
        derivative = (exact_radius(t+h, r0)-exact_radius(t-h, r0))/(2*h)
        value = exact_radius(t, r0)
        assert abs(derivative-(1-value*value)*value) < 2e-9
        np.testing.assert_allclose(exact_radius(.8+.3, r0), exact_radius(.3, exact_radius(.8, r0)), atol=1e-14)


def test_origin_exception_cycle_period_and_radial_direction():
    t = np.linspace(0, 20, 401)
    np.testing.assert_array_equal(exact_radius(t, 0), np.zeros_like(t))
    np.testing.assert_array_equal(exact_radius(t, 1), np.ones_like(t))
    assert np.all(np.diff(exact_radius(t, .2)) >= 0)
    assert np.all(np.diff(exact_radius(t, 2)) <= 0)
    np.testing.assert_allclose(exact_state(2*np.pi, 1, .3), exact_state(0, 1, .3), atol=1e-14)


def test_rk4_fourth_order_against_exact_solution():
    errors = []
    for h in (.05, .025, .0125):
        t, states = integrate(rhs, exact_state(0, .2, .3), 4, h, "rk4")
        errors.append(np.max(np.linalg.norm(states-exact_state(t, .2, .3), axis=1)))
    assert all(14 < a/b < 18 for a, b in zip(errors[:-1], errors[1:]))


def test_origin_has_no_assigned_phase_or_conserved_radius_claim():
    t, states = integrate(rhs, [0, 0], 1, .1)
    result, arrays = metrics(t, states, 0, .3)
    assert result["phase_accuracy_passed"] is None
    assert not result["phase_applicable"] and "phase_error_degrees" not in arrays
    assert exact_radius(1, .2) > .2  # Radius is not a conserved quantity away from the cycle.


def test_formula_domain_is_explicit():
    with pytest.raises(ValueError):
        exact_radius(-1, 2)
    with pytest.raises(ValueError):
        exact_radius(0, -1)
