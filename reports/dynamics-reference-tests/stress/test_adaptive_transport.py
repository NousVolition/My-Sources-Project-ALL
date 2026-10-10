"""Independent analytic checks for nonautonomous transport and changing steps."""
from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adaptive_transport import MODE, T, diagnostics, exact, field, manual_integrate, rk4_step, rotation, speed


@pytest.mark.parametrize("profile", ["slow-fast-slow", "fast-slow-fast"])
def test_rotation_against_independent_quadrature(profile):
    for t in (.13, 2.7, 5.91, T):
        integrated, _ = quad(lambda s: speed(s, profile), 0, t, epsabs=1e-12)
        assert abs(rotation(t, profile)-integrated) < 1e-12


@pytest.mark.parametrize("grid", [64, 128])
def test_resolved_fft_derivative_against_analytic_sine(grid):
    theta, rhs = field(grid, "slow-fast-slow")
    q = np.cos(MODE*theta)
    expected = MODE*speed(3.1, "slow-fast-slow")*np.sin(MODE*theta)
    np.testing.assert_allclose(rhs(3.1, q), expected, rtol=0, atol=3e-12)


def test_aliased_derivative_matches_wrong_mode():
    theta, rhs = field(32, "slow-fast-slow")
    q = np.cos(MODE*theta)
    np.testing.assert_allclose(q, np.cos(15*theta), atol=4e-14)
    wrong = 15*speed(3.1, "slow-fast-slow")*np.sin(15*theta)
    np.testing.assert_allclose(rhs(3.1, q), wrong, atol=3e-12)
    true = MODE*speed(3.1, "slow-fast-slow")*np.sin(MODE*theta)
    assert np.max(np.abs(rhs(3.1, q)-true)) > 50


def test_nonautonomous_rk4_convergence():
    # An independent ODE, q'=t*q, verifies the stage-time implementation.
    errors = []
    for n in (10, 20, 40):
        q = np.array([1.0])
        for j in range(n):
            q = rk4_step(lambda t, y: t*y, j/n, q, 1/n)
        errors.append(abs(q[0]-np.exp(.5)))
    assert all(12 < a/b < 20 for a, b in zip(errors[:-1], errors[1:]))


def test_manual_schedules_equal_work_and_exact_switch_boundaries():
    for order in ("large-small-large", "small-large-small"):
        _, times, _ = manual_integrate(64, "slow-fast-slow", order)
        assert len(times)-1 == 4500
        assert times[-1] == T
        assert 2.0 in times and 6.0 in times
        steps = np.diff(times)
        np.testing.assert_allclose(np.sort(np.unique(np.round(steps, 10))), [.001, .008])


@pytest.mark.parametrize("profile", ["slow-fast-slow", "fast-slow-fast"])
def test_tight_adaptive_reference_and_repeat(profile):
    theta, rhs = field(64, profile)
    inputs = dict(method="DOP853", rtol=1e-9, atol=1e-12, max_step=.1, dense_output=True)
    first = solve_ivp(rhs, (0, T), np.cos(MODE*theta), **inputs)
    second = solve_ivp(rhs, (0, T), np.cos(MODE*theta), **inputs)
    assert first.success and second.success
    np.testing.assert_array_equal(first.t, second.t)
    np.testing.assert_array_equal(first.y, second.y)
    times = np.linspace(0, T, 1601)
    result, _ = diagnostics(theta, times, first.sol(times).T, profile)
    assert result["accuracy_passed"] and result["conservation_passed"]
    assert result["boundedness_guard_passed"]


def test_exact_solution_satisfies_pde_independently():
    theta = np.array([.217, 1.543, 3.118])
    t, h = 2.173, 1e-6
    for profile in ("slow-fast-slow", "fast-slow-fast"):
        dt = (exact(theta, np.array([t+h]), profile)-exact(theta, np.array([t-h]), profile))/(2*h)
        dx = (exact(theta+h, np.array([t]), profile)-exact(theta-h, np.array([t]), profile))/(2*h)
        assert np.max(np.abs(dt+speed(t, profile)*dx)) < 2e-7
