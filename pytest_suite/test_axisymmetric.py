import math

import pytest

from axisymmetric import (
    cylindrical_divergence,
    meridional_velocity,
    phi,
    psi,
    swirl,
    velocity_from_psi,
)

R0 = 1.5
ALPHA = 0.7


def test_phi_peaks_at_bump_center():
    assert phi(R0, 0.0, R0, ALPHA) == 1.0
    assert phi(R0 + ALPHA, 0.0, R0, ALPHA) == pytest.approx(math.exp(-1.0))


def test_psi_vanishes_on_axis_and_midplane():
    assert psi(0.0, 0.4, R0, ALPHA) == 0.0
    assert psi(1.0, 0.0, R0, ALPHA) == 0.0


def test_analytic_matches_finite_difference_velocity():
    for r, z in [(1.0, 0.3), (1.8, -0.2), (0.9, 0.6)]:
        analytic = meridional_velocity(r, z, R0, ALPHA)
        numeric = velocity_from_psi(r, z, R0, ALPHA)
        assert analytic[0] == pytest.approx(numeric[0], abs=1e-10, rel=1e-7)
        assert analytic[1] == pytest.approx(numeric[1], abs=1e-10, rel=1e-7)


def test_meridional_field_is_divergence_free():
    for r, z in [(1.0, 0.3), (1.8, -0.2), (0.9, 0.6)]:
        assert cylindrical_divergence(r, z, R0, ALPHA) == pytest.approx(0.0, abs=1e-8)


def test_swirl_is_r_times_phi():
    r, z = 1.2, -0.4
    assert swirl(r, z, R0, ALPHA) == r * phi(r, z, R0, ALPHA)


def test_field_decays_far_from_bump():
    u_r, u_z = meridional_velocity(8.0, 8.0, R0, ALPHA)
    assert abs(u_r) < 1e-20
    assert abs(u_z) < 1e-20
    assert swirl(8.0, 8.0, R0, ALPHA) < 1e-20
