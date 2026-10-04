from stokes import (
    cylindrical_divergence,
    meridional_velocity,
    phi,
    psi,
    swirl,
    velocity_from_psi,
)

R0 = 1.5
ALPHA = 1.0
R = 1.2
Z = 0.4


def test_checkpoint_matches_central_difference():
    analytic = meridional_velocity(R, Z, R0, ALPHA)
    numeric = velocity_from_psi(R, Z, R0, ALPHA)
    assert abs(analytic[0] - numeric[0]) < 1e-10
    assert abs(analytic[1] - numeric[1]) < 1e-10


def test_written_velocity_uses_the_same_bump():
    bump = phi(R, Z, R0, ALPHA)
    u_r, u_z = meridional_velocity(R, Z, R0, ALPHA)
    expected_r = -R * (1.0 - 2.0 * (Z ** 2) / (ALPHA ** 2)) * bump
    expected_z = Z * (2.0 - 2.0 * R * (R - R0) / (ALPHA ** 2)) * bump
    assert u_r == expected_r
    assert u_z == expected_z


def test_psi_is_zero_on_the_midplane_and_odd_in_z():
    assert psi(R, 0.0, R0, ALPHA) == 0.0
    assert psi(R, Z, R0, ALPHA) == -psi(R, -Z, R0, ALPHA)


def test_swirl_is_added_by_hand():
    assert swirl(R, Z, R0, ALPHA) == R * phi(R, Z, R0, ALPHA)
    u_r, u_z = meridional_velocity(R, Z, R0, ALPHA)
    assert swirl(R, Z, R0, ALPHA) not in (u_r, u_z)


def test_radial_velocity_stays_finite_on_the_axis():
    u_r, _ = meridional_velocity(1e-8, Z, R0, ALPHA)
    assert abs(u_r) < 1e-6


def test_bump_decays_away_from_its_center():
    near = abs(phi(R0, 0.0, R0, ALPHA))
    far = abs(phi(R0 + 4.0, 4.0, R0, ALPHA))
    assert near == 1.0
    assert far < 1e-6


def test_meridional_field_is_divergence_free():
    div = cylindrical_divergence(R, Z, R0, ALPHA)
    assert abs(div) < 1e-8
