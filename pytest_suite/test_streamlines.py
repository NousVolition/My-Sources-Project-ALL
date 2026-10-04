from streamlines import ALPHA, R, R0, Z, meridional_velocity, phi, psi


def test_arrows_follow_curves_of_constant_psi():
    h = 1e-6
    u_r, u_z = meridional_velocity(R, Z, R0, ALPHA)
    dpsi_dr = (psi(R + h, Z, R0, ALPHA) - psi(R - h, Z, R0, ALPHA)) / (2.0 * h)
    dpsi_dz = (psi(R, Z + h, R0, ALPHA) - psi(R, Z - h, R0, ALPHA)) / (2.0 * h)
    along_the_arrow = u_r * dpsi_dr + u_z * dpsi_dz
    assert abs(along_the_arrow) < 1e-8


def test_the_two_cells_are_mirror_images():
    upper = meridional_velocity(R, Z, R0, ALPHA)
    lower = meridional_velocity(R, -Z, R0, ALPHA)
    assert abs(upper[0] - lower[0]) < 1e-12
    assert abs(upper[1] + lower[1]) < 1e-12


def test_the_midplane_is_the_dividing_line():
    for radius in (0.2, 1.0, R0, 3.0):
        assert psi(radius, 0.0, R0, ALPHA) == 0.0


def test_ur_scales_like_r_near_the_axis():
    ratios = []
    for radius in (1e-3, 1e-5, 1e-7):
        u_r, _ = meridional_velocity(radius, Z, R0, ALPHA)
        ratios.append(u_r / radius)
    assert abs(ratios[-1] - ratios[-2]) < 1e-5
    assert abs(ratios[-1]) > 1e-3


def test_without_r_squared_ur_blows_up_on_the_axis():
    """The r^2 in psi is what keeps u_r finite.

    Drop it and the axis fails.
    """
    def psi_without_r_squared(r, z):
        return z * phi(r, z, R0, ALPHA)

    h = 1e-6
    radius = 1e-4
    dpsi_dz = (
        psi_without_r_squared(radius, Z + h)
        - psi_without_r_squared(radius, Z - h)
    ) / (2.0 * h)
    blown_up = abs(-(1.0 / radius) * dpsi_dz)
    kept = abs(meridional_velocity(radius, Z, R0, ALPHA)[0])
    assert blown_up > 100.0
    assert kept < 1e-4


def test_bump_is_centered_on_r0():
    assert phi(R0, 0.0, R0, ALPHA) > phi(R0 + 0.5, 0.0, R0, ALPHA)
    assert phi(R0, 0.0, R0, ALPHA) > phi(R0, 0.5, R0, ALPHA)
    left = phi(R0 - 0.25, 0.0, R0, ALPHA)
    right = phi(R0 + 0.25, 0.0, R0, ALPHA)
    assert abs(left - right) < 1e-12


def test_inner_and_outer_z_send_ur_opposite_ways():
    inside = meridional_velocity(R, 0.0, R0, ALPHA)[0]
    outside = meridional_velocity(R, 1.0, R0, ALPHA)[0]
    assert inside < 0.0
    assert outside > 0.0
