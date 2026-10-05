"""The user's 20 original checks, retained for the original radial Gaussian.

Passing these checks does not establish Cartesian smoothness at r=0.
"""
import math

from navier import PurePythonNavierStokes3D
from stokes import (
    cylindrical_divergence, meridional_velocity, phi, psi,
    swirl, velocity_from_psi,
)

R0, ALPHA, R, Z = 1.5, 1.0, 1.2, 0.4


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
    assert abs(cylindrical_divergence(R, Z, R0, ALPHA)) < 1e-8


def test_arrows_follow_curves_of_constant_psi():
    h = 1e-6
    u_r, u_z = meridional_velocity(R, Z, R0, ALPHA)
    dpsi_dr = (psi(R + h, Z, R0, ALPHA) - psi(R - h, Z, R0, ALPHA)) / (2.0 * h)
    dpsi_dz = (psi(R, Z + h, R0, ALPHA) - psi(R, Z - h, R0, ALPHA)) / (2.0 * h)
    assert abs(u_r * dpsi_dr + u_z * dpsi_dz) < 1e-8


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
    def psi_without_r_squared(r, z):
        return z * phi(r, z, R0, ALPHA)
    h, radius = 1e-6, 1e-4
    dpsi_dz = (psi_without_r_squared(radius, Z + h)
               - psi_without_r_squared(radius, Z - h)) / (2.0 * h)
    blown_up = abs(-dpsi_dz / radius)
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
    assert meridional_velocity(R, 0.0, R0, ALPHA)[0] < 0.0
    assert meridional_velocity(R, 1.0, R0, ALPHA)[0] > 0.0


def _div_at(sim, i, j, k):
    ip, im, jp, jm, kp, km = sim.get_neighbors(i, j, k)
    du = (sim.u[ip][j][k] - sim.u[im][j][k]) / (2.0 * sim.dx)
    dv = (sim.v[i][jp][k] - sim.v[i][jm][k]) / (2.0 * sim.dx)
    dw = (sim.w[i][j][kp] - sim.w[i][j][km]) / (2.0 * sim.dx)
    return du + dv + dw


def test_one_checkpoint_is_not_the_whole_field():
    misses = 0
    for r in (0.4, 0.8, 1.2, 1.5, 2.2):
        for z in (-0.8, -0.2, 0.3, 0.9):
            if abs(cylindrical_divergence(r, z, R0, ALPHA)) > 1e-8:
                misses += 1
    assert misses == 0


def test_stream_function_has_no_time_and_no_pressure():
    assert not hasattr(psi, "step")
    speed = math.hypot(*meridional_velocity(1.2, 0.4, R0, ALPHA))
    assert speed > 0.0
    assert psi(1.2, 0.4, R0, ALPHA) != speed


def test_swirl_can_change_without_moving_the_curves():
    before = psi(1.2, 0.4, R0, ALPHA)
    quiet = swirl(1.2, 0.4, R0, ALPHA)
    loud = 7.0 * quiet
    after = psi(1.2, 0.4, R0, ALPHA)
    assert before == after
    assert loud != quiet


def test_cube_script_starts_flat_not_from_psi():
    sim = PurePythonNavierStokes3D(N=4)
    assert sim.u[0][0][0] == 0.05
    assert sim.v[1][2][3] == 0.01
    assert sim.w[2][2][2] == 0.02
    sample = meridional_velocity(1.2, 0.4, R0, ALPHA)
    assert (0.05, 0.01) != tuple(round(v, 2) for v in sample)


def test_below_the_threshold_the_cube_does_not_kick_w():
    sim = PurePythonNavierStokes3D(N=4)
    before = sim.w[0][0][0]
    sim.step(dt=0.1, P_U=3.0)
    assert abs(sim.w[0][0][0] - before) < 1e-12


def test_the_cube_update_does_not_remove_divergence():
    sim = PurePythonNavierStokes3D(N=8)
    for i in range(sim.N):
        s = math.sin(2.0 * math.pi * i / sim.N)
        for j in range(sim.N):
            for k in range(sim.N):
                sim.u[i][j][k] = s
                sim.v[i][j][k] = 0.0
                sim.w[i][j][k] = 0.0
    before = abs(_div_at(sim, 0, 0, 0))
    sim.step(dt=0.01, P_U=0.0)
    after = abs(_div_at(sim, 0, 0, 0))
    assert before > 0.1
    assert after > before * 0.5
