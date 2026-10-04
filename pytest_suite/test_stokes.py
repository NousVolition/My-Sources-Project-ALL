import math

from navier import PurePythonNavierStokes3D
from stokes import meridional_velocity, swirl
from stokes_box import ALPHA, BOX, R0, max_div, max_grad, sample_stokes


def test_meridional_velocity_on_axis():
    for z in (-0.5, 0.0, 0.7):
        u_r, u_z = meridional_velocity(1e-8, z, R0, ALPHA)
        assert abs(u_r) < 1e-7
        assert math.isclose(u_z, 2.0 * z * math.exp(-(R0 ** 2 + z ** 2)), rel_tol=1e-6)


def test_swirl_matches_radial_profile():
    envelope = math.exp(-((1.2 - R0) ** 2 + 0.3 ** 2) / (ALPHA ** 2))
    u_r, u_z = meridional_velocity(1.2, 0.3, R0, ALPHA)
    assert swirl(1.2, 0.3, R0, ALPHA) == 1.2 * envelope
    assert u_r == -1.2 * (1.0 - 2.0 * 0.3 ** 2) * envelope
    assert u_z == 0.3 * (2.0 - 2.0 * 1.2 * (1.2 - R0)) * envelope


def test_divergence_of_uniform_field_is_zero():
    sim = PurePythonNavierStokes3D(N=4)
    for i in range(sim.N):
        for j in range(sim.N):
            for k in range(sim.N):
                assert sim.divergence_at(i, j, k) == 0.0


def test_divergence_of_radial_field():
    sim = PurePythonNavierStokes3D(N=4, dx=0.5)
    mid = 1.5
    for i in range(sim.N):
        for j in range(sim.N):
            for k in range(sim.N):
                sim.u[i][j][k] = i - mid
                sim.v[i][j][k] = j - mid
                sim.w[i][j][k] = k - mid
    for i in range(1, sim.N - 1):
        for j in range(1, sim.N - 1):
            for k in range(1, sim.N - 1):
                assert sim.divergence_at(i, j, k) == 6.0
    assert sim.divergence_at(0, 0, 0) == -6.0


def test_step_with_pressure_reduces_divergence():
    bare = PurePythonNavierStokes3D(N=16, dx=BOX / 16)
    held = PurePythonNavierStokes3D(N=16, dx=BOX / 16)
    sample_stokes(bare, BOX / 16)
    sample_stokes(held, BOX / 16)
    bare.step(dt=0.02, P_U=0.0)
    held.step_with_pressure(dt=0.02, P_U=0.0)
    assert max_div(held) < max_div(bare)
    assert max_grad(held) < 10.0
