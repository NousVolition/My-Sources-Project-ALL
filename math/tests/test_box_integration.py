"""Checks recovered from PR #1, using the shared sampler and diagnostics."""
import math
from navier import PurePythonNavierStokes3D
from stokes import meridional_velocity, swirl
from box_experiment import sample_stokes, diagnostics

def max_div(sim):
    return diagnostics(sim)["max_div"]

def max_grad(sim):
    return diagnostics(sim)["max_grad"]

def test_meridional_velocity_matches_stream_function():
    u_r, u_z = meridional_velocity(1.5, 0.5, 1.5, 1.0)
    e = math.exp(-0.25)
    assert round(u_r, 5) == round(-1.5 * 0.5 * e, 5)
    assert round(u_z, 5) == round(1.0 * e, 5)
    assert round(swirl(1.5, 0.5, 1.5, 1.0), 5) == round(1.5 * e, 5)


def test_stokes_start_metrics():
    sim = PurePythonNavierStokes3D(N=8, dx=0.75)
    sample_stokes(sim)
    assert round(max_div(sim), 5) == 0.46351
    assert round(max_grad(sim), 5) == 2.02332
    assert sim.mean_S() == 0.0


def test_pressure_projection_reduces_divergence_each_step():
    dx = 0.75
    bare = PurePythonNavierStokes3D(N=8, dx=dx)
    held = PurePythonNavierStokes3D(N=8, dx=dx)
    sample_stokes(bare)
    sample_stokes(held)
    assert round(max_grad(bare), 5) == round(max_grad(held), 5)
    for _ in range(4):
        bare.step(dt=0.02, P_U=0.0)
        held.step_with_pressure(dt=0.02, P_U=0.0, backend="fft")
        assert max_div(held) < max_div(bare)
    # The previous residual fixture described the incompatible Poisson stencil.
    assert max_div(held) < 1e-10


R0, ALPHA = 1.5, 1.0


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
    bare = PurePythonNavierStokes3D(N=16, dx=18.0 / 16)
    held = PurePythonNavierStokes3D(N=16, dx=18.0 / 16)
    sample_stokes(bare)
    sample_stokes(held)
    bare.step(dt=0.02, P_U=0.0)
    held.step_with_pressure(dt=0.02, P_U=0.0, backend="fft")
    assert max_div(held) < max_div(bare)
    assert max_grad(held) < 10.0
