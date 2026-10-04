"""Stokes start on the cube, stepped with pressure and without."""

import math

from navier import PurePythonNavierStokes3D
from stokes import meridional_velocity, swirl

R0 = 1.5
ALPHA = 1.0


def sample_stokes(sim, dx):
    mid = (sim.N - 1) / 2.0
    for i in range(sim.N):
        x = (i - mid) * dx
        for j in range(sim.N):
            y = (j - mid) * dx
            for k in range(sim.N):
                z = (k - mid) * dx
                r = math.hypot(x, y)
                if r < 1e-8:
                    u_r, u_z = meridional_velocity(1e-8, z, R0, ALPHA)
                    u_th = swirl(1e-8, z, R0, ALPHA)
                    cos_t, sin_t = 1.0, 0.0
                else:
                    u_r, u_z = meridional_velocity(r, z, R0, ALPHA)
                    u_th = swirl(r, z, R0, ALPHA)
                    cos_t, sin_t = x / r, y / r
                sim.u[i][j][k] = u_r * cos_t - u_th * sin_t
                sim.v[i][j][k] = u_r * sin_t + u_th * cos_t
                sim.w[i][j][k] = u_z
                sim.S[i][j][k] = 0.0


def max_div(sim):
    worst = 0.0
    for i in range(sim.N):
        for j in range(sim.N):
            for k in range(sim.N):
                worst = max(worst, abs(sim.divergence_at(i, j, k)))
    return worst


def max_grad(sim):
    worst = 0.0
    for i in range(sim.N):
        for j in range(sim.N):
            for k in range(sim.N):
                ip, im, jp, jm, kp, km = sim.get_neighbors(i, j, k)
                comps = (
                    sim.u[ip][j][k] - sim.u[im][j][k],
                    sim.u[i][jp][k] - sim.u[i][jm][k],
                    sim.u[i][j][kp] - sim.u[i][j][km],
                    sim.v[ip][j][k] - sim.v[im][j][k],
                    sim.v[i][jp][k] - sim.v[i][jm][k],
                    sim.v[i][j][kp] - sim.v[i][j][km],
                    sim.w[ip][j][k] - sim.w[im][j][k],
                    sim.w[i][jp][k] - sim.w[i][jm][k],
                    sim.w[i][j][kp] - sim.w[i][j][km],
                )
                mag = math.sqrt(sum((c / (2.0 * sim.dx)) ** 2 for c in comps))
                worst = max(worst, mag)
    return worst


def test_meridional_velocity_matches_stream_function():
    u_r, u_z = meridional_velocity(1.5, 0.5, 1.5, 1.0)
    e = math.exp(-0.25)
    assert round(u_r, 5) == round(-1.5 * 0.5 * e, 5)
    assert round(u_z, 5) == round(1.0 * e, 5)
    assert round(swirl(1.5, 0.5, 1.5, 1.0), 5) == round(1.5 * e, 5)


def test_stokes_start_metrics():
    sim = PurePythonNavierStokes3D(N=8, dx=0.75)
    sample_stokes(sim, 0.75)
    assert round(max_div(sim), 5) == 0.46351
    assert round(max_grad(sim), 5) == 2.02332
    assert sim.mean_S() == 0.0


def test_pressure_projection_reduces_divergence_each_step():
    dx = 0.75
    bare = PurePythonNavierStokes3D(N=8, dx=dx)
    held = PurePythonNavierStokes3D(N=8, dx=dx)
    sample_stokes(bare, dx)
    sample_stokes(held, dx)
    assert round(max_grad(bare), 5) == round(max_grad(held), 5)
    for _ in range(4):
        bare.step(dt=0.02, P_U=0.0)
        held.step_with_pressure(dt=0.02, P_U=0.0)
        assert max_div(held) < max_div(bare)
    assert round(max_div(held), 5) == 0.09750
    assert round(max_grad(held), 5) == 2.04675
