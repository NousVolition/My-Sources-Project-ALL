"""Regressions for grid spacing, convergence, and centered-grid blind modes."""
import copy
import math

import pytest

from navier import PurePythonNavierStokes3D


def sinusoidal_cube(N, dx, quarter_wave=False):
    sim = PurePythonNavierStokes3D(N, dx)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                phase = math.pi / 2 * (i + j + k) if quarter_wave else 2 * math.pi * i / N
                sim.u[i][j][k] = math.sin(phase)
                sim.v[i][j][k] = 0.0
                sim.w[i][j][k] = 0.0
    return sim


@pytest.mark.parametrize("N", [7, 8])
@pytest.mark.parametrize("dx", [0.1875, 0.75, 1.125])
def test_projection_works_at_nonunit_spacing(N, dx):
    sim = sinusoidal_cube(N, dx)
    report = sim.project()
    assert report["before"] > 0.1
    assert report["after"] < 1e-9


def test_damping_removes_mode_that_undamped_jacobi_cannot_converge():
    sim = sinusoidal_cube(8, 0.375, quarter_wave=True)
    assert sim.project()["after"] < 1e-9


def test_nonconvergence_is_reported_without_changing_velocity():
    sim = sinusoidal_cube(8, 0.375)
    before = copy.deepcopy((sim.u, sim.v, sim.w))
    with pytest.raises(RuntimeError, match="did not converge"):
        sim.project(iterations=1, atol=1e-12, rtol=0.0)
    assert (sim.u, sim.v, sim.w) == before


def test_checkerboard_is_invisible_to_centered_divergence():
    sim = PurePythonNavierStokes3D(8, 0.375)
    for i in range(8):
        for j in range(8):
            for k in range(8):
                sim.u[i][j][k] = (-1.0)**i
                sim.v[i][j][k] = sim.w[i][j][k] = 0.0
    before = copy.deepcopy(sim.u)
    report = sim.project()
    assert report["before"] == 0.0
    assert report["after"] == 0.0
    assert sim.u == before


@pytest.mark.parametrize("N", [7, 8])
def test_fft_matches_jacobi_and_is_orthogonal_projection(N):
    np = pytest.importorskip("numpy")
    sim = PurePythonNavierStokes3D(N, 0.375)
    original = np.random.default_rng(192).normal(size=(3, N, N, N))
    sim.u, sim.v, sim.w = original.tolist()
    twin = copy.deepcopy(sim)
    report = sim.project(backend="fft")
    assert report["after"] < 1e-12
    twin.project()
    after = np.array((sim.u, sim.v, sim.w))
    np.testing.assert_allclose(after, (twin.u, twin.v, twin.w), atol=1e-9, rtol=0)
    np.testing.assert_allclose(after.mean(axis=(1, 2, 3)),
                               original.mean(axis=(1, 2, 3)), atol=1e-14)
    assert np.sum(after**2) <= np.sum(original**2)
    correction = original - after
    assert abs(np.sum(after * correction)) < 1e-10
    sim.project(backend="fft")
    np.testing.assert_allclose((sim.u, sim.v, sim.w), after, atol=1e-13, rtol=0)


def test_step_wrapper_projects_actual_updated_velocity():
    sim = sinusoidal_cube(8, 0.375)
    sim.sigma = 0.0
    sim.project()
    sim.u[1][2][3] += 0.1
    report = sim.step_with_pressure(0.02, 0.0)
    assert report["before"] > 0.01
    assert report["after"] < 1e-9
