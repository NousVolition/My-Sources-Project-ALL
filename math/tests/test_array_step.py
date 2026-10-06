"""Verify the accelerated stencils against the independent original loops."""
import copy

import numpy as np
import pytest

from box_experiment import sample_stokes
from navier import PurePythonNavierStokes3D


@pytest.mark.parametrize("points", [7, 8])
@pytest.mark.parametrize("with_force", [False, True])
def test_array_step_matches_original_loops(points, with_force):
    sim = PurePythonNavierStokes3D(points, .375)
    rng = np.random.default_rng(248)
    sim.u, sim.v, sim.w = rng.normal(size=(3, points, points, points)).tolist()
    sim.S = rng.normal(size=(points, points, points)).tolist()
    sim.e_k = [.2, -.3, 1.]
    sim.sigma = .5 if with_force else 0.
    twin = copy.deepcopy(sim)
    for _ in range(3):
        sim.step(.001, .4 if with_force else 0.)
        twin.step(.001, .4 if with_force else 0., backend="numpy")
    np.testing.assert_allclose((sim.u, sim.v, sim.w, sim.S),
                               (twin.u, twin.v, twin.w, twin.S), atol=2e-14, rtol=0)


def test_array_hug_evolution_matches_loop_evolution():
    sim = PurePythonNavierStokes3D(8, .75)
    sim.sigma = 0.
    sample_stokes(sim, "hug")
    sim.project(backend="fft")
    twin = copy.deepcopy(sim)
    for _ in range(4):
        sim.step_with_pressure(.002, 0., backend="fft")
        twin.step_with_pressure(.002, 0., backend="fft", step_backend="numpy")
    np.testing.assert_allclose((sim.u, sim.v, sim.w),
                               (twin.u, twin.v, twin.w), atol=2e-14, rtol=0)
    assert np.max(np.abs(twin.S)) == 0.


def test_array_projection_matches_independent_jacobi():
    sim = PurePythonNavierStokes3D(8, .375)
    sim.u, sim.v, sim.w = np.random.default_rng(18).normal(size=(3, 8, 8, 8))
    twin = copy.deepcopy(sim)
    sim.project(backend="fft")
    twin.u, twin.v, twin.w = [a.tolist() for a in (twin.u, twin.v, twin.w)]
    twin.project()
    np.testing.assert_allclose((sim.u, sim.v, sim.w),
                               (twin.u, twin.v, twin.w), atol=1e-9, rtol=0)
