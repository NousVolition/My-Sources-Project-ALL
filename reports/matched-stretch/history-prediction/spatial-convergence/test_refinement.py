"""Known Fourier fields validate state transfer across nonnested grids."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_refinement import transfer, Flow


def field(n):
    q = np.arange(n)*2*np.pi/n
    x, y, z = np.meshgrid(q, q, q, indexing='ij')
    return np.array([.7+np.cos(3*y-2*z), -.2+np.sin(2*x+z), .1+np.cos(x-3*y)])


def test_transfer_matches_known_non_nested_fourier_field():
    small = Flow(14, .02, workers=1)
    large = Flow(22, .02, workers=1)
    h = small.hat(field(14))
    actual = large.real(transfer(h, 22))
    np.testing.assert_allclose(actual, field(22), atol=8e-15)
    np.testing.assert_allclose(transfer(transfer(h, 22), 14), h, atol=2e-12)


def test_transfer_preserves_energy_mean_and_incompressibility():
    small = Flow(14, .02, workers=1)
    large = Flow(22, .02, workers=1)
    h = small.project(small.hat(field(14)))
    lifted = transfer(h, 22)
    np.testing.assert_allclose(small.inner(h, h), large.inner(lifted, lifted), rtol=4e-15)
    np.testing.assert_allclose(lifted[:, 0, 0, 0]/22**3, [.7, -.2, .1], atol=2e-15)
    divergence = large.real(1j*sum(k*c for k, c in zip(large.k, lifted)))
    assert abs(divergence).max() < 2e-14


def test_restriction_discards_only_unresolvable_modes():
    large = Flow(22, .02, workers=1)
    small = Flow(14, .02, workers=1)
    q = np.arange(22)*2*np.pi/22
    high = np.broadcast_to(np.cos(9*q)[None, :, None], (22, 22, 22))
    u = field(22); u[0] += high
    restricted = small.real(transfer(large.hat(u), 14))
    np.testing.assert_allclose(restricted, field(14), atol=1e-14)
