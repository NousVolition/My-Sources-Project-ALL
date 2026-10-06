"""Catch grid alignment errors that can impersonate lack of convergence."""
import numpy as np

from hug_refinement import compare_fields, sample_on_grid


def wave(n):
    coords = -np.pi+(np.arange(n)+.5)*2*np.pi/n
    x, y, z = np.meshgrid(coords, coords, coords, indexing="ij")
    return np.sin(x)+.3*np.cos(2*y)-.2*np.sin(z)


def test_transfer_matches_physical_centers_with_fourth_order_error():
    target = wave(16)
    errors = [np.max(np.abs(sample_on_grid(wave(n), 16)-target)) for n in (32, 64)]
    assert errors[1] < errors[0]/12
    assert errors[1] < 2e-6


def test_identical_fields_have_zero_comparison_error():
    a = np.array([wave(12), wave(12)*2, -wave(12)])
    result = compare_fields(a, a)
    assert result == dict(velocity_relative_l2=0., gradient_relative_l2=0.)


def test_amplitude_error_is_detected_in_velocity_and_gradient():
    a = np.array([wave(12), wave(12)*2, -wave(12)])
    result = compare_fields(a*.9, a)
    assert abs(result['velocity_relative_l2']-.1) < 1e-14
    assert abs(result['gradient_relative_l2']-.1) < 1e-14
