"""Geometry checks for the corrected field, including the missed axis case."""
import math
from dataclasses import replace

from initial_field import SmoothRingField
from stokes import meridional_velocity


def test_old_bounded_velocity_still_has_an_axis_cusp():
    z, R, alpha, h = 0.4, 1.5, 1.0, 1e-6
    center = meridional_velocity(0.0, z, R, alpha)[1]
    plus = (meridional_velocity(h, z, R, alpha)[1] - center) / h
    minus = (center - meridional_velocity(h, z, R, alpha)[1]) / h
    expected = 6 * R * z / alpha**2 * math.exp(-(R**2 + z**2) / alpha**2)
    assert abs(plus - expected) < 1e-5
    assert abs(minus + expected) < 1e-5
    assert plus - minus > 0.5


def test_smooth_axis_slopes_converge_to_same_value():
    f = SmoothRingField(radial_scale=math.sqrt(3.0))
    center = f.velocity(0.0, 0.0, 0.4)[2]
    slopes = []
    for h in (1e-3, 1e-4, 1e-5):
        plus = (f.velocity(h, 0.0, 0.4)[2] - center) / h
        minus = (center - f.velocity(-h, 0.0, 0.4)[2]) / h
        slopes.append(abs(plus - minus))
    assert slopes[2] < slopes[1] / 9 < slopes[0] / 81
    assert slopes[2] < 1e-5


def test_corrected_bump_is_symmetric_in_squared_radius():
    f = SmoothRingField()
    delta_q = 0.25
    left = math.sqrt(f.ring_radius**2 - delta_q)
    right = math.sqrt(f.ring_radius**2 + delta_q)
    assert abs(f.envelope(left, 0, 0) - f.envelope(right, 0, 0)) < 1e-12
    assert f.envelope(1.25, 0, 0) != f.envelope(1.75, 0, 0)


def test_corrected_streamline_tangency_midplane_and_reflection():
    f = SmoothRingField(radial_scale=math.sqrt(3.0))
    h = 1e-6
    for r in (0.2, 1.0, 1.5, 3.0):
        assert f.stream_function(r, 0.0) == 0.0
        assert f.velocity_cylindrical(r, 0.0)[2] == 0.0
        for z in (0.2, 0.8):
            ur, _, uz = f.velocity_cylindrical(r, z)
            lower = f.velocity_cylindrical(r, -z)
            assert ur == lower[0] and uz == -lower[2]
            dr = (f.stream_function(r+h, z) - f.stream_function(r-h, z)) / (2*h)
            dz = (f.stream_function(r, z+h) - f.stream_function(r, z-h)) / (2*h)
            assert abs(ur * dr + uz * dz) < 1e-8


def test_swirl_change_keeps_meridional_field_and_changes_azimuthal_speed():
    f = SmoothRingField()
    g = replace(f, swirl_rate=7.0)
    u, v = f.velocity_cylindrical(1.2, 0.4), g.velocity_cylindrical(1.2, 0.4)
    assert u[0] == v[0] and u[2] == v[2]
    assert abs(v[1] - 7 * u[1]) < 1e-14
    assert f.stream_function(1.2, 0.4) == g.stream_function(1.2, 0.4)
