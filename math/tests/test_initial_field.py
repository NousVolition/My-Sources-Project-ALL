import math
import unittest

from initial_field import SmoothRingField


class InitialFieldTests(unittest.TestCase):
    def test_axis_values_use_smooth_cartesian_extension(self):
        f = SmoothRingField()
        z = 0.4
        ux, uy, uz = f.velocity(0.0, 0.0, z)
        self.assertEqual((ux, uy), (0.0, 0.0))
        expected = 2.0 * z * math.exp(-1.5**4 - z**2)
        self.assertAlmostEqual(uz, expected, places=14)

    def test_stream_function_derivatives_match_velocity(self):
        f = SmoothRingField(axial_scale=0.8, meridional_rate=0.7)
        r, z, h = 1.1, 0.3, 1e-5
        ur, _, uz = f.velocity_cylindrical(r, z)
        ur_fd = -(f.stream_function(r, z+h) - f.stream_function(r, z-h)) / (2*h*r)
        uz_fd = (f.stream_function(r+h, z) - f.stream_function(r-h, z)) / (2*h*r)
        self.assertAlmostEqual(ur, ur_fd, places=7)
        self.assertAlmostEqual(uz, uz_fd, places=7)

    def test_cartesian_divergence_including_axis(self):
        f = SmoothRingField()
        h = 1e-5
        for point in ((0.0, 0.0, 0.4), (1.2, 0.3, 0.4), (2.0, -0.8, -0.7)):
            divergence = 0.0
            for axis in range(3):
                plus, minus = list(point), list(point)
                plus[axis] += h
                minus[axis] -= h
                divergence += (f.velocity(*plus)[axis] - f.velocity(*minus)[axis]) / (2*h)
            self.assertLess(abs(divergence), 1e-7)

    def test_rotation_about_axis_rotates_velocity(self):
        f = SmoothRingField(swirl_rate=-0.5)
        x, y, z = 0.8, 1.1, 0.4
        ux, uy, uz = f.velocity(x, y, z)
        rotated = f.velocity(-y, x, z)
        for actual, expected in zip(rotated, (-uy, ux, uz)):
            self.assertAlmostEqual(actual, expected, places=13)

    def test_energy_scaling_and_normalization(self):
        f = SmoothRingField()
        doubled = SmoothRingField(meridional_rate=2.0, swirl_rate=2.0)
        self.assertAlmostEqual(doubled.kinetic_energy(), 4.0*f.kinetic_energy(), places=12)
        normalized = f.with_energy(1.0)
        self.assertAlmostEqual(normalized.kinetic_energy(), 1.0, places=13)

    def test_zero_field_energy(self):
        f = SmoothRingField(meridional_rate=0.0, swirl_rate=0.0)
        self.assertEqual(f.kinetic_energy(), 0.0)
        with self.assertRaises(ValueError):
            f.with_energy(1.0)

    def test_invalid_geometry_is_rejected(self):
        for kwargs in ({"ring_radius": -1.0}, {"radial_scale": 0.0}, {"axial_scale": -1.0}):
            with self.assertRaises(ValueError):
                SmoothRingField(**kwargs)


if __name__ == "__main__":
    unittest.main()
