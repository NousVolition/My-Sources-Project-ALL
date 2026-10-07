"""Smooth, axisymmetric initial data on R^3.

This module constructs an initial field and its exact kinetic energy per unit
density. It does not evolve the Navier-Stokes equations.
"""

from dataclasses import dataclass, replace
import math


@dataclass(frozen=True)
class SmoothRingField:
    """Parameters R, a, b, A, B in the accompanying mathematical note.

    ring_radius, radial_scale, axial_scale have units of length.
    meridional_rate and swirl_rate have units of inverse time.
    radial_scale is not the local Gaussian width in r: near R > 0, that width
    is approximately radial_scale**2 / (2 * ring_radius).
    """

    ring_radius: float = 1.5
    radial_scale: float = 1.0
    axial_scale: float = 1.0
    meridional_rate: float = 1.0
    swirl_rate: float = 1.0

    def __post_init__(self):
        values = (
            self.ring_radius, self.radial_scale, self.axial_scale,
            self.meridional_rate, self.swirl_rate,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("All parameters must be finite.")
        if self.ring_radius < 0:
            raise ValueError("ring_radius must be nonnegative.")
        if self.radial_scale <= 0 or self.axial_scale <= 0:
            raise ValueError("Both scales must be positive.")

    def envelope(self, x, y, z):
        q = x * x + y * y
        return math.exp(
            -((q - self.ring_radius**2) ** 2) / self.radial_scale**4
            - z * z / self.axial_scale**2
        )

    def stream_function(self, r, z):
        if r < 0:
            raise ValueError("Cylindrical radius must be nonnegative.")
        return self.meridional_rate * r * r * z * self.envelope(r, 0.0, z)

    def velocity(self, x, y, z):
        """Return Cartesian (u_x, u_y, u_z), including on the symmetry axis."""
        q = x * x + y * y
        phi = self.envelope(x, y, z)
        h = 1.0 - 2.0 * z * z / self.axial_scale**2
        c = 2.0 - 4.0 * q * (q - self.ring_radius**2) / self.radial_scale**4
        a, b = self.meridional_rate, self.swirl_rate
        return (
            (-a * x * h - b * y) * phi,
            (-a * y * h + b * x) * phi,
            a * z * c * phi,
        )

    def velocity_cylindrical(self, r, z):
        """Return (u_r, u_theta, u_z); r=0 uses the smooth extension."""
        if r < 0:
            raise ValueError("Cylindrical radius must be nonnegative.")
        return self.velocity(r, 0.0, z)

    def kinetic_energy(self):
        """Exact E = (1/2) integral over R^3 of |u_0|^2, per unit density."""
        radius2 = self.ring_radius**2
        radial2 = self.radial_scale**2
        radial4 = radial2**2
        axial2 = self.axial_scale**2
        a2 = self.meridional_rate**2
        b2 = self.swirl_rate**2
        z0 = self.axial_scale * math.sqrt(math.pi / 2.0)
        m0 = (
            radial2 * math.sqrt(math.pi) / (2.0 * math.sqrt(2.0))
            * (1.0 + math.erf(math.sqrt(2.0) * radius2 / radial2))
        )
        m1 = radius2 * m0 + radial4 / 4.0 * math.exp(
            -2.0 * radius2**2 / radial4
        )
        return math.pi * z0 / 2.0 * (
            (b2 + 3.0 * a2 / 4.0 + a2 * axial2 * radius2 / radial4) * m1
            + 3.0 * a2 * axial2 / 4.0 * m0
        )

    def with_energy(self, target):
        """Scale both rates together, retaining the ratio of swirl to flow."""
        if not math.isfinite(target) or target < 0:
            raise ValueError("Target energy must be finite and nonnegative.")
        if target == 0:
            return replace(self, meridional_rate=0.0, swirl_rate=0.0)
        current = self.kinetic_energy()
        if current == 0:
            raise ValueError("A zero field cannot be scaled to positive energy.")
        factor = math.sqrt(target / current)
        return replace(
            self,
            meridional_rate=self.meridional_rate * factor,
            swirl_rate=self.swirl_rate * factor,
        )
