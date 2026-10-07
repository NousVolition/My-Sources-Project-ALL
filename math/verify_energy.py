"""Optional independent energy check; requires NumPy.

Run from the repository root with: python math/verify_energy.py
"""

import math

import numpy as np

from initial_field import SmoothRingField


def direct_energy(field, order=180):
    """Integrate the Cartesian speed squared using cylindrical volume measure."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    rmax = math.sqrt(field.ring_radius**2 + 8*field.radial_scale**2)
    zmax = 8*field.axial_scale
    radii = (nodes + 1)*rmax/2
    radial_weights = weights*rmax/2
    zs = nodes*zmax
    z_weights = weights*zmax
    total = 0.0
    for r, wr in zip(radii, radial_weights):
        for z, wz in zip(zs, z_weights):
            components = field.velocity(float(r), 0.0, float(z))
            total += wr*wz*r*sum(u*u for u in components)
    return math.pi*total


def main():
    cases = (
        SmoothRingField(),
        SmoothRingField(ring_radius=0.0),
        SmoothRingField(0.8, 1.2, 0.7, 0.4, -0.9),
        SmoothRingField(3.0, 1.0, 2.0, 0.25, 0.5),
    )
    for field in cases:
        exact = field.kinetic_energy()
        numerical = direct_energy(field)
        relative_error = abs(exact-numerical)/exact
        if relative_error >= 1e-9:
            raise AssertionError((field, relative_error))
        print(f"{field}\n  exact E={exact:.15g}; quadrature E={numerical:.15g}; "
              f"relative error={relative_error:.3g}")


if __name__ == "__main__":
    main()
