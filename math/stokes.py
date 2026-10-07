"""Original radial Gaussian, retained for comparison with the submitted tests.

This constructs an initial field. It neither solves Stokes equations nor
evolves Navier-Stokes. For r0 > 0 its Cartesian velocity is generally NOT C1
at the axis, even though the radial component is bounded there.
Use initial_field.SmoothRingField for the smooth construction.
"""

import math


def phi(r, z, r0, alpha):
    return math.exp(-(((r - r0) ** 2) + (z ** 2)) / (alpha ** 2))


def psi(r, z, r0, alpha):
    return (r ** 2) * z * phi(r, z, r0, alpha)


def meridional_velocity(r, z, r0, alpha):
    bump = phi(r, z, r0, alpha)
    u_r = -r * (1.0 - 2.0 * (z ** 2) / (alpha ** 2)) * bump
    u_z = z * (2.0 - 2.0 * r * (r - r0) / (alpha ** 2)) * bump
    return u_r, u_z


def swirl(r, z, r0, alpha):
    return r * phi(r, z, r0, alpha)


def velocity_from_psi(r, z, r0, alpha, h=1e-6):
    """Finite-difference diagnostic for r > h > 0 only."""
    if not r > h > 0:
        raise ValueError("Require r > h > 0 for this cylindrical diagnostic.")
    dpsi_dz = (psi(r, z + h, r0, alpha) - psi(r, z - h, r0, alpha)) / (2.0 * h)
    dpsi_dr = (psi(r + h, z, r0, alpha) - psi(r - h, z, r0, alpha)) / (2.0 * h)
    return -dpsi_dz / r, dpsi_dr / r


def cylindrical_divergence(r, z, r0, alpha, h=1e-6):
    """Finite-difference diagnostic away from the axis."""
    if not r > h > 0:
        raise ValueError("Require r > h > 0 for this cylindrical diagnostic.")
    def ur(rr, zz):
        return meridional_velocity(rr, zz, r0, alpha)[0]
    def uz(rr, zz):
        return meridional_velocity(rr, zz, r0, alpha)[1]
    d_r_ur = ((r + h) * ur(r + h, z) - (r - h) * ur(r - h, z)) / (2.0 * h)
    d_uz = (uz(r, z + h) - uz(r, z - h)) / (2.0 * h)
    return d_r_ur / r + d_uz
