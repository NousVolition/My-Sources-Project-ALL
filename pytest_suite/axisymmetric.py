"""Initial axisymmetric field from the Stokes stream-function note.

psi builds the starting meridional velocity. It is not the Stokes equations
and it does not step Navier-Stokes.
"""

import math


def phi(r, z, r0, alpha):
    return math.exp(-(((r - r0) ** 2) + (z ** 2)) / (alpha ** 2))


def psi(r, z, r0, alpha):
    """psi_0 = r^2 * z * Gaussian bump centered at r = r0."""
    return (r ** 2) * z * phi(r, z, r0, alpha)


def meridional_velocity(r, z, r0, alpha):
    bump = phi(r, z, r0, alpha)
    u_r = -r * (1.0 - 2.0 * (z ** 2) / (alpha ** 2)) * bump
    u_z = z * (2.0 - 2.0 * r * (r - r0) / (alpha ** 2)) * bump
    return u_r, u_z


def swirl(r, z, r0, alpha):
    """u_theta = r * phi, added by hand. Continuity does not fix it."""
    return r * phi(r, z, r0, alpha)


def velocity_from_psi(r, z, r0, alpha, h=1e-6):
    """Central difference of the stream-function definitions."""
    dpsi_dz = (psi(r, z + h, r0, alpha) - psi(r, z - h, r0, alpha)) / (2.0 * h)
    dpsi_dr = (psi(r + h, z, r0, alpha) - psi(r - h, z, r0, alpha)) / (2.0 * h)
    u_r = -(1.0 / r) * dpsi_dz
    u_z = (1.0 / r) * dpsi_dr
    return u_r, u_z


def cylindrical_divergence(r, z, r0, alpha, h=1e-6):
    """(1/r) d(r u_r)/dr + du_z/dz for the analytic meridional field."""
    def ur(rr, zz):
        return meridional_velocity(rr, zz, r0, alpha)[0]

    def uz(rr, zz):
        return meridional_velocity(rr, zz, r0, alpha)[1]

    d_r_ur = ((r + h) * ur(r + h, z) - (r - h) * ur(r - h, z)) / (2.0 * h)
    d_uz = (uz(r, z + h) - uz(r, z - h)) / (2.0 * h)
    return d_r_ur / r + d_uz
