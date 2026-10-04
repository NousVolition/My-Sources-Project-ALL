"""Analytic axisymmetric initial velocity field from navier_stokes_framework.txt (Section 3B).

The meridional flow is derived from the localized Stokes stream function

    psi(r, z) = r^2 * z * phi(r, z),
    phi(r, z) = exp(-((r - r0)^2 + z^2) / alpha^2),

via the Stokes relations u_r = -(1/r) dpsi/dz and u_z = (1/r) dpsi/dr,
plus an independent swirl component u_theta = r * phi.
"""

import math

# Grid step used for the central-difference recovery in velocity_from_psi.
H = 1e-6


def phi(r, z, r0, alpha):
    """Gaussian bump centered at (r0, 0) with width alpha."""
    return math.exp(-(((r - r0) ** 2) + z ** 2) / (alpha ** 2))


def psi(r, z, r0, alpha):
    """Stokes stream function; zero on the midplane and odd in z."""
    return (r ** 2) * z * phi(r, z, r0, alpha)


def meridional_velocity(r, z, r0, alpha):
    """Written (u_r, u_z) components of the meridional velocity field."""
    bump = phi(r, z, r0, alpha)
    u_r = -r * (1.0 - 2.0 * (z ** 2) / (alpha ** 2)) * bump
    u_z = z * (2.0 - 2.0 * r * (r - r0) / (alpha ** 2)) * bump
    return u_r, u_z


def swirl(r, z, r0, alpha):
    """Azimuthal component u_theta, independent of the meridional flow."""
    return r * phi(r, z, r0, alpha)


def velocity_from_psi(r, z, r0, alpha):
    """Central-difference recovery of (u_r, u_z) from the stream function."""
    dpsi_dz = (psi(r, z + H, r0, alpha) - psi(r, z - H, r0, alpha)) / (2.0 * H)
    dpsi_dr = (psi(r + H, z, r0, alpha) - psi(r - H, z, r0, alpha)) / (2.0 * H)
    return -dpsi_dz / r, dpsi_dr / r


def cylindrical_divergence(r, z, r0, alpha):
    """Divergence of the meridional field in cylindrical coordinates.

    Returns (1/r) d(r u_r)/dr + du_z/dz with the derivatives evaluated
    analytically from the closed-form components.
    """
    bump = phi(r, z, r0, alpha)
    a2 = alpha ** 2
    du_z_dz = (2.0 - 2.0 * r * (r - r0) / a2) * (1.0 - 2.0 * (z ** 2) / a2) * bump
    inv_r_dr_u_r_dr = (
        -2.0
        * (1.0 - 2.0 * (z ** 2) / a2)
        * (1.0 - r * (r - r0) / a2)
        * bump
    )
    return inv_r_dr_u_r_dr + du_z_dz
