"""Localized solenoidal Stokes start from navier_stokes_framework.txt section B.

Meridian velocities derive from the Stokes stream function
psi0(r, z) = r^2 * z * exp(-[(r - r_0)^2 + z^2] / alpha^2)
via u_r = -(1/r) d(psi)/dz and u_z = (1/r) d(psi)/dr, plus an independent
swirl component u_theta.
"""

import math


def _envelope(r, z, r_0, alpha):
    return math.exp(-(((r - r_0) ** 2 + z ** 2) / (alpha ** 2)))


def meridional_velocity(r, z, r_0, alpha):
    e = _envelope(r, z, r_0, alpha)
    u_r = -r * (1.0 - 2.0 * z ** 2 / alpha ** 2) * e
    u_z = z * (2.0 - 2.0 * r * (r - r_0) / alpha ** 2) * e
    return u_r, u_z


def swirl(r, z, r_0, alpha):
    return r * _envelope(r, z, r_0, alpha)
