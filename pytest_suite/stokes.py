"""Stokes stream function start profile from navier_stokes_framework.txt."""

import math


def meridional_velocity(r, z, r0, alpha):
    """Meridional (u_r, u_z) from psi = r^2 z exp(-((r - r0)^2 + z^2) / alpha^2)."""
    envelope = math.exp(-((r - r0) ** 2 + z ** 2) / (alpha ** 2))
    u_r = -r * (1.0 - 2.0 * z ** 2 / (alpha ** 2)) * envelope
    u_z = z * (2.0 - 2.0 * r * (r - r0) / (alpha ** 2)) * envelope
    return u_r, u_z


def swirl(r, z, r0, alpha):
    """Swirling component u_theta = r exp(-((r - r0)^2 + z^2) / alpha^2)."""
    return r * math.exp(-((r - r0) ** 2 + z ** 2) / (alpha ** 2))
