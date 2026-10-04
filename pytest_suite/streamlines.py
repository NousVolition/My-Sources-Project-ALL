"""Meridional Stokes stream function (navier_stokes_framework.txt 3.B)."""

import math

R0 = 1.0
ALPHA = 1.0
R = 1.5
Z = 0.25


def phi(r, z, r0, alpha):
    return math.exp(-((r - r0) ** 2 + z ** 2) / alpha ** 2)


def psi(r, z, r0, alpha):
    return r ** 2 * z * phi(r, z, r0, alpha)


def meridional_velocity(r, z, r0, alpha):
    bump = phi(r, z, r0, alpha)
    u_r = -r * (1.0 - 2.0 * z ** 2 / alpha ** 2) * bump
    u_z = z * (2.0 - 2.0 * r * (r - r0) / alpha ** 2) * bump
    return u_r, u_z
