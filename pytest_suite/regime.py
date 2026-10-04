"""Regime-shift rules from Untitled2.py and Untitled3.py."""

import math


def deterministic_noise(seed):
    return 0.8 + ((math.sin(seed) + 1.0) / 5.0)


def regime_at(t, omega_crit=15.0, sigma=1.2):
    omega = 0.8 * t
    if omega < omega_crit:
        return omega, "Nuanced (Stealth)", 0.2 * omega
    return omega, "Aggressive (Nuance Lost)", (omega * (sigma ** 2)) / 2


def zone_status(zone, omega_crit=15.0):
    sigma = 1.0 + (zone * 0.4)
    system_friction = zone * 1.5
    noise = deterministic_noise(zone)
    omega = (10.0 + (zone * 2.2)) * noise
    if omega < omega_crit:
        x_optimal = 0.2 * omega / (1.0 + (system_friction * 0.1))
        status = "Nuanced (Stealth)"
    else:
        x_optimal = (omega * (sigma ** 2)) / 2
        breaking_point = 150.0 / (1.0 + (zone * 0.05))
        if x_optimal > breaking_point:
            status = "CRASHED"
        else:
            status = "Forced Aggression"
    return omega, x_optimal, status
