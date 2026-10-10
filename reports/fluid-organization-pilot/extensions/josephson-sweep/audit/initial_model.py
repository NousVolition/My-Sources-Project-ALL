"""Dimensionless, noiseless, overdamped RSJ pair with a resistive load.

M phi' = i - sin(phi), M = I + alpha * ones((2,2)).
tau = (2e Ic r / hbar) t; V_total/(Ic r) = phi1' + phi2'.
This implements the same circuit as fluid_pilot/dynamics.py:circuit_rhs.
All phases remain unwrapped. Integrals use the RK4 stage quadrature.
"""
import numpy as np


def rhs(phi, bias, alpha):
    phi = np.asarray(phi, dtype=float)
    b = np.asarray(bias)[..., None] - np.sin(phi)
    a = np.asarray(alpha)[..., None]
    return b - a / (1 + 2 * a) * b.sum(axis=-1, keepdims=True)


def potential(phi):
    return -np.cos(phi).sum(axis=-1)


def dissipation(v, alpha):
    return (v * v).sum(axis=-1) + np.asarray(alpha) * v.sum(axis=-1)**2


def step(phi, bias, alpha, dt):
    k1 = rhs(phi, bias, alpha)
    k2 = rhs(phi + dt/2*k1, bias, alpha)
    k3 = rhs(phi + dt/2*k2, bias, alpha)
    k4 = rhs(phi + dt*k3, bias, alpha)
    p = phi + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    d = dt/6*(dissipation(k1, alpha) + 2*dissipation(k2, alpha)
              + 2*dissipation(k3, alpha) + dissipation(k4, alpha))
    return p, d


def integrate(phi0, bias, alpha, duration, dt, stride=1.0):
    """Batch independent trajectories; return phase and dissipated-energy samples."""
    if dt <= 0 or duration <= 0 or stride <= 0 or np.any(np.asarray(alpha) < 0):
        raise ValueError('Require positive times and nonnegative resistance ratio')
    n, every = round(duration/dt), round(stride/dt)
    if every < 1 or not np.isclose(n*dt, duration) or not np.isclose(every*dt, stride) or n % every:
        raise ValueError('duration and stride must be integral multiples of dt, and duration of stride')
    p = np.asarray(phi0, dtype=float).copy()
    if p.shape[-1] != 2 or not np.isfinite(p).all():
        raise ValueError('Expected finite phase pairs')
    acc = np.zeros(p.shape[:-1]); samples = [p.copy()]; energies = [acc.copy()]
    for j in range(n):
        p, d = step(p, bias, alpha, dt)
        acc += d
        if (j + 1) % every == 0:
            samples.append(p.copy()); energies.append(acc.copy())
    phases, dissipated = np.asarray(samples), np.asarray(energies)
    if not np.isfinite(phases).all():
        raise FloatingPointError('Nonfinite integration')
    return np.arange(len(samples))*stride, phases, dissipated


def voltage(start, end, window):
    if window <= 0:
        raise ValueError('Voltage window must be positive')
    return (np.asarray(end) - np.asarray(start)) / window


def power_residual(phases, dissipated, bias):
    work = np.asarray(bias)*(phases[-1] - phases[0]).sum(axis=-1)
    du = potential(phases[-1]) - potential(phases[0])
    residual = work - du - dissipated[-1]
    scale = 1 + abs(work) + abs(du) + dissipated[-1]
    return residual / scale


def synchronized_voltage(bias, alpha):
    """Infinite-time total voltage on the invariant synchronized branch only."""
    b = np.asarray(bias)
    return 2*np.sign(b)*np.sqrt(np.maximum(b*b - 1, 0))/(1 + 2*np.asarray(alpha))
