"""Small, inspectable solvers for the user-supplied dynamics examples.

All computations are mathematical model tests, not physical observations.
Times and states are dimensionless unless an example explicitly supplies units.
"""
import numpy as np


def advance(fun, y, h, method="rk4"):
    k1 = fun(y)
    if method == "euler":
        return y + h * k1
    if method == "heun":
        return y + h * (k1 + fun(y + h * k1)) / 2
    if method != "rk4":
        raise ValueError(method)
    k2 = fun(y + h * k1 / 2)
    k3 = fun(y + h * k2 / 2)
    k4 = fun(y + h * k3)
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def integrate(fun, initial, duration, h, method="rk4"):
    count = round(duration / h)
    if not np.isclose(count * h, duration, atol=1e-12, rtol=0):
        raise ValueError("duration must be an integer multiple of the step")
    t = np.arange(count + 1) * h
    state = np.empty((count + 1, np.atleast_1d(initial).size))
    state[0] = initial
    for i in range(count):
        state[i + 1] = advance(fun, state[i], h, method)
    return t, state


def logistic_exact(t, initial, rate=1., capacity=1.):
    return capacity / (1 + (capacity / initial - 1) * np.exp(-rate * t))


def lorenz(y, rho=28., sigma=10., beta=8/3):
    x, v, z = y
    return np.array([sigma * (v - x), x * (rho - z) - v, x * v - beta * z])


def jacobian(y, rho=28., sigma=10., beta=8/3):
    x, v, z = y
    return np.array([[-sigma, sigma, 0], [rho-z, -1, -x], [v, x, -beta]])


def variational(y, rho=28., sigma=10., beta=8/3):
    return np.r_[lorenz(y[:3], rho, sigma, beta),
                 (jacobian(y[:3], rho, sigma, beta) @ y[3:].reshape(3, 3)).ravel()]


def lyapunov(initial, h=.005, transient=30., duration=180., block=10.):
    """Benettin QR method; retain block rates, not a fit to one separation curve."""
    state = np.array(initial, dtype=float)
    for _ in range(round(transient/h)):
        state = advance(lorenz, state, h)
    y = np.r_[state, np.eye(3).ravel()]
    qr_interval = .1
    stride = round(qr_interval / h)
    block_stride = round(block / qr_interval)
    block_rates, states, running = [], [], np.zeros(3)
    totals = np.zeros(3)
    for i in range(round(duration / qr_interval)):
        for _ in range(stride):
            y = advance(variational, y, h)
        q, r = np.linalg.qr(y[3:].reshape(3, 3))
        logs = np.log(np.abs(np.diag(r)))
        running += logs
        totals += logs
        y[3:] = q.ravel()
        states.append(y[:3].copy())
        if (i + 1) % block_stride == 0:
            block_rates.append(running / block)
            running[:] = 0
    return totals / duration, np.asarray(block_rates), np.asarray(states)


def delayed_branch(t, delay=0., sign=1.):
    age = np.maximum(np.asarray(t) - delay, 0)
    x = sign * (2 * age / 3)**1.5
    derivative = sign * np.sqrt(2 * age / 3)
    return x, derivative


def wheel_rhs(n, q0=60., q1=28., leak=1., drag=10., torque=10.):
    """Periodic spectral mass transport plus angular-momentum ODE.

    m_t=Q-leak*m-omega*m_theta; omega_t=-drag*omega+torque*a1.
    torque includes the geometric pi*g*radius/inertia factor.
    """
    theta = 2 * np.pi * np.arange(n) / n
    wave = np.fft.fftfreq(n, d=1/n)
    force = q0 + q1 * np.cos(theta)
    def rhs(y):
        m, omega = y[:-1], y[-1]
        derivative = np.fft.ifft(1j * wave * np.fft.fft(m)).real
        a1 = 2 * np.mean(m * np.sin(theta))
        return np.r_[force - leak*m - omega*derivative, -drag*omega + torque*a1]
    return theta, rhs


def wheel_project(states, theta):
    mass_density = states[:, :-1]
    a = 2 * np.mean(mass_density * np.sin(theta), axis=1)
    b = 2 * np.mean(mass_density * np.cos(theta), axis=1)
    # For K=1, drag=torque=10, q1=28: x=omega, y=a, z=28-b.
    return np.c_[states[:, -1], a, 28-b]
