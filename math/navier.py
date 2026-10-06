"""Periodic-grid fluid toy model from Code_.txt.

The original step omits pressure projection. The optional step_with_pressure
adds a compatible centered-difference projection on the periodic cube.
This remains an experimental discretization, with no convergence or stability
claim. Centered differences do not detect Nyquist checkerboard modes.
"""

import math


class PurePythonNavierStokes3D:
    def __init__(self, N=8, dx=1.0):
        self.N = N
        self.dx = dx
        self.nu = 0.01
        self.D = 0.02
        self.sigma = 0.5
        self.rho = 1.0
        self.gamma_crit = 1.5
        self.epsilon = 0.1
        self.e_k = [0.0, 0.0, 1.0]
        self.u = [[[0.05 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        self.v = [[[0.01 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        self.w = [[[0.02 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        self.S = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        self.rho_epsilon = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        mid = N / 2.0
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    r_sq = (i - mid) ** 2 + (j - mid) ** 2 + (k - mid) ** 2
                    self.rho_epsilon[i][j][k] = math.exp(-r_sq / (2.0 * (self.epsilon ** 2)))

    def get_neighbors(self, i, j, k):
        N = self.N
        return (i + 1) % N, (i - 1) % N, (j + 1) % N, (j - 1) % N, (k + 1) % N, (k - 1) % N

    def step(self, dt, P_U, backend="python"):
        if backend == "numpy":
            return self._step_numpy(dt, P_U)
        if backend != "python":
            raise ValueError("step backend must be python or numpy")
        N = self.N
        next_u = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        next_v = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        next_w = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        next_S = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    ip, im, jp, jm, kp, km = self.get_neighbors(i, j, k)
                    grad_Sx = (self.S[ip][j][k] - self.S[im][j][k]) / (2.0 * self.dx)
                    grad_Sy = (self.S[i][jp][k] - self.S[i][jm][k]) / (2.0 * self.dx)
                    grad_Sz = (self.S[i][j][kp] - self.S[i][j][km]) / (2.0 * self.dx)
                    grad_mag_S = math.sqrt(grad_Sx ** 2 + grad_Sy ** 2 + grad_Sz ** 2)
                    dux = (self.u[ip][j][k] - self.u[im][j][k]) / (2.0 * self.dx)
                    duy = (self.u[i][jp][k] - self.u[i][jm][k]) / (2.0 * self.dx)
                    duz = (self.u[i][j][kp] - self.u[i][j][km]) / (2.0 * self.dx)
                    dvx = (self.v[ip][j][k] - self.v[im][j][k]) / (2.0 * self.dx)
                    dvy = (self.v[i][jp][k] - self.v[i][jm][k]) / (2.0 * self.dx)
                    dvz = (self.v[i][j][kp] - self.v[i][j][km]) / (2.0 * self.dx)
                    dwx = (self.w[ip][j][k] - self.w[im][j][k]) / (2.0 * self.dx)
                    dwy = (self.w[i][jp][k] - self.w[i][jm][k]) / (2.0 * self.dx)
                    dwz = (self.w[i][j][kp] - self.w[i][j][km]) / (2.0 * self.dx)
                    adv_u = self.u[i][j][k] * dux + self.v[i][j][k] * duy + self.w[i][j][k] * duz
                    adv_v = self.u[i][j][k] * dvx + self.v[i][j][k] * dvy + self.w[i][j][k] * dvz
                    adv_w = self.u[i][j][k] * dwx + self.v[i][j][k] * dwy + self.w[i][j][k] * dwz
                    adv_S = self.u[i][j][k] * grad_Sx + self.v[i][j][k] * grad_Sy + self.w[i][j][k] * grad_Sz
                    delta_gamma = grad_mag_S - self.gamma_crit
                    switch_profile = 0.5 + 0.5 * math.tanh(delta_gamma / self.epsilon)
                    forcing_scalar = (self.sigma / self.rho) * delta_gamma * switch_profile
                    lap_u = (self.u[ip][j][k] + self.u[im][j][k] + self.u[i][jp][k] + self.u[i][jm][k] + self.u[i][j][kp] + self.u[i][j][km] - 6 * self.u[i][j][k]) / (self.dx ** 2)
                    lap_v = (self.v[ip][j][k] + self.v[im][j][k] + self.v[i][jp][k] + self.v[i][jm][k] + self.v[i][j][kp] + self.v[i][j][km] - 6 * self.v[i][j][k]) / (self.dx ** 2)
                    lap_w = (self.w[ip][j][k] + self.w[im][j][k] + self.w[i][jp][k] + self.w[i][jm][k] + self.w[i][j][kp] + self.w[i][j][km] - 6 * self.w[i][j][k]) / (self.dx ** 2)
                    lap_S = (self.S[ip][j][k] + self.S[im][j][k] + self.S[i][jp][k] + self.S[i][jm][k] + self.S[i][j][kp] + self.S[i][j][km] - 6 * self.S[i][j][k]) / (self.dx ** 2)
                    next_u[i][j][k] = self.u[i][j][k] + dt * (-adv_u + self.nu * lap_u - forcing_scalar * self.e_k[0])
                    next_v[i][j][k] = self.v[i][j][k] + dt * (-adv_v + self.nu * lap_v - forcing_scalar * self.e_k[1])
                    next_w[i][j][k] = self.w[i][j][k] + dt * (-adv_w + self.nu * lap_w - forcing_scalar * self.e_k[2])
                    next_S[i][j][k] = self.S[i][j][k] + dt * (-adv_S + self.D * lap_S + self.rho_epsilon[i][j][k] * P_U)
        self.u, self.v, self.w, self.S = next_u, next_v, next_w, next_S

    def _step_numpy(self, dt, P_U):
        """Array evaluation of the same explicit stencils, including S and force.

        The default loop implementation is retained as an independent reference.
        This backend keeps arrays after stepping to support larger-grid studies.
        """
        import numpy as np
        velocity = np.asarray((self.u, self.v, self.w))
        scalar = np.asarray(self.S)
        def derivative(a, d):
            return (np.roll(a, -1, d)-np.roll(a, 1, d))/(2*self.dx)
        def laplacian(a):
            return (sum(np.roll(a, -1, d)+np.roll(a, 1, d)
                        for d in range(3))-6*a)/self.dx**2
        grad_s = [derivative(scalar, d) for d in range(3)]
        delta = np.sqrt(sum(g*g for g in grad_s))-self.gamma_crit
        force = (self.sigma/self.rho)*delta*(.5+.5*np.tanh(delta/self.epsilon))
        updated = []
        for c in range(3):
            adv = sum(velocity[d]*derivative(velocity[c], d) for d in range(3))
            updated.append(velocity[c]+dt*(-adv+self.nu*laplacian(velocity[c])
                                           -force*self.e_k[c]))
        self.S = scalar+dt*(-sum(velocity[d]*grad_s[d] for d in range(3))
                           +self.D*laplacian(scalar)+np.asarray(self.rho_epsilon)*P_U)
        self.u, self.v, self.w = updated

    def _project_fft(self, atol, rtol):
        """The original centered projection, evaluated with array operations."""
        import numpy as np
        keep_arrays = isinstance(self.u, np.ndarray)
        velocity = np.asarray((self.u, self.v, self.w))
        def derivative(a, d):
            return (np.roll(a, -1, d)-np.roll(a, 1, d))/(2*self.dx)
        rhs = sum(derivative(velocity[d], d) for d in range(3))
        before = float(np.max(np.abs(rhs)))
        tolerance = atol+rtol*before
        wave = np.sin(2*np.pi*np.fft.fftfreq(self.N))/self.dx
        wave[0] = 0.
        if self.N % 2 == 0:
            wave[self.N//2] = 0.
        squared = (wave[:, None, None]**2+wave[None, :, None]**2
                   +wave[None, None, :]**2)
        potential_hat = np.zeros(rhs.shape, dtype=complex)
        np.divide(-np.fft.fftn(rhs), squared, out=potential_hat, where=squared > 0)
        chi = np.fft.ifftn(potential_hat).real
        corrected = [velocity[d]-derivative(chi, d) for d in range(3)]
        after = float(np.max(np.abs(sum(derivative(corrected[d], d) for d in range(3)))))
        if not np.isfinite(after) or after > tolerance:
            raise RuntimeError(f"Projected divergence {after:.3g} exceeds {tolerance:.3g}.")
        self.u, self.v, self.w = corrected if keep_arrays else [a.tolist() for a in corrected]
        self.projection_potential = chi if keep_arrays else chi.tolist()
        return dict(backend="fft", iterations=0, before=before, after=after,
                    tolerance=tolerance)

    def max_w(self):
        return max(max(max(row) for row in plane) for plane in self.w)

    def divergence_at(self, i, j, k):
        ip, im, jp, jm, kp, km = self.get_neighbors(i, j, k)
        return (
            self.u[ip][j][k] - self.u[im][j][k]
            + self.v[i][jp][k] - self.v[i][jm][k]
            + self.w[i][j][kp] - self.w[i][j][km]
        ) / (2.0 * self.dx)

    def project(self, iterations=2000, atol=1e-10, rtol=1e-10,
                backend="jacobi"):
        """Subtract G_h chi, where D_h G_h chi = D_h u.

        The default is damped Jacobi using only the standard library.
        backend='fft' uses optional NumPy and the SAME discrete operator.
        chi is a velocity correction potential; it is not pressure itself.
        For a time-step correction, rho*chi/dt is a pressure approximation.
        A failed Jacobi solve raises before changing the velocity.
        """
        if (not isinstance(iterations, int) or iterations < 1
                or not math.isfinite(self.dx) or self.dx <= 0
                or not math.isfinite(atol) or not math.isfinite(rtol)
                or atol <= 0 or rtol < 0):
            raise ValueError("Require positive dx, iterations, atol and nonnegative rtol.")
        if backend not in ("jacobi", "fft"):
            raise ValueError("backend must be 'jacobi' or 'fft'.")
        if backend == "fft":
            return self._project_fft(atol, rtol)
        N = self.N
        rhs = [[[self.divergence_at(i, j, k) for k in range(N)]
                for j in range(N)] for i in range(N)]
        before = max(abs(v) for plane in rhs for row in plane for v in row)
        tolerance = atol + rtol * before
        scale = 4.0 * self.dx**2
        used = 0
        chi = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        # Undamped Jacobi can alternate forever on quarter-wave modes.
        omega = 2.0 / 3.0
        plus2 = [(i + 2) % N for i in range(N)]
        minus2 = [(i - 2) % N for i in range(N)]
        residual = before
        for used in range(1, iterations + 1):
            if residual <= tolerance:
                used -= 1
                break
            nxt = [[[0.0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
            for i in range(N):
                for j in range(N):
                    for k in range(N):
                        neighbors = (
                            chi[plus2[i]][j][k] + chi[minus2[i]][j][k]
                            + chi[i][plus2[j]][k] + chi[i][minus2[j]][k]
                            + chi[i][j][plus2[k]] + chi[i][j][minus2[k]]
                        )
                        target = (neighbors - scale * rhs[i][j][k]) / 6.0
                        nxt[i][j][k] = (1.0 - omega) * chi[i][j][k] + omega * target
            chi = nxt
            if used % 10 == 0 or used == iterations:
                residual = 0.0
                for i in range(N):
                    for j in range(N):
                        for k in range(N):
                            lap = (
                                chi[plus2[i]][j][k] + chi[minus2[i]][j][k]
                                + chi[i][plus2[j]][k] + chi[i][minus2[j]][k]
                                + chi[i][j][plus2[k]] + chi[i][j][minus2[k]]
                                - 6.0 * chi[i][j][k]
                            ) / scale
                            residual = max(residual, abs(rhs[i][j][k] - lap))
        if residual > tolerance:
            raise RuntimeError(
                f"Projection did not converge: residual={residual:.3g}, "
                f"tolerance={tolerance:.3g}; increase iterations or use fft."
            )
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    ip, im, jp, jm, kp, km = self.get_neighbors(i, j, k)
                    self.u[i][j][k] -= (chi[ip][j][k] - chi[im][j][k]) / (2 * self.dx)
                    self.v[i][j][k] -= (chi[i][jp][k] - chi[i][jm][k]) / (2 * self.dx)
                    self.w[i][j][k] -= (chi[i][j][kp] - chi[i][j][km]) / (2 * self.dx)
        after = max(abs(self.divergence_at(i, j, k))
                    for i in range(N) for j in range(N) for k in range(N))
        self.projection_potential = chi
        if after > tolerance:
            raise RuntimeError(f"Projected divergence {after:.3g} exceeds {tolerance:.3g}.")
        return dict(backend=backend, iterations=used, before=before, after=after,
                    tolerance=tolerance)

    def step_with_pressure(self, dt, P_U, step_backend="python", **projection_options):
        """Apply the original explicit step, then a discrete projection."""
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be positive and finite.")
        if step_backend == "python":
            self.step(dt, P_U)
        else:
            self.step(dt, P_U, backend=step_backend)
        return self.project(**projection_options)

    def mean_S(self):
        total = sum(sum(sum(row) for row in plane) for plane in self.S)
        return total / (self.N ** 3)
