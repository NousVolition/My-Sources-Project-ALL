"""Stokes start on a box of length 18. Step with pressure and without."""

import math

from navier import PurePythonNavierStokes3D
from stokes import meridional_velocity, swirl

R0 = 1.5
ALPHA = 1.0
BOX = 18.0


def sample_stokes(sim, dx):
    mid = (sim.N - 1) / 2.0
    for i in range(sim.N):
        x = (i - mid) * dx
        for j in range(sim.N):
            y = (j - mid) * dx
            for k in range(sim.N):
                z = (k - mid) * dx
                r = math.hypot(x, y)
                rr = r if r >= 1e-8 else 1e-8
                u_r, u_z = meridional_velocity(rr, z, R0, ALPHA)
                u_th = swirl(rr, z, R0, ALPHA)
                cos_t, sin_t = (1.0, 0.0) if r < 1e-8 else (x / r, y / r)
                sim.u[i][j][k] = u_r * cos_t - u_th * sin_t
                sim.v[i][j][k] = u_r * sin_t + u_th * cos_t
                sim.w[i][j][k] = u_z
                sim.S[i][j][k] = 0.0


def max_div(sim):
    return max(
        abs(sim.divergence_at(i, j, k))
        for i in range(sim.N)
        for j in range(sim.N)
        for k in range(sim.N)
    )


def max_grad(sim):
    worst = 0.0
    for i in range(sim.N):
        for j in range(sim.N):
            for k in range(sim.N):
                ip, im, jp, jm, kp, km = sim.get_neighbors(i, j, k)
                comps = (
                    sim.u[ip][j][k] - sim.u[im][j][k],
                    sim.u[i][jp][k] - sim.u[i][jm][k],
                    sim.u[i][j][kp] - sim.u[i][j][km],
                    sim.v[ip][j][k] - sim.v[im][j][k],
                    sim.v[i][jp][k] - sim.v[i][jm][k],
                    sim.v[i][j][kp] - sim.v[i][j][km],
                    sim.w[ip][j][k] - sim.w[im][j][k],
                    sim.w[i][jp][k] - sim.w[i][jm][k],
                    sim.w[i][j][kp] - sim.w[i][j][km],
                )
                mag = math.sqrt(sum((c / (2.0 * sim.dx)) ** 2 for c in comps))
                worst = max(worst, mag)
    return worst


def main():
    for points in (16, 32):
        dx = BOX / points
        bare = PurePythonNavierStokes3D(N=points, dx=dx)
        held = PurePythonNavierStokes3D(N=points, dx=dx)
        sample_stokes(bare, dx)
        sample_stokes(held, dx)
        print(f"{points} points  start  div {max_div(bare):.5f}  grad {max_grad(bare):.5f}")
        for n in range(1, 5):
            bare.step(dt=0.02, P_U=0.0)
            held.step_with_pressure(dt=0.02, P_U=0.0)
            print(
                f"{points} points  step {n}"
                f"  no pressure  div {max_div(bare):.5f}  grad {max_grad(bare):.5f}"
                f"   |   with pressure  div {max_div(held):.5f}  grad {max_grad(held):.5f}"
            )


if __name__ == "__main__":
    main()
