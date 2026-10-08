"""Run the starts from the hug session.

    python run.py hug --n 32 --time 0.4
    python run.py sharper --n 64 --time 2
    python run.py tube --n 48 --time 2
    python run.py tubes --n 48 --time 2
    python run.py compare --n 48 --time 2
"""

import argparse
import json

from solver import (
    advance,
    build,
    energy,
    peak_speed,
    pressure_at_peak,
    pulled_tube,
    ring,
    slope_max,
    terms,
    timestep,
    two_tubes,
)


STARTS = {
    "hug": (ring, 0.01),
    "sharper": (lambda x, y, z: ring(x, y, z, sharp=4.0), 0.01),
    "tube": (pulled_tube, 0.002),
    "tubes": (two_tubes, 0.002),
}


def evolve(name, n, t_end, samples=12):
    psi, nu = STARTS[name]
    u, ops, report = build(n, psi=psi, nu=nu)
    dt = timestep(u, ops, nu)
    steps = max(1, int(round(t_end / dt)))
    dt = t_end / steps
    kx, ky, kz, k2, keep, dx = ops
    series = []
    every = max(1, steps // samples)
    for step in range(steps + 1):
        if step:
            u = advance(u, dt, ops, nu)
        if step % every == 0 or step == steps:
            steep, smooth = terms(u, ops, nu)
            series.append(
                {
                    "t": round(step * dt, 5),
                    "peak": peak_speed(u),
                    "energy": energy(u),
                    "slope": slope_max(u, kx, ky, kz, keep),
                    "steepening": steep,
                    "smoothing": smooth,
                    "peak_split": pressure_at_peak(u, ops, nu),
                }
            )
            print(
                f"t={series[-1]['t']:.2f} peak={series[-1]['peak']:.3f} "
                f"slope={series[-1]['slope']:.2f} energy={series[-1]['energy']:.3f} "
                f"steepening={steep:.6g} smoothing={smooth:.6g} "
                f"peak_carry={series[-1]['peak_split']['carry']:.6g} "
                f"peak_pressure={series[-1]['peak_split']['pressure']:.6g} "
                f"peak_smoothing={series[-1]['peak_split']['smoothing']:.6g}",
                flush=True,
            )
    return {"case": name, "n": n, "nu": nu, "dt": dt, "steps": steps,
            "start": report, "series": series}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=[*STARTS, "compare"])
    parser.add_argument("--n", type=int, default=32)
    parser.add_argument("--time", type=float, default=0.4)
    parser.add_argument("--out", default="")
    args = parser.parse_args()
    name = "tube" if args.case == "compare" else args.case
    out = evolve(name, args.n, args.time)
    if args.out:
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2)
        print("wrote", args.out)


if __name__ == "__main__":
    main()
