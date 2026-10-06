"""Unforced, periodic Navier--Stokes follow-through for the A/B/C example.

Domain: [0, 2*pi)^3. Fourier differentiation, strict 2/3 truncation,
Leray pressure projection, and explicit RK4. This is a separate experiment
from the legacy finite-difference cube, whose scalar forcing is not used.
The names label constructed velocity fields, not people or risk classes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np


class PeriodicFlow:
    def __init__(self, points=32, viscosity=0.05):
        if not isinstance(points, int) or points < 12 or points % 2:
            raise ValueError("points must be an even integer >= 12")
        if not math.isfinite(viscosity) or viscosity <= 0:
            raise ValueError("viscosity must be positive and finite")
        self.n, self.nu = points, viscosity
        k = np.rint(np.fft.fftfreq(points) * points)
        kz = np.rint(np.fft.rfftfreq(points) * points)
        self.k = np.array(np.broadcast_arrays(
            k[:, None, None], k[None, :, None], kz[None, None, :]))
        self.k2 = np.sum(self.k**2, axis=0)
        self.inv_k2 = np.divide(1., self.k2, out=np.zeros_like(self.k2),
                                 where=self.k2 > 0)
        # Strict inequality excludes the aliasing boundary when N/3 is integer.
        self.keep = np.max(np.abs(self.k), axis=0) < points / 3
        self.cutoff = int(np.max(np.abs(self.k[:, self.keep])))
        self.low = self.k2 <= 1
        self.edge = self.keep & (np.max(np.abs(self.k), axis=0) >= .8*self.cutoff)
        self.weights = np.full(points//2+1, 2.)
        self.weights[[0, -1]] = 1.
        q = np.arange(points) * 2*np.pi / points
        self.xyz = np.meshgrid(q, q, q, indexing="ij", sparse=True)

    def fft(self, v):
        return np.fft.rfftn(v, axes=(-3, -2, -1), norm="forward")

    def real(self, h):
        return np.fft.irfftn(h, s=(self.n,)*3, axes=(-3, -2, -1), norm="forward")

    def norm2(self, h):
        """Spatial mean of squared magnitude, by Parseval."""
        return float(np.sum(np.abs(h)**2 * self.weights))

    def project(self, h):
        dot = np.sum(self.k*h, axis=0)
        return (h-self.k*(dot*self.inv_k2)) * self.keep

    def initial(self, field, epsilon=0.5):
        x, y, z = self.xyz
        v = np.zeros((3, self.n, self.n, self.n))
        if field in ("A", "B"):
            v[1] = math.sqrt(6)*np.cos(2*x) * (1 if field == "A" else -1)
        elif field in ("C", "C3D"):
            v[0], v[1], v[2] = np.cos(2*x+y), -2*np.cos(2*x+y), np.cos(2*x)
            if field == "C3D":
                if not math.isfinite(epsilon) or epsilon < 0:
                    raise ValueError("epsilon must be nonnegative and finite")
                v[0] += epsilon*np.cos(2*z)
                v[1] += epsilon*np.cos(2*x)
                v[2] += epsilon*np.cos(2*y)
                # Orthogonal Fourier modes give <|C+epsilon D|^2>=3+1.5 eps^2.
                v *= math.sqrt(3/(3+1.5*epsilon**2))
        else:
            raise ValueError("field must be A, B, C, or C3D")
        return self.project(self.fft(v))

    def curl(self, h):
        return 1j*np.stack((self.k[1]*h[2]-self.k[2]*h[1],
                            self.k[2]*h[0]-self.k[0]*h[2],
                            self.k[0]*h[1]-self.k[1]*h[0]))

    def rhs(self, h):
        v, w = self.real(h), self.real(self.curl(h))
        cross = np.stack((v[1]*w[2]-v[2]*w[1],
                          v[2]*w[0]-v[0]*w[2],
                          v[0]*w[1]-v[1]*w[0]))
        # P(u x curl u) = -P((u . grad)u); P accounts for pressure.
        return self.project(self.fft(cross))-self.nu*self.k2*h

    def dissipation(self, h):
        return self.nu*float(np.sum(self.k2*np.abs(h)**2*self.weights))

    def step(self, h, dt):
        k1 = self.rhs(h)
        b = h+dt*k1/2
        k2 = self.rhs(b)
        c = h+dt*k2/2
        k3 = self.rhs(c)
        d = h+dt*k3
        k4 = self.rhs(d)
        dissipated = dt/6*(self.dissipation(h)+2*self.dissipation(b)
                           +2*self.dissipation(c)+self.dissipation(d))
        return self.project(h+dt/6*(k1+2*k2+2*k3+k4)), dissipated

    def measure(self, h, t, dt, initial_energy, dissipated):
        v = self.real(h)
        grad = self.real(1j*self.k[None, ...]*h[:, None, ...])
        w_hat = self.curl(h)
        w = self.real(w_hat)
        energy = self.norm2(h)/2
        grad_squared = float(np.sum(self.k2*np.abs(h)**2*self.weights))
        stretching = float(np.mean(np.einsum("iabc,ijabc,jabc->abc", w, grad, w)))
        row = dict(time=float(t), energy=energy,
                   low_mode_energy=self.norm2(h*self.low)/2,
                   enstrophy=self.norm2(w_hat)/2,
                   max_speed=float(np.max(np.sqrt(np.sum(v*v, axis=0)))),
                   max_gradient=float(np.max(np.sqrt(np.sum(grad*grad, axis=(0, 1))))),
                   max_vorticity=float(np.max(np.sqrt(np.sum(w*w, axis=0)))),
                   rms_divergence=math.sqrt(self.norm2(1j*np.sum(self.k*h, axis=0))),
                   gradient_squared_mean=grad_squared,
                   stretching_mean=stretching,
                   enstrophy_dissipation=self.nu*float(np.sum(self.k2*np.abs(w_hat)**2*self.weights)),
                   z_derivative_rms=math.sqrt(self.norm2(1j*self.k[2]*h)),
                   edge_energy_fraction=self.norm2(h*self.edge)/max(2*energy, 1e-300),
                   energy_balance_error=energy+dissipated-initial_energy,
                   advective_cfl=dt*float(np.max(np.sum(np.abs(v), axis=0)))/(2*np.pi/self.n),
                   diffusion_number=self.nu*dt*3*self.cutoff**2)
        if not all(math.isfinite(value) for value in row.values()):
            raise FloatingPointError("Nonfinite numerical state; this is not evidence of PDE blowup")
        return row


def run(field, points=32, dt=.01, end=4., viscosity=.05, epsilon=.5, sample=.1):
    if min(dt, end, sample) <= 0 or not all(map(math.isfinite, (dt, end, sample))):
        raise ValueError("time settings must be positive and finite")
    steps = round(end/dt)
    stride = round(sample/dt)
    if steps < 1 or stride < 1 or abs(steps*dt-end) > 1e-10 or abs(stride*dt-sample) > 1e-10:
        raise ValueError("end and sample must be positive integer multiples of dt")
    solver = PeriodicFlow(points, viscosity)
    h = solver.initial(field, epsilon)
    e0 = solver.norm2(h)/2
    dissipated = 0.
    rows = [solver.measure(h, 0, dt, e0, dissipated)]
    started = time.perf_counter()
    for step in range(1, steps+1):
        h, loss = solver.step(h, dt)
        dissipated += loss
        if step % stride == 0 or step == steps:
            rows.append(solver.measure(h, step*dt, dt, e0, dissipated))
            if rows[-1]["advective_cfl"] > 1 or rows[-1]["diffusion_number"] > 2:
                raise RuntimeError("Conservative timestep guard exceeded; reduce dt")
    report = dict(field=field, points=points, dt=dt, end=end, viscosity=viscosity,
                  external_force="f(x,t)=0 identically",
                  epsilon=epsilon if field == "C3D" else 0., retained_component_cutoff=solver.cutoff,
                  sample_interval=sample, wall_seconds=time.perf_counter()-started, samples=rows)
    report["summary"] = dict(
        initial_energy=e0, final_energy=rows[-1]["energy"],
        initial_max_gradient=rows[0]["max_gradient"], final_max_gradient=rows[-1]["max_gradient"],
        peak_sampled_max_gradient=max(r["max_gradient"] for r in rows),
        time_of_sampled_gradient_peak=max(rows, key=lambda r:r["max_gradient"])["time"],
        peak_sampled_low_mode_energy=max(r["low_mode_energy"] for r in rows),
        max_rms_divergence=max(r["rms_divergence"] for r in rows),
        max_abs_energy_balance_error=max(abs(r["energy_balance_error"]) for r in rows),
        max_sampled_edge_energy_fraction=max(r["edge_energy_fraction"] for r in rows),
        max_sampled_advective_cfl=max(r["advective_cfl"] for r in rows),
        max_sampled_diffusion_number=max(r["diffusion_number"] for r in rows))
    return report, h


def relative_field_error(coarse, fine):
    """Full Fourier L2 difference on the finer basis, including omitted modes."""
    nc, nf = coarse.shape[1], fine.shape[1]
    if nc > nf:
        return relative_field_error(fine, coarse)
    lifted = np.zeros_like(fine)
    ids = np.rint(np.fft.fftfreq(nc)*nc).astype(int) % nf
    lifted[:, ids[:, None], ids[None, :], :nc//2+1] = coarse
    weights = np.full(nf//2+1, 2.)
    weights[[0, -1]] = 1
    return float(np.sqrt(np.sum(np.abs(lifted-fine)**2*weights)
                         /np.sum(np.abs(fine)**2*weights)))


def suite(output):
    configs = [("A", 16, .02), ("B", 16, .02), ("C", 32, .01),
               ("C3D", 32, .02), ("C3D", 32, .01), ("C3D", 48, .01)]
    document = dict(equation="unforced incompressible Navier-Stokes", domain="[0,2*pi)^3 periodic",
                    external_force="f(x,t)=0 identically; C3D changes initial velocity only",
                    method="Fourier Galerkin with strict 2/3 truncation, pressure projection, RK4",
                    observation="spherical low-pass |k|<=1; observation only, not the evolution grid",
                    interpretation="finite-time numerical results; no general regularity or blowup proof",
                    numpy_version=np.__version__, runs=[], comparisons={})
    hats = {}
    output.parent.mkdir(parents=True, exist_ok=True)
    for name, n, dt in configs:
        print(f"Running {name}: {n}^3, dt={dt}, nu=.05, t=0..4", flush=True)
        report, h = run(name, n, dt)
        hats[name, n, dt] = h
        document["runs"].append(report)
        output.write_text(json.dumps(document, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(json.dumps(report["summary"]), flush=True)
    document["comparisons"] = dict(
        A_vs_negative_B_final_relative_l2=relative_field_error(hats["A",16,.02], -hats["B",16,.02]),
        C3D_dt_002_vs_001_final_relative_l2=relative_field_error(hats["C3D",32,.02], hats["C3D",32,.01]),
        C3D_N32_vs_N48_final_relative_l2=relative_field_error(hats["C3D",32,.01], hats["C3D",48,.01]))
    output.write_text(json.dumps(document, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(document["comparisons"]), flush=True)
    # Save only the final spectral reference alongside the JSON for any justified refinement.
    np.savez_compressed(output.with_suffix(".reference.npz"), C3D_N48_dt001=hats["C3D",48,.01],
                        C3D_N32_dt001=hats["C3D",32,.01], C_N32_dt001=hats["C",32,.01])
    refine_saved(output)


def refine_saved(output):
    """One targeted refinement if the saved 32/48 field difference exceeds .1%."""
    document = json.loads(output.read_text(encoding="utf-8"))
    comparison = document["comparisons"]["C3D_N32_vs_N48_final_relative_l2"]
    document["refinement_rule"] = "Add N=64 if the N=32/48 final-field relative L2 difference exceeds 0.001."
    if comparison > .001 and not any(r["field"] == "C3D" and r["points"] == 64 for r in document["runs"]):
        reference = next(r for r in document["runs"] if r["field"] == "C3D" and r["points"] == 48)
        with np.load(output.with_suffix(".reference.npz")) as saved:
            hats = {key: saved[key].copy() for key in saved.files}
        print("Refining C3D to 64^3 because the spatial comparison exceeded 0.1%", flush=True)
        report, h = run("C3D", 64, reference["dt"], reference["end"],
                        reference["viscosity"], reference["epsilon"], reference["sample_interval"])
        document["runs"].append(report)
        document["comparisons"]["C3D_N48_vs_N64_final_relative_l2"] = relative_field_error(hats["C3D_N48_dt001"], h)
        hats["C3D_N64_dt001"] = h
        np.savez_compressed(output.with_suffix(".reference.npz"), **hats)
        print(json.dumps(report["summary"]), flush=True)
        print(json.dumps(document["comparisons"]), flush=True)
    output.write_text(json.dumps(document, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", action="store_true")
    parser.add_argument("--refine-saved", action="store_true", help="Resume the single refinement gate without repeating the suite")
    parser.add_argument("--field", choices=("A", "B", "C", "C3D"), default="C3D")
    parser.add_argument("--points", type=int, default=32)
    parser.add_argument("--dt", type=float, default=.01)
    parser.add_argument("--end", type=float, default=4.)
    parser.add_argument("--viscosity", type=float, default=.05)
    parser.add_argument("--epsilon", type=float, default=.5)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent/"results"/"hidden-flow.json")
    args = parser.parse_args()
    if args.refine_saved:
        refine_saved(args.output)
    elif args.suite:
        suite(args.output)
    else:
        report, _ = run(args.field, args.points, args.dt, args.end, args.viscosity, args.epsilon)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
