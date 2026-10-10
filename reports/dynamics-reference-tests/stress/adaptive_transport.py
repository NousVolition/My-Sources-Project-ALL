"""Variable-step transport: exact reference, both speed orders, aliased control."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
PLAN = json.loads((ROOT / "adaptive_transport_plan.json").read_text())
T = PLAN["duration"]
MODE = PLAN["mode"]


def speed(t, profile):
    modulation = np.sin(np.pi * np.asarray(t) / T) ** 2
    low, high = PLAN["speed_minimum"], PLAN["speed_maximum"]
    return low + (high-low)*modulation if profile == "slow-fast-slow" else high - (high-low)*modulation


def rotation(t, profile):
    t = np.asarray(t)
    integral = t/2 - T/(4*np.pi)*np.sin(2*np.pi*t/T)
    low, high = PLAN["speed_minimum"], PLAN["speed_maximum"]
    return low*t + (high-low)*integral if profile == "slow-fast-slow" else high*t - (high-low)*integral


def exact(theta, times, profile):
    return np.cos(MODE*(theta[None, :] - rotation(times, profile)[:, None]))


def field(grid, profile):
    theta = 2*np.pi*np.arange(grid)/grid
    modes = np.fft.fftfreq(grid, d=1/grid)
    def rhs(t, q):
        return -speed(t, profile)*np.fft.ifft(1j*modes*np.fft.fft(q)).real
    return theta, rhs


def rk4_step(rhs, t, q, h):
    # Nonautonomous stages must evaluate the speed at their own stage times.
    k1 = rhs(t, q)
    k2 = rhs(t+h/2, q+h*k1/2)
    k3 = rhs(t+h/2, q+h*k2/2)
    k4 = rhs(t+h, q+h*k3)
    return q+h*(k1+2*k2+2*k3+k4)/6


def manual_integrate(grid, profile, order):
    theta, rhs = field(grid, profile)
    small, large = PLAN["manual_steps"]
    steps = (large, small, large) if order == "large-small-large" else (small, large, small)
    boundaries = (0, *PLAN["manual_switch_times"], T)
    q = np.cos(MODE*theta)
    times, states = [0.0], [q.copy()]
    for a, b, h in zip(boundaries[:-1], boundaries[1:], steps):
        count = round((b-a)/h)
        assert abs(count*h-(b-a)) < 1e-12
        for j in range(count):
            q = rk4_step(rhs, a+j*h, q, h)
            times.append(a+(j+1)*h)
            states.append(q.copy())
    return theta, np.array(times), np.array(states)


def diagnostics(theta, times, states, profile):
    reference = exact(theta, times, profile)
    errors = np.max(np.abs(states-reference), axis=1)
    coefficient = 2*np.mean(states*np.exp(-1j*MODE*theta)[None, :], axis=1)
    amplitude_error = np.abs(np.abs(coefficient)-1)
    # The output interval and maximum accepted interval resolve phase changes here.
    phase_error = np.degrees(np.unwrap(np.angle(coefficient)) + MODE*rotation(times, profile))
    mass_error = 2*np.pi*np.abs(np.mean(states, axis=1)-np.mean(states[0]))
    maximum = float(np.max(np.abs(states)))
    metrics = {
        "maximum_sampled_profile_error": float(errors.max()),
        "endpoint_profile_error": float(errors[-1]),
        "maximum_sampled_amplitude_error": float(amplitude_error.max()),
        "maximum_sampled_accumulated_phase_error_degrees": float(np.max(np.abs(phase_error))),
        "maximum_sampled_mass_error": float(mass_error.max()),
        "maximum_sampled_fluctuation": maximum,
        "accuracy_passed": bool(all(value <= PLAN["accuracy_budgets"][key] for key, value in (
            ("maximum_sampled_profile_error", errors.max()),
            ("maximum_sampled_amplitude_error", amplitude_error.max()),
            ("maximum_sampled_accumulated_phase_error_degrees", np.max(np.abs(phase_error)))))),
        "conservation_passed": bool(mass_error.max() <= PLAN["conservation_budget"]),
        "boundedness_guard_passed": bool(np.isfinite(states).all() and maximum <= PLAN["boundedness_guard"]),
    }
    return metrics, dict(profile_error=errors, amplitude_error=amplitude_error,
                         phase_error_degrees=phase_error, mass_error=mass_error)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run():
    start = time.perf_counter()
    data = ROOT/"data"/"adaptive"
    data.mkdir(parents=True, exist_ok=True)
    sampled = np.linspace(0, T, round(T/PLAN["diagnostic_interval"])+1)
    rows, plots, checks = [], {}, []
    for profile in PLAN["speed_profiles"]:
        for grid in PLAN["grids"]:
            theta, rhs = field(grid, profile)
            for rtol in PLAN["relative_tolerances"]:
                name = f"{profile}_n{grid}_tol{rtol:g}"
                begin = time.perf_counter()
                solution = solve_ivp(rhs, (0, T), np.cos(MODE*theta), method=PLAN["adaptive_method"],
                    rtol=rtol, atol=rtol*PLAN["absolute_tolerance_factor"],
                    max_step=PLAN["maximum_step"], dense_output=True)
                if not solution.success:
                    raise RuntimeError(solution.message)
                states = solution.sol(sampled).T
                metric, arrays = diagnostics(theta, sampled, states, profile)
                steps = np.diff(solution.t)
                mid = (solution.t[1:]+solution.t[:-1])/2
                velocities = speed(mid, profile)
                row = dict(case=name, speed_profile=profile, grid=grid, method="DOP853", schedule="error-controlled",
                    rtol=rtol, accepted_steps=len(steps), rhs_evaluations=solution.nfev,
                    minimum_step=float(steps.min()), maximum_step=float(steps.max()),
                    median_step_at_speed_below_1=float(np.median(steps[velocities < 1])),
                    median_step_at_speed_above_3=float(np.median(steps[velocities > 3])),
                    elapsed_seconds=time.perf_counter()-begin, **metric)
                rows.append(row)
                np.savez_compressed(data/f"{name}.npz", theta=theta, times=sampled, q=states,
                    accepted_times=solution.t, accepted_q=solution.y.T, **arrays)
                if grid == 64 and rtol == min(PLAN["relative_tolerances"]):
                    plots[profile] = (mid, steps, arrays)
        for order in PLAN["manual_orders"]:
            begin = time.perf_counter()
            grid = PLAN["manual_rk4_grid"]
            theta, times, states = manual_integrate(grid, profile, order)
            metric, arrays = diagnostics(theta, times, states, profile)
            steps = np.diff(times)
            velocities = speed((times[1:]+times[:-1])/2, profile)
            name = f"{profile}_{order}"
            rows.append(dict(case=name, speed_profile=profile, grid=grid, method="RK4", schedule=order,
                rtol="", accepted_steps=len(steps), rhs_evaluations=4*len(steps),
                minimum_step=float(steps.min()), maximum_step=float(steps.max()),
                median_step_at_speed_below_1=float(np.median(steps[velocities < 1])),
                median_step_at_speed_above_3=float(np.median(steps[velocities > 3])),
                elapsed_seconds=time.perf_counter()-begin, **metric))
            np.savez_compressed(data/f"{name}.npz", theta=theta, times=times, q=states, **arrays)

    for profile in PLAN["speed_profiles"]:
        own = [r for r in rows if r["speed_profile"] == profile]
        for grid in PLAN["grids"]:
            loose, tight = sorted([r for r in own if r["grid"] == grid and r["method"] == "DOP853"], key=lambda r: -r["rtol"])
            checks.append(dict(name=f"{profile}, grid {grid}: expected accuracy or aliasing", passed=(
                tight["maximum_sampled_profile_error"] > 1 if grid == 32 else
                tight["accuracy_passed"] and tight["maximum_sampled_profile_error"] < loose["maximum_sampled_profile_error"])))
        manual = [r for r in own if r["method"] == "RK4"]
        checks.append(dict(name=f"{profile}: manual schedules have matched evaluation counts", passed=manual[0]["rhs_evaluations"] == manual[1]["rhs_evaluations"]))
    checks.append(dict(name="Mass budgets pass separately", passed=all(r["conservation_passed"] for r in rows)))
    checks.append(dict(name="Finite sampled boundedness guard passes separately", passed=all(r["boundedness_guard_passed"] for r in rows)))
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), ROOT/"adaptive_transport_plan.json")}
    summary = dict(cases=rows, checks=checks, checks_passed=sum(r["passed"] for r in checks), checks_total=len(checks),
        elapsed_seconds=time.perf_counter()-start, source_sha256=hashes,
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        sampling="Adaptive cases: dense output every 0.005. Manual cases: every accepted step, no interpolated claim. Maxima are sampled maxima, not rigorous continuous bounds.",
        classification="Deterministic numerical benchmark, no empirical measurements or statistical confidence intervals.")
    write_csv(data/"cases.csv", rows)
    (data/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for col, profile in enumerate(PLAN["speed_profiles"]):
        mid, steps, arrays = plots[profile]
        axes[0, col].plot(sampled, speed(sampled, profile), color="#267d90")
        axes[0, col].set(title=profile.replace("-", " → "), ylabel="Prescribed angular speed")
        axes[1, col].plot(mid, steps, color="#a74b29", linewidth=1)
        axes[1, col].set(xlabel="Model time", ylabel="Accepted time step")
    fig.suptitle("Adaptive intervals follow the error budget in both directions\nMode 17 on grid 64; DOP853 relative tolerance 10⁻⁹", fontsize=13)
    for ax in axes.flat:
        ax.grid(alpha=.18)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(ROOT/"figures"/f"07_adaptive_intervals.{ext}", dpi=150)
    plt.close(fig)
    print(json.dumps(dict(cases=len(rows), checks_passed=summary["checks_passed"], checks_total=len(checks), elapsed_seconds=summary["elapsed_seconds"])))
    if not all(r["passed"] for r in checks):
        raise AssertionError("Benchmark expectations failed; inspect recorded results without widening budgets.")


if __name__ == "__main__":
    run()
