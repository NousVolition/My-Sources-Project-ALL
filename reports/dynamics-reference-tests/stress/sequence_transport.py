"""Test the user's exact digit sequence at equal RK4 work, with refinement."""
from pathlib import Path
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from adaptive_transport import T, MODE, PLAN as MODEL_PLAN, field, rk4_step, diagnostics, write_csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator, ScalarFormatter

ROOT = Path(__file__).resolve().parent
PLAN = json.loads((ROOT/"sequence_plan.json").read_text())


def decode_pairs(text):
    tokens = text.split()
    if not tokens or any(len(token) != 2 or any(c not in "123" for c in token) for token in tokens):
        raise ValueError("Expected two-digit pairs drawn from 1, 2 and 3")
    return np.array([int(c) for token in tokens for c in token], dtype=float)


def multipliers(schedule):
    original = decode_pairs(PLAN["user_pairs"])
    if schedule == "user":
        return original
    if schedule == "reverse":
        return original[::-1]
    if schedule == "uniform":
        return np.full(len(original), np.mean(original))
    raise ValueError(schedule)


def schedule_times(schedule, cycles):
    if not isinstance(cycles, int) or cycles <= 0:
        raise ValueError("cycles must be a positive integer")
    weights = multipliers(schedule)
    offsets = np.r_[0., np.cumsum(weights)]/weights.sum()
    # Anchor each cycle independently; do not accumulate thousands of rounded h values.
    times = ((np.arange(cycles)[:, None]+offsets[:-1][None, :])*(T/cycles)).ravel()
    return np.r_[times, T]


def integrate_case(profile, schedule, cycles, grid):
    theta, rhs = field(grid, profile)
    times = schedule_times(schedule, cycles)
    states = np.empty((len(times), grid))
    states[0] = np.cos(MODE*theta)
    for j, h in enumerate(np.diff(times)):
        states[j+1] = rk4_step(rhs, times[j], states[j], h)
    return theta, times, states


def run():
    started = time.perf_counter()
    data = ROOT/"data"/"sequence"
    data.mkdir(parents=True, exist_ok=True)
    rows, checks, selected = [], [], {}
    largest = max(PLAN["cycle_counts"])
    for cycles in PLAN["cycle_counts"]:
        grids = [PLAN["primary_grid"]] + (PLAN["finest_level_control_grids"] if cycles == largest else [])
        for grid in grids:
            for profile in PLAN["speed_profiles"]:
                for schedule in PLAN["schedules"]:
                    begin = time.perf_counter()
                    theta, times, states = integrate_case(profile, schedule, cycles, grid)
                    metric, arrays = diagnostics(theta, times, states, profile)
                    step = np.diff(times)
                    name = f"{profile}_{schedule}_c{cycles}_n{grid}"
                    rows.append(dict(case=name, speed_profile=profile, schedule=schedule, cycles=cycles, grid=grid,
                        base_interval=T/(35*cycles), steps=len(step), rhs_evaluations=4*len(step),
                        minimum_step=float(step.min()), maximum_step=float(step.max()),
                        elapsed_seconds=time.perf_counter()-begin, **metric))
                    indices = np.arange(0, len(times), 20)
                    np.savez_compressed(data/f"{name}.npz", times=times, theta=theta,
                        snapshot_indices=indices, q_at_cycle_boundaries=states[indices], **arrays)
                    if grid == 64 and cycles == largest and schedule == "user":
                        selected[profile] = (times.copy(), states.copy())
                    print(f"Completed {name}: profile error {metric['maximum_sampled_profile_error']:.5g}", flush=True)

    repeats = []
    for profile, (old_times, old_states) in selected.items():
        begin = time.perf_counter()
        _, times, states = integrate_case(profile, "user", largest, 64)
        passed = bool(np.array_equal(times, old_times) and np.array_equal(states, old_states))
        repeats.append(dict(speed_profile=profile, grid=64, cycles=largest, schedule="user",
                            exact_same_environment_repeat=passed, elapsed_seconds=time.perf_counter()-begin))
        checks.append(dict(name=f"{profile}: selected full-state repeat", passed=passed))
    for profile in PLAN["speed_profiles"]:
        for cycles in PLAN["cycle_counts"]:
            own = [r for r in rows if r["grid"] == 64 and r["cycles"] == cycles and r["speed_profile"] == profile]
            checks.append(dict(name=f"{profile}, {cycles} cycles: matched work", passed=len({r["rhs_evaluations"] for r in own}) == 1))
        for schedule in PLAN["schedules"]:
            own = sorted([r for r in rows if r["grid"] == 64 and r["schedule"] == schedule and r["speed_profile"] == profile], key=lambda r: r["cycles"])
            checks.append(dict(name=f"{profile}, {schedule}: error decreases on each time refinement", passed=all(
                b["maximum_sampled_profile_error"] < a["maximum_sampled_profile_error"] for a, b in zip(own[:-1], own[1:]))))
            fine = next(r for r in rows if r["grid"] == 128 and r["schedule"] == schedule and r["speed_profile"] == profile)
            coarse = next(r for r in rows if r["grid"] == 32 and r["schedule"] == schedule and r["speed_profile"] == profile)
            checks.append(dict(name=f"{profile}, {schedule}: finest resolved accuracy and aliased control", passed=bool(
                own[-1]["accuracy_passed"] and fine["accuracy_passed"] and coarse["maximum_sampled_profile_error"] > 1)))
    checks.append(dict(name="All mass-budget checks pass independently", passed=all(r["conservation_passed"] for r in rows)))
    checks.append(dict(name="All finite sampled boundedness guards pass independently", passed=all(r["boundedness_guard_passed"] for r in rows)))
    audit_files = [Path(__file__), ROOT/"sequence_plan.json", ROOT/"adaptive_transport.py", ROOT/"adaptive_transport_plan.json"]
    result = dict(cases=rows, checks=checks, checks_passed=sum(c["passed"] for c in checks), checks_total=len(checks),
        decoded_sequence=multipliers("user").astype(int).tolist(), repeats=repeats,
        elapsed_seconds=time.perf_counter()-started,
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in audit_files},
        classification="Deterministic numerical results, not empirical evidence. No optimization claim. All errors are sampled; no continuous-time bound.")
    write_csv(data/"cases.csv", rows)
    (data/"summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    draw(rows)
    print(json.dumps(dict(cases=len(rows), checks_passed=result["checks_passed"], checks_total=len(checks), elapsed_seconds=result["elapsed_seconds"])))
    if not all(c["passed"] for c in checks):
        raise AssertionError("Recorded benchmark checks failed; preserve results and investigate.")


def draw(rows):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    colors = {"user": "#14798c", "reverse": "#ae6b20", "uniform": "#595c96"}
    for schedule in PLAN["schedules"]:
        axes[0].step(np.arange(1, 21), multipliers(schedule), where="mid", label=schedule, color=colors[schedule])
    axes[0].set(title="One 20-step cycle", xlabel="Step number", ylabel="Multiple of the base interval", xticks=[1, 5, 10, 15, 20], yticks=[1, 1.75, 2, 3])
    axes[0].legend(fontsize=8)
    for ax, profile in zip(axes[1:], PLAN["speed_profiles"]):
        for schedule in PLAN["schedules"]:
            own = sorted([r for r in rows if r["grid"] == 64 and r["speed_profile"] == profile and r["schedule"] == schedule], key=lambda r: r["steps"])
            ax.loglog([r["steps"] for r in own], [r["maximum_sampled_profile_error"] for r in own], marker="o" if schedule != "reverse" else "x", linestyle="--" if schedule == "reverse" else "-", color=colors[schedule], label=schedule)
        ax.axhline(MODEL_PLAN["accuracy_budgets"]["maximum_sampled_profile_error"], color=".5", ls=":", label="profile target")
        ax.set(title=profile.replace("-", " → "), xlabel="RK4 steps over the same duration", ylabel="Maximum sampled profile error")
        ax.set_xticks([1600, 6400, 25600])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(NullLocator())
        ax.legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.15)
    fig.suptitle("Your sequence works when sufficiently refined; compare equal work", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(ROOT/"figures"/f"08_sequence_comparison.{ext}", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    run()
