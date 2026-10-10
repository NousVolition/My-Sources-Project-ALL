"""Exact circular limit-cycle benchmark from the supplied screenshot."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from core import integrate
PLAN = json.loads((ROOT/"limit_cycle_plan.json").read_text())


def rhs(state):
    x, y = state
    growth = 1-x*x-y*y
    return np.array([growth*x-y, x+growth*y])


def exact_radius(t, r0):
    t = np.asarray(t, dtype=float)
    if r0 < 0 or np.any(t < 0):
        raise ValueError("This implementation specifies nonnegative radius and forward time")
    if r0 == 0:
        return np.zeros_like(t)
    # This form is valid on both sides of r=1 and avoids exp(+2t) overflow.
    return r0/np.sqrt(r0*r0+(1-r0*r0)*np.exp(-2*t))


def exact_state(t, r0, theta0):
    angle = theta0+np.asarray(t)
    radius = exact_radius(t, r0)
    return np.stack((radius*np.cos(angle), radius*np.sin(angle)), axis=-1)


def metrics(t, states, r0, theta0):
    radius = np.linalg.norm(states, axis=1)
    radial_error = np.abs(radius-exact_radius(t, r0))
    state_error = np.linalg.norm(states-exact_state(t, r0, theta0), axis=1)
    phase_error = None if r0 == 0 else np.degrees(np.unwrap(np.arctan2(states[:, 1], states[:, 0]))-theta0-t)
    result = dict(maximum_sampled_radial_error=float(radial_error.max()),
        maximum_sampled_phase_error_degrees=None if phase_error is None else float(np.abs(phase_error).max()),
        maximum_sampled_state_error=float(state_error.max()), endpoint_radius=float(radius[-1]),
        maximum_sampled_radius=float(radius.max()), phase_applicable=r0 > 0,
        radial_accuracy_passed=bool(radial_error.max() <= PLAN["accuracy_budgets"]["maximum_sampled_radial_error"]),
        phase_accuracy_passed=None if phase_error is None else bool(np.abs(phase_error).max() <= PLAN["accuracy_budgets"]["maximum_sampled_phase_error_degrees"]),
        boundedness_guard_passed=bool(np.isfinite(states).all() and radius.max() <= PLAN["boundedness_guard_radius"]))
    arrays = dict(radius=radius, radial_error=radial_error, state_error=state_error)
    if phase_error is not None:
        arrays["phase_error_degrees"] = phase_error
    return result, arrays


def run():
    started = time.perf_counter()
    data = ROOT/"data"/"limit-cycle"
    data.mkdir(parents=True, exist_ok=True)
    rows, traces = [], {}
    for duration, radii in ((PLAN["duration"], PLAN["initial_radii"]),
                            (PLAN["phase_drift_duration"], [PLAN["phase_drift_initial_radius"]])):
        for r0 in radii:
            for method in PLAN["methods"]:
                for h in PLAN["steps"]:
                    theta0 = PLAN["initial_angle"]
                    initial = exact_state(0, r0, theta0)
                    t, states = integrate(rhs, initial, duration, h, method)
                    measure, arrays = metrics(t, states, r0, theta0)
                    name = f"{method}_r{r0:g}_T{duration:g}_h{h:g}"
                    row = dict(case=name, method=method, initial_radius=r0, initial_angle=theta0, duration=duration, step=h, steps=len(t)-1, **measure)
                    rows.append(row)
                    np.savez_compressed(data/f"{name}.npz", t=t, state=states, **arrays)
                    if duration == PLAN["phase_drift_duration"] and h == max(PLAN["steps"]):
                        traces[method] = (t, arrays)
    checks = []
    checks.append(dict(name="All origin starts remain exactly at the origin", passed=all(r["maximum_sampled_radius"] == 0 for r in rows if r["initial_radius"] == 0)))
    checks.append(dict(name="All finite sampled boundedness guards", passed=all(r["boundedness_guard_passed"] for r in rows)))
    for r0 in PLAN["initial_radii"]:
        fine = next(r for r in rows if r["initial_radius"] == r0 and r["duration"] == PLAN["duration"] and r["method"] == "rk4" and r["step"] == min(PLAN["steps"]))
        checks.append(dict(name=f"Fine RK4 exact-state error r0={r0}", value=fine["maximum_sampled_state_error"],
            passed=fine["maximum_sampled_state_error"] <= PLAN["fine_rk4_state_error_budget"]))
    repeat_initial = exact_state(0, 1, PLAN["initial_angle"])
    t_repeat, states_repeat = integrate(rhs, repeat_initial, PLAN["phase_drift_duration"], min(PLAN["steps"]), "rk4")
    with np.load(data/"rk4_r1_T100_h0.0125.npz") as stored:
        identical = bool(np.array_equal(t_repeat, stored["t"]) and np.array_equal(states_repeat, stored["state"]))
    checks.append(dict(name="Selected fine RK4 100-unit full-state repeat", passed=identical))
    with (data/"cases.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = dict(cases=rows, checks=checks, checks_passed=sum(c["passed"] for c in checks), checks_total=len(checks),
        elapsed_seconds=time.perf_counter()-started, repeat_scope="One fine RK4 100-unit trajectory; not all 84 cases.",
        environment=dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__),
        source_sha256={str(p.relative_to(ROOT.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), ROOT/"limit_cycle_plan.json", ROOT.parent/"core.py")})
    (data/"summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    draw(traces)
    print(json.dumps(dict(cases=len(rows), checks_passed=result["checks_passed"], checks_total=len(checks), elapsed_seconds=result["elapsed_seconds"])))
    if not all(c["passed"] for c in checks):
        raise AssertionError("A declared benchmark expectation failed; preserve results and investigate")


def draw(traces):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    t = np.linspace(0, 12, 1801)
    for r0, angle in ((.2, .3), (.8, 1.5), (1.5, 3), (2, 4.4)):
        states = exact_state(t, r0, angle)
        axes[0].plot(states[:, 0], states[:, 1], label=f"r₀={r0}", linewidth=1.2)
        axes[0].scatter(*states[0], s=18)
    phi = np.linspace(0, 2*np.pi, 361)
    axes[0].plot(np.cos(phi), np.sin(phi), "k--", linewidth=1, label="r=1 cycle")
    axes[0].scatter(0, 0, c="k", s=25, label="origin stays put")
    axes[0].set(aspect="equal", xlabel="x", ylabel="y", title="Exact paths approach the circle")
    axes[0].legend(fontsize=7, loc="lower left")
    for method, (times, arrays) in traces.items():
        axes[1].plot(times, arrays["radius"], label=method)
        axes[2].plot(times, arrays["phase_error_degrees"], label=method)
    axes[1].axhline(1, c="k", ls=":", label="exact radius")
    axes[1].set(xlabel="Model time", ylabel="Radius", title="Numerical radius, h=0.1")
    axes[2].axhline(0, c="k", ls=":")
    axes[2].set(xlabel="Model time", ylabel="Accumulated phase error (degrees)", title="Being on a circle is not enough")
    for ax in axes[1:]:
        ax.legend(fontsize=8)
        ax.grid(alpha=.15)
    fig.suptitle("The screenshot's oscillator: geometry and timing need separate checks", fontsize=13)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(ROOT/"figures"/f"09_limit_cycle.{ext}", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    run()
