"""Synthetic ground-to-pond coupling, separated from numerical integration error."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adaptive_transport import rk4_step

ROOT = Path(__file__).resolve().parent
PLAN = json.loads((ROOT/"pond_vibration_plan.json").read_text())
L, H, G = PLAN["length_m"], PLAN["depth_m"], PLAN["gravity_m_s2"]
DRAG, RHO = PLAN["linear_drag_s_inverse"], PLAN["density_kg_m3"]
K = np.pi/L
OMEGA = np.sqrt(G*H)*K
START = PLAN["bed_pulse_start_s"]


def bed(t, amplitude, width):
    t = np.asarray(t, dtype=float)
    inside = (t >= START) & (t <= START+width)
    s = np.clip((t-START)/width, 0, 1)
    a = np.pi*s
    b = np.where(inside, amplitude*np.sin(a)**6, 0.)
    db = np.where(inside, amplitude*6*np.pi/width*np.sin(a)**5*np.cos(a), 0.)
    return b, db


def rhs(t, y, amplitude, width):
    # eta=A cos(kx), u=U sin(kx), eta_t+H*u_x=b_t, u_t+g*eta_x=-drag*u.
    a, u = y[:2]
    _, db = bed(t, amplitude, width)
    work_rate = RHO*L*G*a*db/2
    dissipation_rate = RHO*L*H*DRAG*u*u/2
    return np.array([db-H*K*u, G*K*a-DRAG*u, work_rate, dissipation_rate])


def reference(t, amplitude, width, nodes=64):
    """Analytic damped-wave Green function integrated by Gauss quadrature."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    z, weights = leggauss(nodes)
    upper = np.maximum(START, np.minimum(t, START+width))
    half = (upper-START)/2
    s = START+half[:, None]*(1+z[None, :])
    tau = np.maximum(t[:, None]-s, 0)
    alpha = DRAG/2
    wd = np.sqrt(OMEGA*OMEGA-alpha*alpha)
    decay = np.exp(-alpha*tau)
    sine = np.sin(wd*tau)
    g00 = decay*(np.cos(wd*tau)+alpha/wd*sine)
    g10 = decay*G*K/wd*sine
    _, db = bed(s, amplitude, width)
    a = half*np.sum(weights[None, :]*g00*db, axis=1)
    u = half*np.sum(weights[None, :]*g10*db, axis=1)
    return np.stack((a, u), axis=1)


def integrate(amplitude, width, step):
    count = round(PLAN["duration_s"]/step)
    assert np.isclose(count*step, PLAN["duration_s"])
    t = np.arange(count+1)*step
    state = np.zeros((count+1, 4))
    for j in range(count):
        state[j+1] = rk4_step(lambda s, y: rhs(s, y, amplitude, width), t[j], state[j], step)
    return t, state


def energy(state):
    return RHO*L/4*(G*state[:, 0]**2+H*state[:, 1]**2)


def run():
    started = time.perf_counter()
    data = ROOT/"data"/"pond-vibration"
    data.mkdir(parents=True, exist_ok=True)
    rows, response, plots = [], [], {}
    dense_t = np.linspace(0, PLAN["duration_s"], round(PLAN["duration_s"]/PLAN["reference_diagnostic_step_s"])+1)
    for amplitude in PLAN["assumed_bed_mode_amplitudes_m"]:
        for width in PLAN["pulse_durations_s"]:
            ref = reference(dense_t, amplitude, width, PLAN["reference_quadrature_nodes"])
            ref32 = reference(dense_t, amplitude, width, 32)
            peak = float(np.max(np.abs(ref[:, 0])))
            peak_energy = float(energy(ref).max())
            b, db = bed(dense_t, amplitude, width)
            name = f"B{amplitude:g}_width{width:g}"
            np.savez_compressed(data/f"reference_{name}.npz", time_s=dense_t, surface_mode_m=ref[:, 0], velocity_mode_m_s=ref[:, 1], bed_mode_m=b, bed_velocity_mode_m_s=db)
            response.append(dict(amplitude_m=amplitude, width_s=width, peak_surface_mode_m=peak,
                peak_surface_mode_after_pulse_m=float(np.max(np.abs(ref[dense_t >= START+width, 0]))),
                peak_velocity_mode_m_s=float(np.max(np.abs(ref[:, 1]))),
                quadrature_32_vs_64_max_surface_difference_m=float(np.max(np.abs(ref[:, 0]-ref32[:, 0])))))
            if amplitude == max(PLAN["assumed_bed_mode_amplitudes_m"]):
                plots[width] = (dense_t, b, ref[:, 0])
            for h in PLAN["rk4_steps_s"]:
                t, state = integrate(amplitude, width, h)
                exact = reference(t, amplitude, width)
                error = np.max(np.abs(state[:, 0]-exact[:, 0]))
                energy_residual = energy(state)-state[:, 2]+state[:, 3]
                # Modal volume conservation is analytical, not a general PDE grid test.
                x = np.linspace(0, L, 65)
                bed_now, _ = bed(t, amplitude, width)
                volume = np.trapezoid((state[:, 0]-bed_now)[:, None]*np.cos(K*x)[None, :], x, axis=1)
                row = dict(amplitude_m=amplitude, width_s=width, step_s=h,
                    steps=len(t)-1, maximum_sampled_surface_error_m=float(error),
                    error_over_reference_peak=float(error/peak),
                    maximum_sampled_energy_residual_J_per_m=float(np.abs(energy_residual).max()),
                    energy_residual_over_reference_peak=float(np.abs(energy_residual).max()/peak_energy),
                    maximum_sampled_volume_error_m2=float(np.abs(volume).max()),
                    accuracy_passed=bool(error/peak <= PLAN["fine_relative_surface_error_budget"]),
                    energy_budget_passed=bool(np.abs(energy_residual).max()/peak_energy <= PLAN["fine_relative_energy_residual_budget"]))
                rows.append(row)
                np.savez_compressed(data/f"rk4_{name}_h{h:g}.npz", time_s=t, state=state,
                    reference_state=exact, energy_residual_J_per_m=energy_residual, volume_error_m2=volume)
    zero_t, zero_state = integrate(0, min(PLAN["pulse_durations_s"]), min(PLAN["rk4_steps_s"]))
    np.savez_compressed(data/"zero_motion_control.npz", time_s=zero_t, state=zero_state)
    fine = [r for r in rows if r["step_s"] == min(PLAN["rk4_steps_s"])]
    checks = [dict(name="Zero bed motion leaves resting water exactly at rest", passed=bool(np.all(zero_state == 0))),
        dict(name="Finest time step meets each surface-error budget", passed=all(r["accuracy_passed"] for r in fine)),
        dict(name="Finest time step meets each energy budget", passed=all(r["energy_budget_passed"] for r in fine)),
        dict(name="All modal volume budgets close within 1e-15 square metres", passed=all(r["maximum_sampled_volume_error_m2"] < 1e-15 for r in rows)),
        dict(name="Independent quadrature refinement below 1e-12 metres", passed=all(r["quadrature_32_vs_64_max_surface_difference_m"] < 1e-12 for r in response))]
    result = dict(status=PLAN["status"], cases=rows, responses=response, checks=checks,
        checks_passed=sum(c["passed"] for c in checks), checks_total=len(checks),
        natural_period_s=float(2*np.pi/OMEGA), depth_wavenumber=float(K*H),
        shallow_frequency_relative_difference_from_finite_depth=float(np.sqrt(K*H/np.tanh(K*H))-1),
        elapsed_seconds=time.perf_counter()-started, environment=dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__),
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), ROOT/"pond_vibration_plan.json", ROOT/"adaptive_transport.py")})
    with (data/"cases.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (data/"summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for width, (t, b, a) in plots.items():
        axes[0].plot(t, b*1e6, label=f"{width:g} s pulse")
        axes[1].plot(t, a*1e6, label=f"{width:g} s pulse")
        own = [r for r in rows if r["width_s"] == width and r["amplitude_m"] == max(PLAN["assumed_bed_mode_amplitudes_m"])]
        axes[2].loglog([r["step_s"] for r in own], [r["error_over_reference_peak"] for r in own], "o-", label=f"{width:g} s pulse")
    axes[0].set(xlim=(.8, 3.3), xlabel="Time (s)", ylabel="Bed-mode amplitude (µm)", title="Assumed movement at the pond")
    axes[1].set(xlabel="Time (s)", ylabel="Surface-mode amplitude (µm)", title="Predicted response in this model")
    axes[2].axhline(PLAN["fine_relative_surface_error_budget"], color="k", ls=":", label="0.1% error target")
    axes[2].set(xlabel="RK4 time step (s)", ylabel="Numerical error / reference peak", title="Short pulses demand closer timing")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=.15)
    fig.suptitle("Illustrative pond: no footstep or rock strength has been calibrated", fontsize=13)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(ROOT/"figures"/f"10_pond_vibration.{ext}", dpi=160)
    plt.close(fig)
    print(json.dumps(dict(cases=len(rows), checks_passed=result["checks_passed"], checks_total=len(checks), elapsed_seconds=result["elapsed_seconds"])))
    if not all(c["passed"] for c in checks):
        raise AssertionError("A predeclared pond benchmark check failed; inspect recorded results")


if __name__ == "__main__":
    run()
