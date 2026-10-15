"""Build the authorized local higher-viscosity comparison from saved results."""
import html
import os
from pathlib import Path

import numpy as np

from diagnostics import ROOT, read, save

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "mpl-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def build():
    p, r = read(ROOT / "protocol.json"), read(ROOT / "result.json")
    rows = r["series"]
    last = rows[-1]
    times = np.array([x["t"] for x in rows])
    amps = np.array([x["amplification"] for x in rows])
    controls = []
    for label, path, key in [
        ("0.01", Path(p["viscosity01_directory"]) / "result.json", "D_relative_initial_reference_l2"),
        ("0.001", Path(p["original_study"]) / "runs/perturb-n128-a0.001-s202/result.json", "D_relative_initial_baseline_l2")]:
        raw = read(path)["series"]
        available = [x for x in raw if x["t"] <= last["t"] + 1e-12]
        controls.append((label, available, key, raw[0][key]))
    fig, ax = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    ax[0].plot(times, amps, marker=".", label="Viscosity 0.1")
    for label, series, key, norm in controls:
        ax[0].plot([x["t"] for x in series], [x[key]/norm for x in series], marker=".", label=f"Viscosity {label}")
    ax[0].axhline(1, color="black", ls=":", label="Starting disturbance")
    if len(rows) == 1:
        ax[0].set_ylim(.98, 1.02)
    ax[0].set(title="Disturbance relative to its starting size", ylabel="D(t) / D(0)")
    for key, label in [("difference_energy_transfer_from_reference", "Transfer from reference flow"),
                       ("difference_energy_viscous_loss_rate", "Viscous removal"),
                       ("difference_energy_net_rate", "Net rate")]:
        ax[1].plot(times, [x[key] for x in rows], marker=".", label=label)
    ax[1].axhline(0, color="black", ls=":")
    ax[1].set(title="Why the disturbance changes", ylabel="Difference-energy rate")
    for name, label in [("reference", "Undisturbed"), ("disturbed", "Disturbed")]:
        series = [x["reference_diagnostics"] if name == "reference" else x for x in rows]
        ax[2].plot(times, [x["energy"]/series[0]["energy"] for x in series], marker=".", label=label)
    ax[2].set(title="Flow energy is a separate measurement", ylabel="Energy / initial energy")
    for axis in ax:
        axis.set_xlabel("Model time")
        axis.grid(alpha=.2)
        axis.legend(fontsize=8)
    fig.suptitle("Viscosity 0.1: same 0.1% disturbance, seed 202")
    fig.savefig(ROOT / "comparison.png", dpi=140)
    plt.close(fig)
    direction = "No evolved measurement yet" if len(rows) == 1 else ("smaller than at the start" if amps[-1] < 1 else "larger than at the start")
    current = "shrinking" if last["difference_energy_net_rate"] < 0 else "growing"
    summary = {"status": r["status"], "saved_time": last["t"], "saved_pairs": len(rows),
               "amplification": float(amps[-1]), "direction_since_start": direction,
               "instantaneous_direction": current, "instantaneous_log_norm_rate": last["instantaneous_log_norm_rate"],
               "failed_disturbed_spatial_screens": sum(not x["resolved_screen"] for x in rows),
               "failed_reference_spatial_screens": sum(not x["reference_diagnostics"]["resolved_screen"] for x in rows),
               "max_difference_energy_budget_residual": max(abs(x["difference_energy_budget_relative_residual"]) for x in rows),
               "comparison_at_latest_time": {label: series[-1][key]/norm for label, series, key, norm in controls
                                             if series and abs(series[-1]["t"]-last["t"]) < 1e-12},
               "interpretation": "Single finite-amplitude seed; no grid/time convergence or physical/asymptotic stability claim."}
    save(ROOT / "summary.json", summary)
    text = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Viscosity 0.1 disturbance comparison</title>
    <style>body{{font:17px/1.6 system-ui;color:#193447;max-width:1150px;margin:35px auto;padding:0 22px;background:#f6f9fc}}section{{background:white;border:1px solid #ccd9e3;padding:22px;margin:18px 0;border-radius:10px}}img{{width:100%}}a{{color:#08778a}}</style></head><body>
    <h1>Does higher viscosity make the disturbance shrink?</h1>
    <section><p><b>{html.escape(r['status'])}: saved through {last['t']:.2f} of 0.40</b>, with {len(rows)} paired outputs.</p>
    <p>This test uses viscosity <b>0.1</b>, ten times 0.01, with the same starting flow and a <b>0.1% disturbance, seed 202</b>. Both new flows use grid 128, the same mean, domain, zero force, Heun method and time step as the earlier pair. They run in one process; the existing tests continue.</p>
    <p><b>{direction}.</b> The disturbance is {amps[-1]:.6g} times its starting size. Below one means net shrinking since the start; above one means net growth. At this saved instant its calculated tendency is <b>{current}</b>. An initial tendency alone does not establish a persistent trend.</p></section>
    <section><img src="comparison.png" alt="Disturbance comparison, difference-energy rates and flow energy"><p>The difference is the L2 velocity separation divided by the initial undisturbed-flow norm. Transfer minus viscous removal gives the instantaneous difference-energy rate. Existing lower-viscosity curves are shown only where saved times overlap.</p></section>
    <section><h2>What this test can establish</h2><p>It measures this one disturbance over the existing finite time window. Raw strain is nonsmooth across periodic joins and initial maxima are outside the central tube. {summary['failed_disturbed_spatial_screens']} disturbed and {summary['failed_reference_spatial_screens']} reference spatial screens fail among the {len(rows)} saved pairs. Passing individual screens is not a convergence test. Finer-grid and half-step pairs have not been run for this higher-viscosity case. These results do not establish a physical stability threshold, an asymptotic Lyapunov exponent or pond behavior.</p>
    <p>Both fields, their hashes and energy budgets are saved. The disturbance norm and transfer identity are checked independently, and nonfinite fields, excessive CFL, mean drift, divergence and budget failures stop this pair while preserving evidence.</p>
    <p><a href="protocol.json">Protocol</a> · <a href="summary.json">Current summary</a> · <a href="result.json">Saved measurements</a> · <a href="initial-audit.json">Initial checks</a></p></section></body></html>'''
    tmp = ROOT / "report.tmp"
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, ROOT / "report.html")


if __name__ == "__main__":
    build()
