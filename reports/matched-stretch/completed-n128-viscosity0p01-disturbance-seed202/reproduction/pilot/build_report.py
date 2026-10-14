"""Refresh a local report from saved observations; no trajectory integration."""
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
    t = np.array([x["t"] for x in rows])
    A = np.array([x["amplification"] for x in rows])
    D = np.array([x["D_relative_initial_reference_l2"] for x in rows])
    low = read(Path(p["original_study"]) / "runs/perturb-n128-a0.001-s202/result.json")
    matched = low["series"][:len(rows)]
    assert len(matched) == len(rows)
    low_A = np.array([x["D_relative_initial_baseline_l2"] / low["series"][0]["D_relative_initial_baseline_l2"] for x in matched])
    fig, ax = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    ax[0].semilogy(t, A, marker=".", label="Viscosity 0.01")
    ax[0].semilogy(t, low_A, marker=".", label="Viscosity 0.001")
    ax[0].axhline(1, color="black", ls=":", label="Starting disturbance")
    ax[0].set(title="Difference from its matching reference", ylabel="D(t) / D(0)")
    ax[1].plot(t, [x["difference_energy_transfer_from_reference"] for x in rows], label="Transfer from reference flow")
    ax[1].plot(t, [x["difference_energy_viscous_loss_rate"] for x in rows], label="Viscous removal")
    ax[1].plot(t, [x["difference_energy_net_rate"] for x in rows], label="Net difference-energy rate")
    ax[1].axhline(0, color="black", ls=":")
    ax[1].set(title="Why the difference changes", ylabel="Instantaneous model rate")
    ax[2].plot(t, [x["Wmax"] for x in rows], label="Disturbed flow")
    ax[2].plot(t, [x["reference_Wmax"] for x in rows], label="Reference flow")
    ax[2].set(title="Peak spin is a separate measurement", ylabel="Peak vorticity")
    for a in ax:
        a.set(xlabel="Model time")
        a.grid(alpha=.2)
        a.legend(fontsize=8)
    fig.suptitle("Viscosity 0.01: matched disturbance pilot; physical stability not established")
    fig.savefig(ROOT / "comparison.png", dpi=140)
    plt.close(fig)
    last = rows[-1]
    direction = "No evolved comparison yet" if len(rows) == 1 else ("smaller than at the start" if A[-1] < 1 else "larger than at the start")
    summary = {"status": r["status"], "saved_time": last["t"], "saved_fields": len(rows),
        "initial_D": float(D[0]), "final_saved_D": float(D[-1]), "amplification": float(A[-1]),
        "direction": direction, "low_viscosity_amplification_at_same_time": float(low_A[-1]),
        "endpoint_finite_window_log_rate": float(np.log(A[-1]) / t[-1]) if t[-1] else None,
        "failed_spatial_screens": sum(not x["resolved_screen"] for x in rows),
        "max_physical_norm_relative_error": max(x["norm_relative_error"] for x in rows),
        "max_difference_energy_identity_relative_error": max(x["difference_energy_identity_relative_error"] for x in rows),
        "interpretation": "One finite-amplitude seed on an underresolved grid; no asymptotic Lyapunov, continuum stability or pond claim."}
    save(ROOT / "summary.json", summary)
    body = f"""<!doctype html><html><head><meta charset="utf-8"><title>Viscosity 0.01 disturbance comparison</title>
    <style>body{{font:17px/1.6 system-ui;color:#193447;max-width:1150px;margin:35px auto;padding:0 22px;background:#f6f9fc}}section{{background:white;border:1px solid #ccd9e3;padding:22px;margin:18px 0;border-radius:10px}}img{{width:100%}}a{{color:#08778a}}</style></head><body>
    <h1>Does a small disturbance shrink at viscosity 0.01?</h1><section><p><b>{html.escape(r['status'])}: saved through {last['t']:.2f} of 0.40.</b> {len(rows)} paired saved times are available.</p>
    <p>A 0.1% low-mode disturbance with seed 202 is compared with the existing undisturbed viscosity 0.01 flow. Grid 128, the initial flow, mean, domain, zero force, Heun method and timestep match. The original 48-job queue is unchanged; its undisturbed run is reused.</p>
    <p>The saved difference is <b>{direction}</b>: D={D[-1]:.8g}, giving <b>{A[-1]:.6g} times its starting size</b>. A value below one means this difference has shrunk over the recorded interval. A value above one means it has grown. The prior matched viscosity 0.001 case gives {low_A[-1]:.6g} at the same time.</p></section>
    <section><img src="comparison.png"><p>D is the L2 velocity difference divided by the initial reference-flow L2 norm. Transfer minus viscous removal gives the instantaneous rate of difference energy. This is different from measuring peak spin in either flow alone.</p></section>
    <section><h2>Limits of this pilot</h2><p>The raw strain is nonsmooth across periodic joins and initial maxima are outside the central tube. {summary['failed_spatial_screens']} of {len(rows)} disturbed-field spatial screens fail. The same grid and timestep remove a comparison mismatch, but do not establish resolution or time convergence. One 0.1% seed is a finite-amplitude, finite-time test; it does not establish infinitesimal or asymptotic stability, an asymptotic Lyapunov exponent, or behavior in a pond.</p>
    <p>The original sources are frozen by hashes. Physical-space and Fourier norms are cross-checked, and the difference-energy transfer identity is checked independently. Mean, divergence, finite fields, speed CFL and energy budgets have separate stops. Failed checkpoints and logs are preserved.</p>
    <p><a href="protocol.json">Protocol</a> · <a href="result.json">Saved measurements</a> · <a href="summary.json">Current summary</a> · <a href="initial-audit.json">Initial checks</a></p></section></body></html>"""
    target = ROOT / "report.html"
    temp = target.with_suffix(".tmp")
    temp.write_text(body, encoding="utf-8")
    os.replace(temp, target)


if __name__ == "__main__":
    build()
