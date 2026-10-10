"""Compare identical held-out material labels and frozen models across grids."""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(BASE))
import numpy as np
from features import dataset, future_target, index_at, wrap
from simulate import sha, save_json
from export_models import predict_model
from evaluate import summarize, paired


def file_for(folder, seed, n, dt=.005):
    return folder / f"s{seed}_n{n}_nu0.02_dt{dt}.npz"


def checked(path):
    meta = json.loads(path.with_suffix(".json").read_text())
    assert sha(path) == meta["sha256"]
    expected = {"simulate.py": sha(BASE/"simulate.py"),
                "numerics.py": sha(BASE.parent/"numerics.py"),
                "protocol.json": sha(BASE/"protocol.json")}
    # Windows executed the parent solver with CRLF; the published parent has LF.
    # Admit only byte hashes of these newline-equivalent forms of current source.
    sources = {"simulate.py": BASE/"simulate.py", "numerics.py": BASE.parent/"numerics.py",
               "protocol.json": BASE/"protocol.json"}
    assert set(meta["sources"]) == set(expected)
    for name, source_path in sources.items():
        assert source_hash_matches(source_path.read_bytes(), meta["sources"][name])
    with np.load(path, allow_pickle=False) as archive:
        data = dict(archive)
    assert all(np.isfinite(v).all() for v in data.values())
    return data, meta


def source_hash_matches(raw, recorded):
    lf = raw.replace(b"\r\n", b"\n")
    return recorded in {hashlib.sha256(x).hexdigest() for x in (raw, lf, lf.replace(b"\n", b"\r\n"))}


def material_comparison(a, b, p, cutoff):
    """Same initial material labels and observation times, never spatial nearest fit."""
    assert np.array_equal(a["time"], b["time"]), "Unmatched observation times"
    assert np.array_equal(a["positions"][0], b["positions"][0]), "Unmatched material labels"
    ya, yb, volume = [], [], []
    for t in p["anchors"]:
        i = index_at(a["time"], t+p["gap"])
        j = index_at(a["time"], t+p["gap"]+p["primary_horizon"])
        ya.extend(future_target(a["tangent"][i], a["tangent"][j], p["primary_horizon"]))
        yb.extend(future_target(b["tangent"][i], b["tangent"][j], p["primary_horizon"]))
        volume.extend(abs(np.linalg.det(b["tangent"][j])/np.linalg.det(b["tangent"][i])-1))
    ya, yb = np.array(ya), np.array(yb)
    ea, eb = ya >= cutoff, yb >= cutoff
    union = int(np.sum(ea | eb))
    rms = float(np.sqrt(np.mean((ya-yb)**2)))
    path_diff = np.linalg.norm(wrap(a["positions"]-b["positions"]), axis=-1)
    return {"target_RMSE": rms, "relative_target_RMSE": rms/float(np.sqrt(np.mean(yb**2))),
            "target_event_Jaccard": float(np.sum(ea & eb)/union) if union else 1.,
            "target_event_flip_fraction": float(np.mean(ea != eb)),
            "coarse_event_count": int(ea.sum()), "fine_event_count": int(eb.sum()),
            "path_RMS": float(np.sqrt(np.mean(path_diff**2))), "path_max": float(path_diff.max()),
            "fine_target_volume_error_p95": float(np.quantile(volume, .95)),
            "fine_target_volume_error_max": float(np.max(volume))}


def mean_rows(rows):
    return {k: float(np.mean([r[k] for r in rows])) for k in rows[0] if k != "seed"}


def matched_snapshot(folder):
    """Compare all three original grids on exactly the same saved time interval."""
    parsed, hashes = {}, {}
    for n in (64, 128, 256):
        path = folder/"runs"/f"baseline-n{n}-base"/"result.json"
        raw = path.read_bytes()
        parsed[n] = json.loads(raw)
        hashes[str(n)] = hashlib.sha256(raw).hexdigest()
    rows = {n: {round(x["t"], 8): x for x in r["series"]} for n, r in parsed.items()}
    times = sorted(set(rows[64]) & set(rows[128]) & set(rows[256]))
    comparisons = {}
    for lo, hi in ((64, 128), (128, 256)):
        comparisons[f"{lo}-{hi}"] = {}
        for key in ("Wmax", "I", "energy", "enstrophy"):
            x = np.array([rows[lo][t][key] for t in times])
            y = np.array([rows[hi][t][key] for t in times])
            comparisons[f"{lo}-{hi}"][key] = {"relative_curve_l2": float(np.linalg.norm(x-y)/np.linalg.norm(y)),
                                               "absolute_curve_RMS": float(np.sqrt(np.mean((x-y)**2)))}
    return {"through": times[-1], "observation_count": len(times), "result_json_sha256": hashes,
            "saved_times": {str(n): r["series"][-1]["t"] for n, r in parsed.items()},
            "comparisons_same_interval": comparisons,
            "limits": "Partial original matrix. Nonsmooth periodic joins and off-central initial maxima retained. These global vorticity comparisons do not validate marker FTLE or the smooth-flow prediction pilot."}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matched-root", type=Path)
    args = ap.parse_args()
    p = json.loads((BASE/"protocol.json").read_text())
    followup = json.loads((HERE/"protocol.json").read_text())
    execution = json.loads((HERE/"execution.json").read_text())
    assert execution["complete"] and len(execution["records"]) == followup["new_simulations"]
    assert execution["protocol_sha256"] == sha(HERE/"protocol.json")
    assert execution["runner_sha256"] == sha(HERE/"run_refinement.py")
    models_path = BASE/"results"/"primary_models.json"
    models_sha = sha(models_path)
    models = json.loads(models_path.read_text())
    cutoff = models["event_target_cutoff"]
    data_root = BASE/"data" if (BASE/"data").is_dir() else BASE/"recorded-data"
    new_root = HERE/"recorded-data"
    results = {"built_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "protocol_sha256": sha(HERE/"protocol.json"), "frozen_models_sha256": models_sha,
               "original_protocol_sha256": sha(BASE/"protocol.json"), "train_event_cutoff": cutoff,
               "independent_test_groups": 8, "new_simulations": 12, "models_refitted": False,
               "grid_models": {}, "spatial": {}, "temporal": {}, "files": []}
    results["analysis_sources_sha256"] = {str(path.relative_to(BASE)): sha(path) for path in
        (Path(__file__), BASE/"features.py", BASE/"export_models.py", BASE/"evaluate.py")}
    predictions = {}
    for n in (26, 38, 56):
        paths = [file_for(new_root if n == 56 else data_root, seed, n) for seed in p["test_seeds"]]
        for path in paths:
            raw_data, meta = checked(path)
            trace = np.trace(raw_data["gradient"], axis1=-2, axis2=-1)
            results["files"].append({"path": str(path.relative_to(BASE)), "sha256": meta["sha256"],
                "seed": meta["seed"], "config": meta["config"], "max_CFL": meta["max_cfl"],
                "max_energy_residual": max(abs(x["energy_budget_residual"]) for x in meta["diagnostics"]),
                "max_high_band_energy": max(x["high_band_energy_fraction"] for x in meta["diagnostics"]),
                "max_spectral_field_divergence": max(x["divergence_max"] for x in meta["diagnostics"]),
                "recorded_source_sha256": meta["sources"],
                "interpolated_divergence_RMS": float(np.sqrt(np.mean(trace**2))),
                "interpolated_divergence_max": float(abs(trace).max())})
        d = dataset(paths, p)
        e = (d["target"] >= cutoff).astype(int)
        grid = {}
        predictions[f"n{n}_target"] = d["target"]
        predictions[f"n{n}_position"] = d["position"]
        for key in ("group", "particle", "anchor"):
            if key in predictions:
                assert np.array_equal(predictions[key], d[key])
            predictions[key] = d[key]
        for name, X in (("A_full", d["present"]), ("B_history", np.column_stack([d["present"], d["history"]]))):
            model = models["models"][name]
            pred = predict_model(model["continuous"], X)
            prob = predict_model(model["event"], X)
            grid[name] = summarize(d["target"], e, pred, prob, d["group"], model["decision_threshold"])
            predictions[f"n{n}_{name}_prediction"] = pred
            predictions[f"n{n}_{name}_probability"] = prob
        grid["paired"] = paired(grid["A_full"], grid["B_history"], p)
        results["grid_models"][str(n)] = grid
    # Reused coarse-grid predictions must reproduce the published pilot exactly.
    with np.load(BASE/"results"/"predictions.npz", allow_pickle=False) as original:
        for name in ("A_full", "B_history"):
            for kind in ("prediction", "probability"):
                assert np.max(abs(original[f"{name}_{kind}"]-predictions[f"n26_{name}_{kind}"])) < 1e-10
    for lo, hi in ((26, 38), (38, 56)):
        rows = []
        for seed in p["test_seeds"]:
            a, ma = checked(file_for(data_root, seed, lo))
            b, mb = checked(file_for(new_root if hi == 56 else data_root, seed, hi))
            row = {"seed": seed, **material_comparison(a, b, p, cutoff)}
            ea = np.array([x["energy"] for x in ma["diagnostics"]])
            eb = np.array([x["energy"] for x in mb["diagnostics"]])
            row["energy_curve_relative_l2"] = float(np.linalg.norm(ea-eb)/np.linalg.norm(eb))
            mask = predictions["group"] == seed
            for name in ("A_full", "B_history"):
                diff = predictions[f"n{lo}_{name}_prediction"][mask]-predictions[f"n{hi}_{name}_prediction"][mask]
                row[name+"_prediction_drift_RMSE"] = float(np.sqrt(np.mean(diff**2)))
            rows.append(row)
        results["spatial"][f"{lo}-{hi}"] = {"per_run": rows, "mean": mean_rows(rows)}
    for n in followup["half_step_grids"]:
        rows = []
        for seed in followup["half_step_seeds"]:
            a, _ = checked(file_for(new_root if n == 56 else data_root, seed, n))
            path = file_for(new_root, seed, n, .0025)
            b, meta = checked(path)
            rows.append({"seed": seed, **material_comparison(a, b, p, cutoff)})
            results["files"].append({"path": str(path.relative_to(BASE)), "sha256": meta["sha256"],
                                     "seed": seed, "config": meta["config"]})
        results["temporal"][str(n)] = {"per_run": rows, "mean": mean_rows(rows)}
    coarse = results["spatial"]["26-38"]["per_run"]
    fine = results["spatial"]["38-56"]["per_run"]
    screens = followup["screens"]
    verdict = {
        "all_seed_grid_difference_contracts": all(b["target_RMSE"] < a["target_RMSE"] for a, b in zip(coarse, fine)),
        "finest_relative_target_RMSE_below_1pct_all_seeds": all(r["relative_target_RMSE"] <= screens["finest_pair_target_relative_RMSE_max"] for r in fine),
        "finest_event_Jaccard_at_least_95pct_all_seeds": all(r["target_event_Jaccard"] >= screens["finest_pair_event_Jaccard_min"] for r in fine),
        "time_error_below_tenth_grid_error_checked_seeds": all(t["target_RMSE"] < .1*next(r["target_RMSE"] for r in fine if r["seed"] == t["seed"])
            for v in results["temporal"].values() for t in v["per_run"])}
    results["screens"] = verdict
    fine_gain = abs(results["grid_models"]["56"]["paired"]["RMSE"]["mean_B_minus_A"])
    results["finest_target_change_to_AB_error_difference_ratio"] = results["spatial"]["38-56"]["mean"]["target_RMSE"]/max(fine_gain, 1e-30)
    results["numerical_stability_screens_pass"] = all(verdict.values())
    results["continuum_convergence_established"] = False
    results["interpretation"] = "These comparisons assess numerical stability on the declared finite window, not a continuum proof. A/B scores use frozen coarse-trained fits; this is not a new independent replication or a fine-grid-trained optimal model comparison. No Eulerian first-event location or exact onset target is tested."
    if args.matched_root:
        results["original_matched_stretch_snapshot"] = matched_snapshot(args.matched_root)
    assert models_sha == sha(models_path)
    save_json(HERE/"results.json", results)
    np.savez_compressed(HERE/"predictions.npz", **predictions)
    with (HERE/"per-run-comparisons.csv").open("w", newline="") as stream:
        fields = ["comparison"] + list(fine[0])
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for name, val in results["spatial"].items():
            writer.writerows({"comparison": name, **r} for r in val["per_run"])
    make_report(results)
    print(json.dumps({"screens": verdict, "spatial": {k: v["mean"] for k, v in results["spatial"].items()},
                      "models": {n: v["paired"]["RMSE"] for n, v in results["grid_models"].items()}}, indent=2))


def make_report(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    grids = [26, 38, 56]
    for seed in range(8):
        axes[0, 0].plot([38, 56], [r["spatial"][p]["per_run"][seed]["target_RMSE"] for p in ("26-38", "38-56")], "o-", alpha=.7)
        axes[0, 1].plot([38, 56], [r["spatial"][p]["per_run"][seed]["target_event_Jaccard"] for p in ("26-38", "38-56")], "o-", alpha=.7)
    axes[0, 0].set(title="Target disagreement between successive grids", ylabel="FTLE RMSE", xlabel="Finer grid size")
    axes[0, 1].axhline(.95, color="black", linestyle="--", linewidth=1)
    axes[0, 1].set(title="Strong-event set agreement", ylabel="Jaccard overlap", xlabel="Finer grid size", ylim=(0, 1.03))
    for ax in axes[0]:
        ax.set_xticks([38, 56], ["26 → 38", "38 → 56"])
        ax.set_xlabel("Successive grid comparison")
    for ax, metric, title in ((axes[1, 0], "RMSE", "B − A prediction error (negative is better)"),
                             (axes[1, 1], "AUPRC_AP", "B − A event average precision (positive is better)")):
        vals = [r["grid_models"][str(n)]["paired"][metric] for n in grids]
        y = np.array([v["mean_B_minus_A"] for v in vals])
        ci = np.array([v["ci95"] for v in vals])
        ax.errorbar(grids, y, yerr=np.stack([y-ci[:, 0], ci[:, 1]-y]), fmt="o-", capsize=5)
        ax.axhline(0, color="black", linestyle="--", linewidth=1)
        ax.set(title=title, xlabel="Grid size; same 8 test flows")
        ax.set_xticks(grids)
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Marker history: three-grid numerical follow-up · frozen models", fontsize=14)
    fig.savefig(HERE/"convergence.png", dpi=160)
    fig.savefig(HERE/"convergence.svg")
    plt.close(fig)
    lines = ["# Marker-history spatial-resolution follow-up", "",
        "**Spatial convergence remains unestablished.** This bounded follow-up reuses the original 26³/38³ data and adds 56³ for the same eight held-out initial conditions. It also halves the timestep on grids 38³ and 56³ for two prechosen seeds. All 12 new runs retain the same smooth initial fields, viscosity, duration, marker labels, observation times and unchanged Navier–Stokes solver. No predictor was refitted.", "",
        "These are refinements of the smooth-flow prediction pilot, separate from the original sharp strain/tube matrix. The original test seeds have already been examined; this follow-up is exploratory numerical validation, not a fresh independent statistical confirmation.", "",
        "A uses current velocity, the full present local velocity gradient and current neighbor geometry. B uses exactly those inputs plus causal marker-history features. Neither predictor observes future positions or persistent marker identity.", "",
        "| Grid | A RMSE | B RMSE | B − A RMSE (95% paired run interval) | A event AP | B event AP |", "|---|---:|---:|---|---:|---:|"]
    for n, v in r["grid_models"].items():
        a, b, d = v["A_full"]["macro"], v["B_history"]["macro"], v["paired"]["RMSE"]
        lines.append(f"| {n}³ | {a['RMSE']:.6f} | {b['RMSE']:.6f} | {d['mean_B_minus_A']:.6f} [{d['ci95'][0]:.6f}, {d['ci95'][1]:.6f}] | {a['AUPRC_AP']:.6f} | {b['AUPRC_AP']:.6f} |")
    lines += ["", "![Convergence checks](convergence.png)", "", "Each colored line in the upper panels is one held-out flow. Lower-panel error bars are 95% intervals from 2,000 paired bootstraps of the same eight independent seed groups, not bootstraps of individual markers.", "", "| Successive grids | Mean target RMSE | Mean relative target RMSE | Mean event-set Jaccard |", "|---|---:|---:|---:|"]
    for pair, v in r["spatial"].items():
        x = v["mean"]
        lines.append(f"| {pair} | {x['target_RMSE']:.6f} | {100*x['relative_target_RMSE']:.3f}% | {x['target_event_Jaccard']:.4f} |")
    lines += ["", "Events use the original training-derived threshold at every grid, never a new per-grid quantile. Jaccard measures overlap of the sets of material markers classified as strong future deformation; the future Eulerian location and exact first onset are not predicted here.", "", "## Predeclared numerical screens", ""]
    labels = ["Target differences decrease for all eight flows", "Finest-pair target change is below 1% for every flow",
              "Finest-pair event-set overlap is at least 95% for every flow", "Timestep change is below one tenth of grid change in the checked flows"]
    lines += [f"- {label}: **{'pass' if v else 'fail'}**." for label, v in zip(labels, r["screens"].values())]
    lines += ["", "These tolerances are operational screens, not a proof. Decreasing differences alone do not establish an asymptotic regime. No Richardson order or extrapolated continuum answer is claimed. Three-grid checks and their assumptions follow [NASA's spatial convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).", "", "## Timestep checks", ""]
    for n, v in r["temporal"].items():
        lines.append(f"- Grid {n}³, two prechosen seeds: mean target change on halving timestep = {v['mean']['target_RMSE']:.6g}.")
    lines += ["", "## Interpolation diagnostic", "", "The underlying spectral fluid is divergence-free to numerical roundoff, but the velocity interpolated between cells need not be. The gradient saved for each marker is the exact derivative of that interpolant. The following are means over the eight runs; the divergence RMS uses all saved times and markers.", "", "| Grid | Maximum spectral-field divergence | Interpolated divergence RMS |", "|---|---:|---:|"]
    for n in grids:
        rows = [f for f in r["files"] if f["config"]["n"] == n and f["config"]["dt"] == .005]
        lines.append(f"| {n}³ | {np.mean([f['max_spectral_field_divergence'] for f in rows]):.3g} | {np.mean([f['interpolated_divergence_RMS'] for f in rows]):.6g} |")
    lines += ["", "This diagnoses an interpolation contribution, not a complete decomposition of all target error. Changing the interpolation would require a separately documented measurement/integration variant and fresh validation.", ""]
    lines += ["", "## What can be concluded", "", "The original small continuous-error benefit was out of sample with respect to initial-condition seeds. Its physical interpretation still depends on numerical accuracy. Marker relabeling invariance demonstrates that the code uses geometry rather than persistent identity; it does not by itself demonstrate predictive usefulness, causal influence, or a necessary role for the direction of time. The original reversed-history control retained a similar benefit.", "", "Trilinear interpolation and its piecewise spatial derivative contribute to the grid dependence of the numerical FTLE. Frozen coarse-trained models may also shift under refined observations. Stronger evidence would require a sufficiently stable velocity/gradient interpolation and trajectory/tangent calculation, additional resolved grid checks, then independently held-out model comparisons with numerical uncertainty small relative to the claimed benefit.", ""]
    lines += [f"The finest-pair target change is {r['finest_target_change_to_AB_error_difference_ratio']:.1f} times the absolute A/B RMSE difference at grid 56. This is a scale comparison, not a statistical error bar on the paired A/B difference; numerical errors can be correlated across predictors.", ""]
    lines += ["Reversing a past history preserves its collection of geometries and transforms several features only by sign. Because the reversed-history predictor was trained on those transformed features, it can compensate. This is a limited control: its similar score does not prove that chronology is physically irrelevant. The study measures stretching/deformation, not mechanical power, and does not establish that the location of power matters more than its organization.", ""]
    if "original_matched_stretch_snapshot" in r:
        s = r["original_matched_stretch_snapshot"]
        lines += ["## Separate original matched-stretch matrix", "", f"Snapshot comparison uses exactly {s['observation_count']} common samples through t={s['through']:.2f} for grids 64³, 128³ and 256³. It does not compare different time windows.", ""]
        for pair, v in s["comparisons_same_interval"].items():
            lines.append(f"- {pair}: maximum-vorticity curve relative difference {100*v['Wmax']['relative_curve_l2']:.3f}%; energy {100*v['energy']['relative_curve_l2']:.3f}%; enstrophy {100*v['enstrophy']['relative_curve_l2']:.3f}%.")
        lines += ["", s["limits"], ""]
    lines += ["## Reproduce and audit", "", "From this directory, using the original study environment:", "", "```text", "python -B run_refinement.py --pilot", "python -B run_refinement.py", "python -B analyze.py", "python -B -m pytest test_analysis.py -q", "```", "", "The optional `--matched-root` argument takes the local original matched-stretch-study directory and reads a snapshot only. Numerical hashes and per-run metrics are in `results.json`, raw new trajectories with provenance in `recorded-data/`, predictions in `predictions.npz`, and the bounded design in `protocol.json`. Reused raw data are in the parent study's `recorded-data/` (local execution may use `data/`). No original source or completed result is overwritten.", ""]
    lines += ["To audit the published data directly, run `python -B analyze.py`; new simulations are unnecessary. Analysis verifies every data SHA-256 and source hashes, permitting only LF/CRLF variants of otherwise identical current source. The original simulation cache remains stricter and requires exact executed source bytes, including line endings; its recorded Windows solver SHA-256 is `16d11a080c1559d428b33f803c1f83cfabbc3336e0bab9256989b9a2cef40d53`, while the published LF parent is `22b169971aaf2efd56d47f685495a575e1622edd56c199b07a785fc97e952f46`. Their normalized text was verified identical. Run regeneration in a separate fresh study checkout when its source bytes differ from the archived cache.", ""]
    (HERE/"README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
