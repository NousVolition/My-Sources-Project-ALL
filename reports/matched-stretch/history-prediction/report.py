"""Generate a readable report, data tables and static scientific figures."""
import base64
import csv
import html
import json
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parent/"mpl-cache"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from matplotlib.ticker import MaxNLocator
from simulate import ROOT,save_json

BLUE="#2463a6"; ORANGE="#cc6523"; GREY="#6b7280"
LABELS={"A_full":"A: current full gradient","B_history":"B: + marker history","A_weak":"Current weak baseline",
        "B_weak":"Weak baseline + history","A_matched_current":"Equal-size current-only terms","history_only":"History alone",
        "B_past_only":"+ past affine deformation only","B_no_past":"+ history excluding affine deformation",
        "B_no_exchange":"+ history excluding exchanges","B_reversed":"+ reversed past history","B_irrelevant":"+ distant-marker history",
        "B_permuted_labels":"Shuffled training targets","B_shuffled_91":"+ shuffled history 1","B_shuffled_92":"+ shuffled history 2",
        "B_shuffled_93":"+ shuffled history 3","tree_A_full":"Tree A: current","tree_B_history":"Tree B: + history"}


def table(headers,rows):
    return "<div class='scroll'><table><thead><tr>"+"".join(f"<th>{html.escape(str(x))}</th>" for x in headers)+"</tr></thead><tbody>"+"".join("<tr>"+"".join(f"<td>{html.escape(str(x))}</td>" for x in row)+"</tr>" for row in rows)+"</tbody></table></div>"


def figsave(fig,name):
    folder=ROOT/"results"/"figures";folder.mkdir(exist_ok=True)
    fig.savefig(folder/(name+".png"),dpi=160,bbox_inches="tight")
    fig.savefig(folder/(name+".svg"),bbox_inches="tight")
    svg=folder/(name+".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n",encoding="utf-8")
    plt.close(fig)


def image(name,caption):
    data=base64.b64encode((ROOT/"results"/"figures"/(name+".png")).read_bytes()).decode()
    return f"<figure><img src='data:image/png;base64,{data}' alt='{html.escape(caption)}'><figcaption>{html.escape(caption)}</figcaption></figure>"


def main():
    out=ROOT/"results"
    primary=json.loads((out/"primary.json").read_text()); sensitivities=json.loads((out/"sensitivity.json").read_text())
    verify=json.loads((out/"verification.json").read_text()); old=json.loads((out/"existing_data_audit.json").read_text())
    protocol=json.loads((ROOT/"protocol.json").read_text())
    with np.load(out/"predictions.npz",allow_pickle=False) as d: pred=dict(d)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                         "axes.labelcolor":"#283342","axes.titleweight":"bold","figure.facecolor":"white"})
    models=primary["models"]; A=models["A_full"];B=models["B_history"]
    delta=primary["paired"]["B_history_vs_A_full"]
    pct=-100*delta["RMSE"]["mean_B_minus_A"]/A["macro"]["RMSE"]
    # Primary and controls, paired by independent test run.
    chosen=["B_history","B_past_only","B_no_past","B_no_exchange","B_reversed","B_irrelevant","B_shuffled_91","B_shuffled_92","B_shuffled_93"]
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for ax,metric,title in zip(axs,["RMSE","AUPRC_AP"],["Continuous error: left is better","Event ranking: right is better"]):
        for i,name in enumerate(chosen):
            r=primary["paired"][name+"_vs_A_full"][metric]; val=r["mean_B_minus_A"];lo,hi=r["ci95"]
            ax.errorbar(val,i,xerr=[[val-lo],[hi-val]],fmt="o",color=ORANGE if name=="B_history" else BLUE,capsize=3)
        ax.axvline(0,color=GREY,lw=1);ax.set_yticks(range(len(chosen)),[LABELS[n] for n in chosen] if ax==axs[0] else [])
        ax.invert_yaxis();ax.set_title(title);ax.set_xlabel("Difference from A (95% paired run bootstrap)");ax.grid(axis="x",alpha=.15)
    figsave(fig,"controls")
    # Pooled PR and reliability; macro numbers are in the table, not these curves.
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout="constrained")
    for name,color in [("A_full",BLUE),("B_history",ORANGE)]:
        precision,recall,_=precision_recall_curve(pred["event"],pred[name+"_probability"])
        axs[0].plot(recall,precision,color=color,label=LABELS[name])
        bins=models[name]["calibration"]
        axs[1].plot([b["probability"] for b in bins],[b["observed"] for b in bins],"o-",color=color,label=LABELS[name])
    axs[0].axhline(pred["event"].mean(),color=GREY,ls="--",label="Test prevalence")
    axs[0].set(xlabel="Recall",ylabel="Precision",title="Held-out precision–recall (pooled)",xlim=(0,1),ylim=(0,1.02))
    axs[1].plot([0,1],[0,1],"--",color=GREY);axs[1].set(xlabel="Predicted probability",ylabel="Observed event fraction",title="Calibration (10 fixed bins)",xlim=(0,1),ylim=(0,1))
    for ax in axs:ax.legend(fontsize=8);ax.grid(alpha=.15)
    figsave(fig,"events")
    fig,axs=plt.subplots(1,2,figsize=(12,4),layout="constrained")
    for ax,metric,title in zip(axs,["RMSE","AUPRC_AP"],["Error by independent test seed","Event ranking by independent test seed"]):
        for name,color in [("A_full",BLUE),("B_history",ORANGE)]:
            by=models[name]["per_run"]; keys=sorted(by)
            ax.plot(range(len(keys)),[by[k][metric] for k in keys],"o-",label=LABELS[name],color=color)
        ax.set_xticks(range(len(keys)),keys,rotation=45);ax.set_title(title);ax.set_ylabel(metric);ax.grid(alpha=.15);ax.legend(fontsize=8)
    figsave(fig,"per_run")
    fig,axs=plt.subplots(1,2,figsize=(12,7),layout="constrained")
    for ax,metric,title in zip(axs,["RMSE","AUPRC_AP"],["Change in continuous error","Change in event ranking"]):
        for i,(name,res) in enumerate(sensitivities.items()):
            r=res["paired"]["B_history_vs_A_full"][metric];v=r["mean_B_minus_A"];lo,hi=r["ci95"]
            ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt="o",color=BLUE,capsize=3)
        ax.axvline(0,color=GREY);ax.set_yticks(range(len(sensitivities)),list(sensitivities) if ax==axs[0] else []);ax.invert_yaxis()
        ax.set_title(title);ax.set_xlabel("B − A, paired 95% interval");ax.grid(axis="x",alpha=.15)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        if metric=="RMSE": ax.ticklabel_format(axis="x",style="sci",scilimits=(-3,-3))
    figsave(fig,"sensitivity")
    # Prespecified example: first test seed, middle anchor. Do not choose best case.
    use=(pred["group"]==1024)&np.isclose(pred["anchor"],.5)
    pos=(pred["position"][use]+3)%6-3
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout="constrained",sharex=True,sharey=True)
    vals=[pred["event"][use],pred["A_full_probability"][use],pred["B_history_probability"][use]]
    for ax,z,title in zip(axs,vals,["Observed future event","A: current information","B: current + history"]):
        sc=ax.scatter(pos[:,0],pos[:,1],c=z,cmap="viridis",vmin=0,vmax=1,s=12)
        ax.set(title=title,xlabel="Marker x at prediction time",aspect="equal",xlim=(-3,3),ylim=(-3,3))
    axs[0].set_ylabel("Marker y at prediction time");fig.colorbar(sc,ax=axs,label="Event / predicted probability",shrink=.75)
    figsave(fig,"location")
    rows=[]
    for name,m in models.items():
        a=m["macro"]; rows.append([LABELS[name],m["input_dimensions"]]+[f"{a[k]:.5f}" for k in ["RMSE","MAE","R2","AUPRC_AP","Brier","ECE10","precision","recall"]])
    headers=["Model","Inputs","RMSE ↓","MAE ↓","R² ↑","AP ↑","Brier ↓","ECE ↓","Precision","Recall"]
    with (out/"metrics.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(rows)
    srows=[]
    for name,res in sensitivities.items():
        da=res["paired"]["B_history_vs_A_full"]
        srows.append([name,f"{res['models']['A_full']['macro']['RMSE']:.5f}",f"{res['models']['B_history']['macro']['RMSE']:.5f}",
                      f"{da['RMSE']['mean_B_minus_A']:.6f}",str([round(x,6) for x in da['RMSE']['ci95']]),
                      f"{da['AUPRC_AP']['mean_B_minus_A']:.5f}","frozen primary fits" if res.get("frozen_primary_models") else "train/validation refit"])
    with (out/"sensitivity.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["Variant","A RMSE","B RMSE","Delta RMSE","95% CI","Delta AP","Fit policy"]);w.writerows(srows)
    max_energy=max(r["max_abs_energy_residual"] for r in verify["files"])
    max_div=max(r["max_divergence"] for r in verify["files"])
    grid=float(np.mean([r["target_RMSE"] for r in verify["refinement"]["grid"]]))
    time=float(np.mean([r["target_RMSE"] for r in verify["refinement"]["time"]]))
    det=max(r["refined_target_volume_error_max"] for r in verify["refinement"]["grid"])
    maxvolume=max(r["max_volume_error_from_time_zero"] for r in verify["files"])
    text=f"""<p class='eyebrow'>REPRODUCIBLE RESEARCH PILOT · 3D NAVIER–STOKES</p>
<h1>Do marker histories forecast deformation?</h1>
<p class='lead'>A small held-out improvement in continuous deformation prediction. Event results depend on the threshold, and the finer-grid result remains uncertain.</p>
<div class='cards'><div><strong>{pct:.2f}%</strong><span>lower mean per-run RMSE with history</span></div><div><strong>{delta['AUPRC_AP']['mean_B_minus_A']:+.4f}</strong><span>change in mean event average precision</span></div><div><strong>8</strong><span>independent test initial conditions</span></div></div>
<h2>What the result supports</h2>
<p>The primary regularized linear comparison gives RMSE {A['macro']['RMSE']:.5f} → {B['macro']['RMSE']:.5f}. The paired difference is {delta['RMSE']['mean_B_minus_A']:.6f}, with a 95% run-bootstrap interval [{delta['RMSE']['ci95'][0]:.6f}, {delta['RMSE']['ci95'][1]:.6f}]. Seven of eight test seeds improve. This is a small out-of-sample predictive association beyond a full present local gradient, within this smooth-flow ensemble.</p>
<p>Event average precision changes from {A['macro']['AUPRC_AP']:.5f} to {B['macro']['AUPRC_AP']:.5f}; the difference interval [{delta['AUPRC_AP']['ci95'][0]:.5f}, {delta['AUPRC_AP']['ci95'][1]:.5f}] spans zero. MAE improvement also has an interval spanning zero. The same-capacity tree models show a small continuous-error reduction and no clear event-ranking gain. These results do not establish a practically reliable warning for the next strong deformation.</p>
<p>The prechosen secondary top-5% event test is positive: AP 0.80112 → 0.83454, paired difference +0.03342 with interval [0.00742, 0.07276]. It needs confirmation: it is one of many exploratory comparisons and the largest per-run AP change comes from a seed with only two events. At the finer 38³ grid, the continuous-error difference is −0.000471 with interval [−0.000974, +0.000048], spanning zero. A robust resolution-independent benefit has not been established.</p>
<p>Reversed past histories retain a similar gain. The data therefore do not isolate a direction-sensitive reorganization mechanism or identify one special arrangement. Passive marker labels do not enter any predictor. Consistently permuting particle storage order preserves real-run features to {max(r['max_feature_error_after_ID_permutation'] for r in verify['identity_invariance']):.2g}. This checks freedom from identity artifacts; it is not evidence that marker arrangements cause the flow.</p>
<h2>Data, equations and the locked split</h2>
<p>Source: <a href='https://github.com/NousVolition/My-Sources-Project-ALL/tree/17c14e89940d4c96ec9df315b1a689430d8c5ddf/reports/matched-stretch'>NousVolition/My-Sources-Project-ALL, commit 17c14e89940d4c96ec9df315b1a689430d8c5ddf</a>. The parent <code>numerics.Flow</code> solver is reused unchanged: unforced incompressible 3D Navier–Stokes on a periodic box of side 6, Fourier projection, spectral rotational nonlinearity and Heun stepping. No Hénon map, particle force, feedback, molecular model or unrelated surrogate dynamics is used.</p>
<p>New initial conditions are explicitly different from the original sharp strain/tube start: independent random low-mode Fourier coefficients and phases, projected divergence-free, zero mean and RMS speed 1. Thirty-two seed groups supply 16 training, 8 validation and 8 test runs. The main grid is 26³, viscosity 0.02, timestep 0.005, saved spacing 0.025 and duration 1.2. The nominal initial U<sub>rms</sub>L/ν is 300 (the viscosity probes are 600 and 150); this does not establish fully developed turbulence. Each run has 512 passive markers. All grid, timestep and viscosity variants of one seed stay in the same group.</p>
<p>There are 64 simulations in total: 32 main and 32 dependent test-seed sensitivity runs. The predeclared split is seeds 1000–1015 / 1016–1023 / 1024–1031. Five fixed prediction times (0.3, 0.4, 0.5, 0.6, 0.7) give {primary['n_train']:,} training, {primary['n_validation']:,} validation and {primary['n_test']:,} test rows. Overlapping observations remain entirely inside their run. Uncertainty resamples eight runs, never these rows. No test result selects models, features, thresholds, horizons or examples.</p>
<h2>Exactly what is predicted</h2>
<p>For a marker at time t, the main target is delayed forward FTLE, <code>log σmax[M(t+0.025+0.2) M(t+0.025)⁻¹] / 0.2</code>. M obeys the variational equation along the numerical marker flow. Prediction inputs end at t. Target deformation starts at t+0.025, so the integration intervals do not overlap. The prior tangent history cancels algebraically; an automated test verifies this basis invariance. No future position, tangent matrix or velocity is supplied to a predictor.</p>
<p>The event is target ≥ {primary['train_event_cutoff']:.6f}, the training-only 90th percentile. Test prevalence is {100*A['pooled']['prevalence']:.2f}%; test labels are never re-quantiled. The target identifies deformation of parcels currently at each marker location in a fixed future window. It does not predict an exact onset time, a first-passage time, or a future Eulerian location. Reported FTLE values have inverse model-time units.</p>
<h2>Matched information and models</h2>
<p>A has 31 current inputs: 3 velocity components; all 9 gradient entries; 3 strain eigenvalues; strain and rotation norms; speed; four distance summaries; six entries of the normalized neighbor-direction covariance; strain alignment; and pair-angle mean and standard deviation. B contains those identical 31 inputs plus 23 history summaries over 0.2 time units, using 8 current nearest neighbors and their known past trajectories. History covers distance changes and paths, angular changes, neighborhood exchanges, alignment changes, fitted past affine deformation and shape changes. Full feature names and standardized coefficients are saved.</p>
<p>The full gradient is the exact derivative of the same piecewise trilinear velocity interpolant used by the markers, not a claim of exact continuum measurement. The primary regression is standardized ridge; events use standardized logistic regression. Scaling is fitted only on training runs. Validation selects ridge α from [0.1, 1, 10, 100, 1000] by mean per-run RMSE and logistic C from [0.01, 0.1, 1, 10] by mean per-run log loss. Validation also chooses the probability cutoff for precision/recall by F1; AP and calibration do not use that cutoff. There is no separate test-fitted calibration.</p>
<p>Equal-size current-only controls append 23 preselected quadratic current terms. Trees in A and B have identical depth 3, up to 7 leaves, 80 boosting steps and validation-selected L2 from [1, 10], with early stopping disabled. The impoverished baseline uses velocity, largest strain eigenvalue, strain/rotation norms and speed. Current position, time, seed, viscosity and marker ID are excluded from all predictor inputs.</p>
<h2>Held-out metrics</h2><p>Values below are unweighted means over the eight independent test runs. AP is non-interpolated average precision (a precision–recall area summary). ECE uses 10 fixed probability bins. Pooled values, per-run values, calibration bins, selected hyperparameters, coefficients and predictions are in the machine-readable results.</p>
{table(headers,rows)}
{image('controls','Paired differences from the full-current baseline. Intervals resample whole independent test runs, 2,000 replicates.')}
{image('events','Pooled curves are descriptive; the primary comparisons and uncertainty use independent-run macro metrics.')}
{image('per_run','All eight held-out runs, without selecting successful cases. Rare event counts make individual AP values variable.')}
<h2>Controls and interpretation</h2>
<p>Three history permutations independently shuffle history rows within each run and prediction time, preserving current inputs and the run/time distribution. Distant-marker histories retain a real but spatially irrelevant history. Both largely remove the small RMSE gain. The reversed-history control reverses only already-observed past snapshots: it still contains useful deformation information and is not guaranteed to have zero skill. Similar reversed performance prevents a claim that forward temporal ordering is necessary.</p>
<p>Training-label permutations are also within run/time. They destroy local geometry-to-target matching but intentionally retain run/time marginal structure, so they need not produce prevalence-level AP. Validation and test labels remain true. The large performance loss is a leakage/association diagnostic. Identity permutations consistently rename/reorder a marker throughout its trajectory and target; shuffling correspondence independently at each time would instead destroy trajectories and would not be a mere renaming.</p>
<p>The 95% intervals are paired percentile cluster bootstraps, conditional on the fixed fitted models. They do not include training-ensemble uncertainty. Unadjusted two-sided paired sign-flip p-values are also saved (primary RMSE {delta['RMSE']['two_sided_signflip_p_unadjusted']:.5f}). Many secondary comparisons are exploratory, without multiplicity correction. Eight test groups are a small inferential sample.</p>
<h2>Sensitivity: prechosen variations</h2>
{table(['Variant','A RMSE','B RMSE','B − A','95% paired interval','Δ AP','Fit policy'],srows)}
{image('sensitivity','Lookback 0.1/0.3; horizon 0.1/0.4; neighbors 4/16; markers 256; history observation stride 2; event quantiles 0.8/0.95. Grid 38³, half timestep and viscosities 0.01/0.04 use frozen main models on the same test seeds. These are one-factor checks, not a full factorial sweep.')}
<h2>Where is the forecast assigned?</h2>
{image('location','Fixed example: seed 1024, prediction time 0.5, target interval 0.525–0.725. Each point is a current marker location; z is collapsed in this projection. These panels are not forecasts of future marker positions.')}
<h2>Original marker-stress data reused</h2>
<p>All eight existing marker-stress trajectory archives were hashed and reprocessed through the geometry extractor: reference, coarse/medium tracking, fluid half-step, three 128-grid variants, and coarse saved times. Only their first 250 primary markers were used; companions and exact clones were excluded. These are one initial-condition family, not eight independent training/test runs. They therefore contribute a compatibility/sensitivity audit, not an out-of-sample claim.</p>
<p>The audit uses a separate finite-neighborhood affine logarithmic stretch target (not tangent FTLE), with lookback 0.04, gap 0.02, horizon 0.04 and anchors 0.10/0.16/0.22. The reference target has standard deviation {old['records']['reference']['target_std']:.3f}; the refined 128-grid target differs from reference by RMSE {old['records']['grid128_refined']['target_RMSE_vs_reference']:.3f}. This large sensitivity is why the original sharp-start trajectories cannot support a robust physical prediction claim here. Every reused archive hash and its exact original path is saved in <code>existing_data_audit.json</code>.</p>
<h2>Numerical checks, failures and limits</h2>
<p>Checks cover a manufactured affine stretch, rigid rotation, constant velocity, periodic interpolation, finite-difference verification of the sampled full gradient, a decaying Taylor–Green solution and Heun order, initial-field agreement across grids, divergence projection, disjoint groups, future-mutation invariance and marker relabeling. Across 64 runs, maximum normalized energy-budget residual is {max_energy:.3g} and maximum grid divergence is {max_div:.3g}. File hashes and source hashes verify successfully.</p>
<p>Mean paired target RMSE from grid 26³ to 38³ is {grid:.5f}; from timestep 0.005 to 0.0025 it is {time:.5f}. The grid target change exceeds the primary model's tiny RMSE advantage. Although the spectral fluid is divergence-free, trilinear interpolation is not exactly volume preserving: maximum accumulated |det M−1| over all runs is {maxvolume:.3f}, and the maximum refined-grid target-window volume error is {det:.3f}. The target is therefore FTLE of a numerical interpolated flow. Two spatial resolutions are a sensitivity check, not demonstrated continuum convergence. Small timestep changes also cross interpolation cell faces, where the gradient is discontinuous.</p>
<p>The ensemble is smooth, periodic, unforced and short. It does not establish generality to the user's original high-strain start, fully developed turbulence, walls, forcing, longer times or extreme Reynolds numbers. Uniform random markers and symmetric summaries may miss more specific arrangements. History-only features perform much worse than current strain information. A stronger nonlinear current-state predictor could absorb some history benefit. No predictive correlation establishes marker leadership, causal influence, regularity, singularity formation or a Navier–Stokes breakthrough.</p>
<p>For deterministic dynamics, the complete current velocity field with parameters and boundary conditions determines future evolution while a solution exists uniquely. A local velocity gradient is an incomplete observation of that state. History may proxy unobserved spatial structure or reduce model approximation error; this pilot cannot distinguish those explanations.</p>
<h2>Reproduce and inspect</h2>
<pre>python -m pip install -r requirements.txt
python -m pytest test_history_prediction.py -q
python simulate.py --pilot
python simulate.py --sensitivity
python evaluate.py
python audit_existing.py
python verify_results.py
python export_models.py
python report.py</pre>
<p>Run from <code>reports/matched-stretch/history-prediction</code>. Dependencies, source, locked protocol, per-run metadata, hashes, metrics, full predictions and PNG/SVG figures are supplied. The recorded-data folder contains all 64 numerical trajectory archives and metadata; the existing marker-stress folder contains the eight reused original archives. Fresh reproductions use data/ and preserve the immutable published recordings.</p>
<p class='small'>Method references: <a href='https://www.math.hkust.edu.hk/~masyleung/Reprints/leu11.pdf'>Leung (2011), finite-time Lyapunov exponent and flow map</a>; <a href='https://scikit-learn.org/stable/modules/cross_validation.html'>scikit-learn grouped cross-validation</a>; <a href='https://scikit-learn.org/1.8/common_pitfalls.html'>training-only preprocessing and leakage prevention</a>. These support definitions and evaluation practice; the numerical findings are from this study.</p>"""
    css="""body{font:16px/1.6 system-ui,-apple-system,Segoe UI,sans-serif;color:#202c3b;background:#f5f7fa;margin:0}main{max-width:1160px;margin:40px auto;background:white;padding:48px 56px;border:1px solid #dce3ec;border-radius:12px}h1{font-size:44px;line-height:1.14;letter-spacing:-1.5px;margin:16px 0}h2{font-size:25px;margin-top:44px;line-height:1.3}.lead{font-size:22px;color:#46566d;max-width:900px}.eyebrow{letter-spacing:1.5px;color:#2463a6;font-size:12px;font-weight:750}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:28px 0}.cards div{background:#eef4fa;border-radius:8px;padding:20px}.cards strong{display:block;font-size:32px;color:#2463a6}.cards span{font-size:14px}table{border-collapse:collapse;font-size:12px;width:100%}td,th{padding:9px 8px;text-align:right;border-bottom:1px solid #dde4ed;white-space:nowrap}th{background:#edf2f7}td:first-child,th:first-child{text-align:left}.scroll{overflow-x:auto}figure{margin:24px 0}img{width:100%;height:auto}figcaption,.small{font-size:13px;color:#68768b}a{color:#2463a6}code{font-size:13px;background:#eff2f6;padding:2px 4px}pre{padding:20px;background:#142337;color:#e8edf3;border-radius:8px;overflow:auto;font-size:13px}@media(max-width:700px){main{padding:24px;margin:0}.cards{grid-template-columns:1fr}h1{font-size:32px}}"""
    (out/"report.html").write_text("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Marker history prediction pilot</title><style>"+css+"</style><main>"+text+"</main></html>",encoding="utf-8")
    save_json(out/"headline.json",{"relative_RMSE_reduction_percent":pct,"A":A["macro"],"B":B["macro"],
                                   "paired":delta,"grid_target_RMSE":grid,"time_target_RMSE":time})
    print("Report and five figures generated",flush=True)


if __name__=="__main__":main()
