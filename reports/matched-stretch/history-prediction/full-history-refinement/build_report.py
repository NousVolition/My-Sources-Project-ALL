"""Report the frozen predictor study, keeping numerical and predictive claims separate."""
from pathlib import Path
import json
from collections import Counter
HERE=Path(__file__).resolve().parent


def main():
    r=json.loads((HERE/'results.json').read_text());p=r['protocol'];grids=r['grids'];fine=str(grids[-1])
    v=r['grid_models'][fine]; d=v['paired']['RMSE'];ap=v['paired']['AUPRC_AP']
    numerical='All predeclared operational numerical screens pass.' if r['all_screens_pass'] else 'Full-history numerical stability remains unestablished at the frozen tolerances.'
    gain=100*(-d['mean_B_minus_A'])/v['A_full']['macro']['RMSE']
    continuous=('The paired stretch-error interval includes zero; a predictive advantage remains unresolved.' if d['ci95'][0]<=0<=d['ci95'][1] else
        'The reused-run interval favors lower stretch error with history.' if d['ci95'][1]<0 else 'The reused-run interval favors the present-only predictor for stretch error.')
    lines=['# Full-history refinement with frozen predictors','',f'**{numerical}**',
        '',continuous+' These previously examined eight runs are exploratory numerical validation, not fresh statistical confirmation.','',
        '## What this tests','',
        'The earlier mirror/direction study refined short future branches from common saved states. This study instead integrates every marker history from the same original smooth initial velocity polynomial and initial marker positions, over the entire interval 0–1.2. Present features, past geometry and later target trajectories all change with resolution. No original source, predictor coefficient, training scaler, regularization, event threshold, or decision threshold was changed.','',
        'A uses 31 present features: velocity, the full 3×3 local velocity gradient, strain eigenvalues/norm, rotation norm, speed and present neighbor geometry. B uses exactly A plus the original 23 causal historical descriptors. Both are the original coarse-grid ridge/logistic fits, evaluated without refitting. This tests numerical robustness and transfer of those fitted predictors, not the best attainable fine-grid-trained predictor.','',
        'The unchanged parent incompressible Navier–Stokes solver uses joint Heun steps, periodic side-6 domain, viscosity 0.02, zero mean/force and the same trilinear velocity sampler and exact piecewise Jacobian. The smooth independent Fourier-mode seeds 1024–1031, all 512 initial markers and original labels are reused. They are distinct from the original unresolved sharp strain/tube field.','',
        'The main step is 0.0025, with two prechosen seeds checked at 0.00125 on the fine grids. Two existing grid-56/step-0.0025 histories are reused; the other six grid-56 main histories are new. The main ladder in this record is '+', '.join(str(n)+'³' for n in grids)+'. Conditional finer stages follow the [frozen protocol](protocol.json), with a maximum of 44 new simulations through 224³ and no automatic expansion beyond that limit. The grid-80/seed-1024 pilot appears once through cache reuse.','',
        'Prediction anchors are 0.3, 0.4, 0.5, 0.6 and 0.7; each history ends at its anchor and spans 0.2. The target starts strictly later, after a gap of 0.025, and accumulates logarithmic tangent stretch for 0.2. The fixed training-derived high-event cutoff is 1.0106081570657628. Locations are material marker start locations, not exact Eulerian peak location or first event onset.','',
        '## Frozen model results','',
        '| Main grid | A RMSE | B RMSE | B − A RMSE, paired 95% interval | A average precision | B average precision |',
        '|---|---:|---:|---|---:|---:|']
    for n in grids:
        x=r['grid_models'][str(n)];a=x['A_full']['macro'];b=x['B_history']['macro'];delta=x['paired']['RMSE']
        lines.append(f'| {n} | {a["RMSE"]:.6f} | {b["RMSE"]:.6f} | {delta["mean_B_minus_A"]:.6f} [{delta["ci95"][0]:.6f}, {delta["ci95"][1]:.6f}] | {a["AUPRC_AP"]:.6f} | {b["AUPRC_AP"]:.6f} |')
    lines+=['',f'At the finest grid, the mean relative stretch-error reduction is **{gain:.3f}%**. The event average-precision difference B−A is **{ap["mean_B_minus_A"]:.6f}**, interval [{ap["ci95"][0]:.6f}, {ap["ci95"][1]:.6f}]. Continuous and binary targets are separate outcomes. A stretch-error reduction cannot substitute for evidence of strong-event prediction benefit.','',
        'Intervals use 2000 paired bootstraps of the eight entire independent seed groups. All anchors and markers within a run stay together. Grids and time-step controls add no independent test runs. Sign-flip values in results are exploratory and unadjusted; no result-based model selection or threshold recalibration was performed. MAE, R², precision, recall, Brier score, log loss and probability calibration bins are retained in the full results.','',
        '## Numerical differences and the prediction advantage','',
        '| Adjacent main grids | Per-run target relative RMS range | Minimum event Jaccard | Present standardized RMS range | History standardized RMS range | Change in mean B−A RMSE |',
        '|---|---:|---:|---:|---:|---:|']
    for key,rec in r['spatial'].items():
        rows=rec['per_run']; ranges=[]
        for field,scale in [('relative_target_RMSE',100),('present_standardized_RMS',1),('history_standardized_RMS',1)]:
            vals=[scale*x[field] for x in rows];ranges.append(f'{min(vals):.4f}–{max(vals):.4f}'+('%' if scale==100 else ''))
        lines.append(f'| {key} | {ranges[0]} | {min(x["event_Jaccard"] for x in rows):.5f} | {ranges[1]} | {ranges[2]} | {rec["paired_metric_drift"]["RMSE"]["mean_fine_minus_coarse"]:.6g} |')
    lines+=['','Standardized descriptor changes use the original training scalers and include changes in neighborhood membership. They are descriptive measurements, not additional tuned pass/fail tests. Each descriptor’s change and each predictor’s continuous/probability output drift are saved per run.','',
        'The finest-pair target relative RMS must be at most 1% in every run, high-event Jaccard at least 0.95, and both predictors’ output relative RMS differences at most 1%. Successive target RMS differences must contract. Because the old mean A/B difference was only about 0.000364, a 1% target screen alone is inadequate: the absolute change in mean B−A RMSE must also be at most 0.00009 (approximately one quarter of that previously observed difference), and each run’s change at most 0.00018. Mean event-AP contrast drift must be at most 0.002 and each run’s at most 0.01.','',
        'The finest-grid half-step target difference must be at most one tenth of its same-run adjacent spatial difference. Its A/B RMSE contrast change must be at most 0.00001 and AP contrast change at most 0.001 for each checked run. These bounds were frozen before the pilot; they bound measured resolution sensitivity and do not provide a continuum error bound, significance test, or practical-benefit guarantee.','',
        '| Frozen screen | Passed / evaluated |','|---|---:|']
    counts=Counter(x['kind'] for x in r['screens']);passed=Counter(x['kind'] for x in r['screens'] if x['pass'])
    for key,count in counts.items():lines.append(f'| {key} | {passed[key]} / {count} |')
    if r['failures']:
        lines+=['','Every failed screen is preserved; none was relaxed.','',
                '| Screen and identifier | Value | Required |','|---|---:|---:|']
        for x in r['failures']:lines.append(f'| {x["kind"]}: {x["identifier"]} | {x["value"]:.6g} | {x["comparison"]} {x["limit"]:.6g} |')
    lines+=['','## Time-step and interpolation diagnostics','',
        '| Half-step grid | Maximum target RMS difference | Maximum absolute per-run B−A RMSE change | Maximum absolute per-run B−A AP change |',
        '|---|---:|---:|---:|']
    for n,x in r['temporal'].items():lines.append(f'| {n} | {max(y["target_RMSE"] for y in x["per_run"]):.6g} | {max(abs(y) for y in x["paired_metric_drift"]["RMSE"]["per_run_fine_minus_coarse"]):.6g} | {max(abs(y) for y in x["paired_metric_drift"]["AUPRC_AP"]["per_run_fine_minus_coarse"]):.6g} |')
    lines+=['','| Main grid | Maximum CFL | Maximum absolute relative energy residual | Mean interpolated divergence RMS | Maximum accumulated probe volume error |',
        '|---|---:|---:|---:|---:|']
    for n in grids:
        rows=[x for x in r['files'] if x['config']['n']==n]
        lines.append(f'| {n} | {max(x["max_CFL"] for x in rows):.5f} | {max(x["max_energy_residual"] for x in rows):.3e} | {sum(x["interpolated_divergence_RMS"] for x in rows)/8:.6g} | {100*max(x["max_probe_det_error"] for x in rows):.3f}% |')
    lines+=['','The spectral field remains divergence-free within the saved numerical checks, while its trilinear interpolant is not exactly divergence-free and its Jacobian jumps at cell faces. Accumulated probe-volume error and history-feature changes remain visible. Shrinking grid differences do not justify an unsupported Richardson extrapolation or formal order.','',
        '## Frozen historical-input sensitivity controls','',
        '| Input supplied to the existing B predictor | Mean RMSE | Mean average precision |','|---|---:|---:|']
    for key,x in r['frozen_history_controls'].items():lines.append(f'| {key} | {x["macro"]["RMSE"]:.6f} | {x["macro"]["AUPRC_AP"]:.6f} |')
    lines+=['','These controls keep present inputs fixed and change only historical descriptors: reverse the past sequence, use the spatially farthest marker’s history, or permute rows within each run/anchor. The original B coefficients are reused without adapting to these inputs. This measures sensitivity of an already fitted representation; it is not an equally refitted negative-control comparison and cannot establish chronology’s necessity or causality. The parent consistent-label, causality and fitting-leakage tests remain applicable.','',
        '## Provenance and reproduction','',
        '[input-verification.json](input-verification.json) records 63 source/model/archive/metadata files matched to GitHub before the pilot. [launch-verification.json](launch-verification.json) records frozen local source hashes. Per-run metadata records executed sources, initial-condition seed, configuration, array checksum and numerical diagnostics. The recorder’s complete outputs match the original parent simulation bit-for-bit in an independent small-run test; a separate analytic transfer test verifies the same initial polynomial across non-nested grids. All 34 parent and new analytic/regression tests passed, recorded in [tests.xml](tests.xml).','',
        'The [execution record](execution.json), [results](results.json), [per-run comparisons](per-run-comparisons.csv), [saved predictors’ outputs and feature arrays](predictions-and-features.npz), [raw histories](recorded-data), stage decisions and logs preserve completed and failed screens. [publication-manifest.json](publication-manifest.json) lists published checksums. Recovery checkpoints are local continuation evidence, excluded from bulk publication. Numerical sources and fitted models were not edited while active.','',
        'From the repository root, with the parent requirements installed:','',
        '```text','python -B reports/matched-stretch/history-prediction/full-history-refinement/analyze_history.py',
        'python -B reports/matched-stretch/history-prediction/full-history-refinement/build_report.py','```','',
        'Analysis reuses saved histories; it does not launch simulations or train models. Simulation caches require exact executed study-source hashes. Parent analysis accepts only LF/CRLF-equivalent source bytes; actual source edits are rejected. Preserve partial checkpoints and inspect process ownership before any restart. Stage snapshots and their failures are historical evidence, not extra independent observations.','',
        '## What remains necessary','',
        'A numerically stable frozen-model score does not establish predictive usefulness. The eight test seeds have already been examined. A new confirmation would need independently prepared training/validation/test runs at resolved accuracy, fixed comparison and practical-benefit criteria, matched model families and negative controls, and run-level uncertainty on untouched test results. No such fresh study was silently substituted here. If the interval includes zero, that is unresolved evidence; ruling out a useful benefit requires an upper confidence bound below a predeclared useful-effect size.','',
        'History may help an incomplete observation by summarizing omitted spatial information. It adds no fundamental information beyond a complete deterministic present fluid state. Passive marker records do not change fluid stress or establish material memory, leadership, or a Navier–Stokes breakthrough. The original nonsmooth periodic strain and off-central initial maxima remain separate limitations of the untouched original matrix.','',
        '![Full-history target stability and frozen A/B scores](convergence.png)','']
    (HERE/'README.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':main()
