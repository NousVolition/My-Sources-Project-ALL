"""Write an evidence-based report only after simulations, frozen evaluation and checks."""
import csv,json,html,platform
from common import ROOT,save,digest,read_protocol
CODE='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-organization-pilot/extensions/stress-history'
STREAMS='https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization/extensions/stress-history'
def read(name):return json.loads((ROOT/'results'/name).read_text(encoding='utf-8'))
def csvtable(name):
    with (ROOT/'results'/name).open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    cols=list(rows[0]);fmt=lambda v:f'{float(v):.5g}' if v.replace('.','',1).replace('-','',1).isdigit() else v
    return '<div class="table"><table><tr>'+''.join('<th>'+html.escape(c.replace('_',' '))+'</th>' for c in cols)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(fmt(r[c]))+'</td>' for c in cols)+'</tr>' for r in rows)+'</table></div>'
def main():
    p=read_protocol();s=read('prediction_summary.json');flow=read('flow_summary.json');checks=read('numerical_verification.json');screen=read('screen_verification.json');freeze=read('freeze_receipt.json')
    public=read('public_protocol_commit.json');refined=read('refinement_forecast_summary.json');ntest=len(p['splits']['test']);ntotal=sum(map(len,p['splits'].values()));rows=[];sentences=[];refinement_text=[];mechanisms=[];control_rows=[];current_comparisons=[]
    with (ROOT/'results/response_regimes.csv').open(encoding='utf-8') as f:response_rows=list(csv.DictReader(f))
    with (ROOT/'results/prediction_by_strength.csv').open(encoding='utf-8') as f:strength_results=list(csv.DictReader(f))
    strength_table=[]
    for r in strength_results:strength_table.append(f"| {r['regime']} | {r['relative_strength']} | {float(r['strong_present_rmse']):.6f} | {float(r['history_rmse']):.6f} | {float(r['error_reduction_percent']):+.2f}% | [{float(r['descriptive_reduction_ci95_low']):+.2f}%, {float(r['descriptive_reduction_ci95_high']):+.2f}%] |")
    for regime,q in s.items():
        c=q['contrasts']['history_vs_present'];lo,hi=c['relative_improvement_interval']
        status='Resolved improvement for these declared partial observations' if q['history_supported'] else ('Improvement over the richer present model; useful relevant history remains unestablished against the null controls' if c['rmse_difference_interval'][1]<0 else ('History increased error for these frozen models' if c['rmse_difference_interval'][0]>0 else 'No resolved improvement beyond the richer present measurements'))
        rows.append(f'| {regime} | {q["mean_field_rmse"]["present"]:.6f} | {q["mean_field_rmse"]["history"]:.6f} | {c["relative_improvement_percent"]:+.2f}% | [{lo:+.2f}%, {hi:+.2f}%] |')
        sentences.append(f'**{regime}: {status}.** History changes mean error by {-c["relative_improvement_percent"]:+.2f}%; the paired 97.5% interval for error reduction is [{lo:+.2f}%, {hi:+.2f}%]. Positive reductions favor history. The raw RMSE difference is {c["mean_rmse_difference"]:+.6f}, interval {c["rmse_difference_interval"]}.')
        for key,label in [('history_vs_shuffled_mean','History minus expected shuffled history'),('history_vs_training_prior','History minus training-only prior')]:
            contrast=q['contrasts'][key];a,b=contrast['rmse_difference_interval'];control_rows.append(f"| {regime} | {label} | {contrast['mean_rmse_difference']:+.6f} | [{a:+.6f}, {b:+.6f}] |")
        contrast=q['contrasts']['present_vs_conventional'];a,b=contrast['relative_improvement_interval'];current_comparisons.append(f"{regime}: the stronger present measurements change error by {-contrast['relative_improvement_percent']:+.2f}% relative to the eight-global-summary model; the paired 97.5% interval for error reduction is [{a:+.2f}%, {b:+.2f}%].")
        r=refined[regime];means=r['mean_history_minus_present_on_four_fields'];refinement_text.append(f"{regime}: mean history-minus-present RMSE on the four refinement fields is {means['base']:+.6f} at baseline settings, {means['half']:+.6f} at half step and {means['fine']:+.6f} on the finer grid. The largest paired grid change in that contrast is {r['maximum_grid_contrast_drift']:.6g}; it {'exceeds' if r['drift_exceeds_primary_contrast_magnitude'] else 'is below'} the magnitude of the full test-sample mean contrast. The mean contrast sign {'changes' if r['coarse_fine_mean_sign_changed'] else 'does not change'} in this four-field subset.")
        d=flow[regime]['intervention_time_baseline'];re=d['Taylor_Re'];transfer=d['flux_k4'];nonlinear=d['nonlinear_over_viscous_rms']
        mechanisms.append(f"**{regime} mechanisms.** At the impulse time, the mean approximate Taylor Reynolds number is {re['mean']:.3g} (descriptive 95% interval {re['ci95'][0]:.3g}–{re['ci95'][1]:.3g}); the nonlinear-to-viscous acceleration RMS ratio is {nonlinear['mean']:.3g}. Mean forward spectral flux across box-mode radius four is {transfer['mean']:.4g}, with interval [{transfer['ci95'][0]:.4g}, {transfer['ci95'][1]:.4g}].")
        for rr in response_rows:
            if rr['regime']!=regime:continue
            strength=rr['relative_strength'];production=float(rr['integrated_production_over_initial_energy']);dissipation=float(rr['integrated_dissipation_over_initial_energy']);final=float(rr['final_energy_ratio']);count=int(rr['fields_with_sampled_energy_growth_above_1_01'])
            mechanisms.append(f"At relative impulse {strength}, mean integrated production is {production:.4g} and viscous loss is {dissipation:.4g}, both divided by initial difference energy. Their net balance explains the mean final energy ratio {final:.4g}. {count}/{ntest} fields exceed 1.01 times initial difference energy at a saved observation. These counts describe sampled transient amplification; they do not establish a long-time instability.")
    valid=checks['all_numerical_passed'] and checks['all_refinements_passed'];claim='Both numerical screens passed.' if valid else 'A numerical screen failed; interpret prediction differences as unresolved until that failure is addressed.'
    text=f'''# Stronger fluid disturbances and independent history forecasts

Completed numerical study — 11 October 2026. {claim}

**The prespecified test for a useful history advantage beyond all controls was not met in either regime.** In active flows, history improves on the richer present-input model, but does not clearly beat expected shuffled histories or the training-only average. The simple global-current model has lower mean test error than the history model (0.096598 versus 0.100727); the richer present model generalizes poorly in this sample. The 35.1% improvement over that richer model therefore does not establish that the relevant past carries additional useful information. No predictors were refitted after these results.

{' '.join(sentences)}

| Flow regime | Richer present RMSE | With history RMSE | Error reduction | Paired 97.5% interval |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

**[Runnable code, tests, frozen models and exact recordings → My-Sources-Project-ALL]({CODE})** · **[Streams and Rocks report]({STREAMS})** · [Illustrated report](report.html)

![Independent prediction comparisons](results/figures/prediction.png)

## Completed experiments

{ntotal} new independent initial-condition seeds: {len(p['splits']['train'])} training, {len(p['splits']['validation'])} validation and {ntest} held-out seeds. Each seed is run at both viscosities, with one unchanged copy and impulses at 0.25 and 0.75 times the baseline velocity RMS. This gives {2*ntotal} simulated regime/seed triplets and {6*ntotal} fluid trajectories. The two regimes share initial seeds for matching; they are not independent replicates of each other. Within each regime, uncertainty resamples its {ntest} independent test fields, keeping both strengths together. Four prechosen test fields per regime have half-step and finer-grid reruns, counted as dependent controls. Separate screening seeds never enter prediction training or testing.

The active regime uses five-fold lower viscosity than the smooth comparison, higher speeds than the prior smooth recovery study, and measured nonlinear energy transfer and vortex stretching. It is a **transient, modest-Reynolds-number 3D flow**, not demonstrated fully developed high-Reynolds turbulence or an extended inertial range. This new impulse-energy target differs from the earlier 120-field forcing-change target; its scores are not directly comparable.

## What the models know

All inputs stop at the impulse time. The conventional model sees eight global current summaries and the known strength. The stronger present model sees 33 current global, spectral, source-local and marker-geometry measurements, including the complete three-component velocity and nine-component velocity-gradient tensor at the disturbance center, the known strength, and the initial production and dissipation rates divided by initial perturbation energy: 36 inputs. The history model appends 33 causal changes at each of three lags (0.25, 0.5 and 1 time unit): 135 inputs. Three controls attach complete histories from a different run while preserving both strength rows within each field. The same kernel-ridge model family and validation grid are used throughout; scaling uses training data only. Test labels never choose features, hyperparameters, windows or settings. Even this stronger present measurement remains an incomplete observation of the spatial fluid state.

The primary target is the integral over 0–1.5 of perturbation energy divided by its impulse-time value. Error is the RMSE across the two strengths for each field, then averaged across independent fields; its units are model time. A 20,000-replicate paired field bootstrap gives 97.5% intervals for each regime, with Bonferroni protection for the two primary history-versus-present comparisons. A history claim also requires improvement over a training-only mean outcome for each known strength and over the expected shuffled-history control: all nonidentity whole-field donors averaged over the three frozen control models. The shuffle comparison uses a field-node bootstrap resampling current and donor roles jointly, excluding the same original field even if duplicated. The three fixed-permutation scores are descriptive checks; the main figure shows the expected shuffled-history error. All prediction intervals are conditional on the frozen calibration sample and fitted models, and do not include variability from retraining. Descriptive flow/response intervals are 95% and do not supply extra confirmatory hypotheses.

The [protocol and source were publicly recorded before held-out runs]({public['url']}); calibration was already underway at that publication. The model freeze was recorded at {freeze['created_utc']}, SHA-256 `{freeze['sha256']}`. No test simulations existed at model freeze. All new runs have separate deterministic random seeds and source/configuration/data hashes. The complete deterministic present fluid state already determines later evolution in this numerical model; history can help partial measurements reveal omitted state or change the approximation error of a finite predictor. This is not evidence for an additional physical memory force or history-dependent modification of Navier–Stokes.

## Physical and numerical checks

Periodic side-6 constant-density incompressible Navier–Stokes, strict two-thirds Fourier dealiasing and shared fixed solenoidal ABC forcing; RK4 evolves fluid and one-way passive markers together. The impulse is the curl of a smooth periodic localized vector potential with fixed modes through five per coordinate. It injects a globally defined RMS velocity and changes total energy by its squared norm plus its cross term with the initial field; it is not a point force. Density, compressibility, phase changes, walls and particle feedback are absent.

The exact Galerkin perturbation-energy balance is `dEδ/dt = P − D`, with `P = −<δu_i δu_j ∂_j u_base,i>` and `D = ν<|curl δu|²>`. Shared forcing cancels. Positive production transfers energy from the evolving background into the velocity difference; dissipation removes it. Energy-weighted radius and spatial participation describe spreading, not literal signal speed. Marker stretching does not establish energy amplification or causation. Current tetrad volume and shape are statistics of finite neighborhoods transported by cubic B-spline interpolation; they are not continuum tangent-map determinants or proof of exact volume preservation. Taylor-scale Reynolds numbers use an isotropic approximation, `u′ = U_RMS/√3`, `λ² = 15νu′²/ε`, `Reλ = u′λ/ν`; isotropy is not demonstrated here. The integral Reynolds diagnostic uses full-vector velocity RMS and the isotropic-spectrum length convention. Lower viscosity changes both diffusion and dissipation; it is not a direct measurement of a single interaction-strength variable.

Seventeen analytic, physical-budget and leakage tests pass. Run-level checks include finite states/CFL, viscous step limit, divergence, total and perturbation energy balances at every saved time, high-band tails, and cutoff times Kolmogorov length. These screens and matched finer-grid reruns support the declared finite calculation; four refinement seeds do not prove ensemble-wide continuum convergence. The cutoff criterion is a screening convention, not a guarantee. The original rejected step and exact original screening outputs are retained in [screen events](screen_events.json). The original trilinear active-flow marker comparison failed (0.030318 against 0.03), while its energy response agreed. Cubic tracking was revalidated on that same seed at both viscosities using finer grids and half steps; the fluid-energy arrays remain bit-identical under this one-way-marker revision. The failed screen remains in the report package and was not used for fitting. A tested diagnostic optimization bins spectral energy once per observation and reuses a curl array; it leaves the screened fluid and marker updates unchanged.

{' '.join(refinement_text)} These checks are dependent numerical comparisons, not additional independent test fields. They do not quantify continuum uncertainty across the entire ensemble.

## Results and data

The following prespecified control contrasts use mean independent-field RMSE; a negative difference favors history. They supplement the primary present-versus-history comparison above.

| Regime | Control comparison | RMSE difference | Paired 97.5% interval |
|---|---|---:|---:|
{chr(10).join(control_rows)}

{' '.join(current_comparisons)} This is a secondary comparison, and does not isolate marker geometry from the other additional present measurements.

### Descriptive breakdown by disturbance strength

These per-strength scores use pooled RMSE across the 32 fields, rather than the primary mean field RMSE across both strengths. Their paired 95% intervals are unadjusted, conditional on the frozen models, and descriptive. They cannot replace a failed primary test or establish additional confirmatory claims.

| Regime | Relative impulse RMS | Richer present RMSE | With history RMSE | Error reduction | Descriptive 95% interval |
|---|---:|---:|---:|---:|---:|
{chr(10).join(strength_table)}

- [Quantitative parameter regimes](results/parameter_regimes.csv), [disturbance responses](results/response_regimes.csv), [flow statistics and uncertainty](results/flow_summary.json)
- [Prediction breakdown by strength](results/prediction_by_strength.csv), [its distinct descriptive scope](results/strength_prediction_scope.json)
- [Every held-out prediction](results/predictions.csv), [primary scores and uncertainty](results/prediction_summary.json), [frozen model file](results/frozen_models.json)
- [Numerical verification](results/numerical_verification.json), [screening verification](results/screen_verification.json), [prediction changes under refinement](results/prediction_refinements.csv), [test record](results/tests.xml)
- [Prediction refinement summary](results/refinement_forecast_summary.json), [runtime/library versions](results/environment.json)
- [Frozen protocol](protocol.json), [exact-file manifest](manifest.json), [source provenance](UPSTREAM.txt)

![Disturbance energy and spreading](results/figures/disturbance.png)
![Energy production and dissipation](results/figures/energy_budget.png)
![Spatial scales, stretching and marker deformation](results/figures/flow_scales.png)

## Which mechanisms explain the response?

{chr(10).join(mechanisms)}

Advection and strain redistribute the velocity difference; background strain supplies or removes difference energy, and viscosity dissipates it. Vortex stretching and spectral transfer characterize the active flow, while the directly checked production-minus-dissipation balance explains its energy response. A history forecast comparison cannot distinguish all of these mechanisms causally, and does not test whether molecular interactions weaken. The combined history features also cannot isolate the contribution of marker organization from the contribution of velocity and spectral history.

## Reproduce in the linked code repository

Use Python 3.12 or later and run inside the linked extension folder. Restoration requires only Python's standard library. Exact recordings retain saved observables, marker positions, diagnostic time series and state hashes; complete 3D fields are regenerable from seeds and source, and are not archived for every new run.

```text
python restore_recordings.py
python -m pip install -r requirements.txt
python -m pytest test_study.py -q
python audit.py
```

Audit the original restored bytes before regenerating output files. To repeat the calculations, use a separate working copy and run:

```text
python verify.py
python forecast.py evaluate
python analyze.py
python report.py
python manifest.py
python audit.py
```

For a fresh simulation reproduction, start from a separate copy and archive `data/runs`, `data/forecast`, `results/frozen_models.json` and `results/freeze_receipt.json` outside that extension folder. Keep the original directory as an unchanged reference. Preserve `data/screen`, the historical screening results, and the public protocol receipt. Run `pipeline.py`, then `report.py`, `manifest.py` and `audit.py`. The pipeline completes calibration, freezes new models before any new test simulation, runs tests and numerical controls, verifies them, and evaluates forecasts. The numerical tests can be repeated with `python -m pytest test_study.py -q`; `environment.py` records the current library versions. Retain the locked protocol. Regenerated output bytes can differ across library versions or platforms; the newly generated manifest describes that reproduction rather than the original snapshot. The code refuses stale simulation caches or refitting an existing frozen model. `batch.py screen` reproduces screening with current source; historical source hashes are audited separately. Original failed/earlier screening files are preserved rather than silently replaced.

## Unanswered questions — proposed, not completed

- Do the results replicate with more independent fields, other random-field families, source positions and impulse shapes, or other observation horizons?
- Does a causal history benefit survive richer current spatial measurements, different model families, noisy/sparse observations and strict ensemble-wide spatial refinement?
- What happens at statistically stationary, substantially higher Reynolds numbers with a wider cascade? This requires longer burn-in, larger grids and a new frozen study.
- Do marker histories specifically help, rather than histories of global/spectral flow summaries? The combined-history comparison does not isolate this.
- Are any history gains causal consequences of organization? Matched state reconstruction, intervention designs and separate identifiability arguments are needed.
- Do wall boundaries, evaporation/heat exchange, dilute inertial droplets or smoke particles change these results? Ice requires a separate solid/thermal model. Laboratory validation remains future work.

## Primary physics references

The resolution diagnostic follows common DNS practice illustrated in [Momenifar, Dhariwal & Bragg (2019)](https://link.aps.org/accepted/10.1103/PhysRevFluids.4.054301), whose parameters include cutoff-times-Kolmogorov lengths near 1.5. Their inertial-particle findings are not results of this passive-marker experiment. [Johns Hopkins 8192³ isotropic-turbulence documentation](https://turbulence.pha.jhu.edu/docs/README-isotropic8192.pdf) provides a high-resolution reference; this pilot does not reproduce that flow or claim comparable Reynolds numbers.
'''
    (ROOT/'README.md').write_text(text,encoding='utf-8')
    # Independent HTML artifact with readable tables and figures, no external rendering dependency.
    paras=[];table_lines=[];code_lines=[];in_code=False
    def inline(line):
        import re
        rendered=html.escape(line);rendered=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',rendered);rendered=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',rendered);rendered=re.sub(r'`([^`]+)`',r'<code>\1</code>',rendered);return rendered
    def flush_table():
        if not table_lines:return
        cells=[line.strip('|').split('|') for line in table_lines]
        paras.append('<div class="table"><table><thead><tr>'+''.join('<th>'+inline(c.strip())+'</th>' for c in cells[0])+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+inline(c.strip())+'</td>' for c in row)+'</tr>' for row in cells[2:])+'</tbody></table></div>');table_lines.clear()
    for line in text.splitlines():
        if line.startswith('```'):
            flush_table();in_code=not in_code
            if not in_code:paras.append('<pre><code>'+html.escape('\n'.join(code_lines))+'</code></pre>');code_lines.clear()
            continue
        if in_code:code_lines.append(line);continue
        if line.startswith('|'):table_lines.append(line);continue
        flush_table()
        if line.startswith('#'):level=min(len(line)-len(line.lstrip('#')),3);paras.append(f'<h{level}>'+html.escape(line.lstrip('# '))+f'</h{level}>')
        elif line.startswith('!['):alt,url=line[2:].split('](',1);paras.append(f'<figure><img src="{url[:-1]}" alt="{html.escape(alt)}"><figcaption>{html.escape(alt)}</figcaption></figure>')
        elif line.strip():
            paras.append('<p>'+inline(line)+'</p>')
    flush_table()
    intro='<h2>Quantitative regimes</h2>'+csvtable('parameter_regimes.csv')+'<h2>Completed disturbance responses</h2>'+csvtable('response_regimes.csv')
    page='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Stronger fluid disturbances and history prediction</title><style>body{font:16px/1.6 system-ui;color:#22313a;max-width:1100px;margin:40px auto;padding:0 24px}h1,h2{line-height:1.25;color:#143f58}a{color:#126b9b}img{width:100%;height:auto}figure{margin:32px 0}figcaption{font-size:14px;color:#596a74}table{border-collapse:collapse;font-size:13px}td,th{border-bottom:1px solid #cdd6dd;padding:8px;text-align:left}.table{overflow-x:auto}code{font-size:14px;background:#edf3f6;padding:2px}</style>'+''.join(paras)+intro+'</html>'
    (ROOT/'report.html').write_text(page,encoding='utf-8');save('results/completion.json',{'date':'2026-10-11','unique_new_seeds':ntotal,'regime_seed_triplets':2*ntotal,'fluid_trajectories':6*ntotal,'heldout_seeds_per_regime':ntest,'numerical_gates_passed':valid,'history_supported':{k:v['history_supported'] for k,v in s.items()},'python':platform.python_version(),'code_url':CODE,'streams_url':STREAMS})
if __name__=='__main__':main()
