"""Write completed findings from recorded values; never substitutes predictions."""
import csv
from common import ROOT,read,protocol,save,now
from html_report import render

CODE='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-organization-pilot/extensions/higher-re-history'
STREAMS='https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization/extensions/higher-re-history'

def main():
    p=protocol();r=read('results/prediction_summary.json');s=read('results/test_stationarity.json');cal=read('results/calibration_stationarity.json');physics=read('results/physical_summary.json');checks=read('results/numerical_checks.json');ref=read('results/refinement_prediction_summary.json');freeze=read('results/freeze_receipt.json');execution=read('results/execution.json')
    valid=checks['all_passed'] and s['passed'] and cal['passed'];support=valid and r['history_supported']
    conclusion=('Matched history improved prediction beyond every declared current-state and null comparator.' if support else 'The declared history features did not establish a useful prediction advantage beyond all the current-state and null comparators.')
    if not valid:conclusion='Numerical or stationarity gates failed, so the prediction comparison is exploratory. '+conclusion
    primary=r['contrasts']['present'];lo,hi=primary['reduction_interval'];rl,rh=s['re_lambda_interval']
    comparators=[]
    for name,label in [('present','Richer current measurements'),('global_present','Global current measurements'),('wrong_history','Wrong-run histories'),('training_prior','Training average by strength'),('event_prior','Training average by time and strength')]:
        q=r['contrasts'][name];e=q.get('null_mean_run_rmse',r['mean_run_rmse'].get(name));a,b=q['reduction_interval'];comparators.append(f"| {label} | {e:.6f} | {q['error_reduction_percent']:+.2f}% [{a:+.2f}%, {b:+.2f}%] |")
    station=[]
    for key,label in [('energy','Kinetic energy'),('dissipation','Dissipation'),('integral_length','Integral length')]:
        q=s[key];a,b=q['drift_interval'];station.append(f"| {label} | {q['early_mean']:.5f} | {q['late_mean']:.5f} | {100*q['relative_drift']:+.2f}% [{100*a:+.2f}%, {100*b:+.2f}%] |")
    pilots=[]
    for q in checks['screening']['long_grid_pairs']:
        z=q['statistics'];pilots.append(f"| {q['seed']} | {z['energy']['base']:.4f} to {z['energy']['fine']:.4f} | {100*z['energy']['relative_change']:.2f}% | {100*z['dissipation']['relative_change']:.2f}% | {100*z['re_lambda']['relative_change']:.2f}% |")
    physical=[]
    for q in physics['responses']:physical.append(f"| {q['relative_impulse']:.2f} | {q['mean_final_energy_ratio']:.4f} [{q['final_ratio_ci_low']:.4f}, {q['final_ratio_ci_high']:.4f}] | {q['mean_production_integral']:.4f} | {q['mean_dissipation_integral']:.4f} | {q['mean_forcing_difference_integral']:.4f} |")
    with (ROOT/'results/prediction_by_strength.csv').open(encoding='utf-8') as f:subgroups=list(csv.DictReader(f))
    subrows=[f"| {float(q['strength']):.2f} | {float(q['reduction_percent']):+.2f}% [{float(q['ci_low']):+.2f}%, {float(q['ci_high']):+.2f}%] |" for q in subgroups]
    grid=max(q['max_curve_relative_change'] for q in checks['refinements'] if q['label']=='fine');step=max(q['max_curve_relative_change'] for q in checks['refinements'] if q['label']=='half')
    old_re=p['upstream']['previous_measured_re_lambda']
    text=f'''# History prediction at increased Reynolds number

**Completed finding: {conclusion}** The experiment has {len(p['splits']['test'])} independent test flows, separate from 24 training and 8 validation flows. Four intervention times and two strengths produce 256 dependent targets. Matched-history mean run RMSE is **{r['mean_run_rmse']['history']:.6f}**; the richer current-state RMSE is **{r['mean_run_rmse']['present']:.6f}**. History reduces error by **{primary['error_reduction_percent']:+.2f}%**, with a paired 95% interval **[{lo:+.2f}%, {hi:+.2f}%]**. Positive reductions favor history.

The measured approximate Taylor Reynolds number is **{s['re_lambda_mean']:.2f} [95% interval {rl:.2f}, {rh:.2f}]**, compared with **{old_re:.2f}** in the [earlier stationary study](../stationary-history/README.md), a factor of {s['re_lambda_mean']/old_re:.2f}. The new study reduces viscosity from 0.025 to {p['nu']} while retaining the forcing, initial-speed distribution, observations and interventions. New independent seeds supply every training and test observation. This controlled extension reaches a higher Reynolds number; its spectrum and finite resolution do not establish high-Reynolds asymptotic behavior or a broad inertial range.

[Runnable code, tests and restorable recordings]({CODE}) · [Streams and Rocks]({STREAMS})

## Completed experiment

| Item | Executed setting |
| --- | --- |
| Equations and boundaries | Three dimensional incompressible Navier–Stokes in a periodic cube of side 6 |
| Viscosity and power | {p['nu']} and {p['power']} in dimensionless units |
| Main and finer grid | {p['n']} cubed and {p['fine_n']} cubed, strict two-thirds Fourier cutoff |
| Integration | Double precision RK4; baseline step {p['spin_dt']}, response step {p['response_dt']} |
| Initial conditions | 64 independent low-mode, divergence-free Gaussian draws, each normalized to RMS speed {p['initial_speed']} |
| Disturbances | Localized divergence-free impulses, ratios 0.25 and 0.75 of baseline RMS speed, at times 32, 36, 40 and 44 |
| Response | Identical initial copies; impulse added to one; horizon 1.5 and saved spacing 0.05 |
| Histories | Eulerian measurements at the disturbance location and global quantities; lags 0.5, 1 and 2 |
| Full-duration grid pilots | Two separate initial fields, each integrated on both grids through time 44 |
| Test refinement | Four preselected flows at time 44; matched restarts at time 42, half steps and finer grid at those same smaller steps |

The disturbance is a fixed Fourier polynomial derived from a smooth periodic vector potential. Its strength is a global RMS velocity ratio, not a fraction of total kinetic energy. Constant-power negative damping forces the large-scale modes with mode length at most 2. Each copy evaluates the same forcing law on its own state, so their forcing fields can differ. No particle interactions, phase changes or calibrated material properties enter this idealized fluid experiment.

The [acquisition record](results/execution.json) identifies the immutable workflow commit, frozen artifact, completed runs and platforms. The {execution['local_test_count']} local tests and {execution['linux_test_count']} Linux tests passed. Original integration bytes and hashes are retained; cached copies are never counted as extra independent realizations.

## Resolution and stationarity checks

All declared numerical gates **{'passed' if checks['all_passed'] else 'did not pass'}**. The minimum sampled cutoff times Kolmogorov length among test baselines and disturbed copies is **{physics['minimum_keta']:.4f}**, against 1.5. The cutoff is the inscribed axis cutoff `(ceil(N/3)-1) times 2 pi/6`; Kolmogorov length is `(viscosity cubed/dissipation)^(1/4)`. The maximum response CFL is {physics['maximum_response_cfl']:.4f}, and the maximum normalized difference-energy budget residual is {physics['maximum_difference_budget_residual']:.3g}. Every training, validation and test flow is checked, along with the controls.

Across the four predeclared final-event controls, the maximum normalized response-curve change is **{100*step:.5f}%** for halved steps and **{100*grid:.5f}%** for the finer grid at the same smaller steps. The maximum change in the frozen history-minus-current RMSE contrast is {ref['max_absolute_step_contrast_change']:.6g} for step and {ref['max_absolute_grid_contrast_change']:.6g} for grid. Other intervention times have numerical diagnostics but no matched restart controls. These checks test short-time response resolution, not population-wide continuum extrapolation.

![Matched time and grid controls](results/figures/resolution.png)

The full-duration pilots provide separate statistical grid comparisons. The table gives absolute relative differences between means over time 20 to 44, for two paired initial fields. Changes in both grid and step are combined here; chaotic trajectories need not remain pointwise close over the whole integration. The two pilot pairs do not provide a precise population confidence interval.

| Pilot field | Mean energy at 128 then 160 | Energy difference | Dissipation difference | Taylor Reynolds difference |
| --- | --- | --- | --- | --- |
{chr(10).join(pilots)}

Calibration stationarity **{'passed' if cal['passed'] else 'failed'}** before model freezing, and independent test stationarity **{'passed' if s['passed'] else 'failed'}**. The first intervention follows at least {s['minimum_spin_turnovers']:.2f} accumulated integral turnover times. The windows 20 to 32 and 32 to 44 compare statistics across independent flows; time samples are not independent replicates.

| Quantity | Earlier window | Later window | Relative drift and paired 95% interval |
| --- | --- | --- | --- |
{chr(10).join(station)}

The declared drift-equivalence limit is ±15%. Mean dissipation/input mismatch is {100*s['power_balance']['relative_mismatch']:+.2f}%, with its 95% interval required to lie within ±10%. The normalized-spectrum total variation is {100*s['spectrum']['total_variation']:.2f}%, and its 95% upper bound is {100*s['spectrum']['upper_95_percent_bound']:.2f}% against 10%. Forward flux through mode length 4 is {s['flux_mean']:.5f}, with a 95% interval [{s['flux_interval'][0]:.5f}, {s['flux_interval'][1]:.5f}]; the lower bound must exceed 5% of power input. The checks support approximate finite-window statistical stationarity.

![Stationarity and spectrum across independent flows](results/figures/stationarity.png)

The Taylor Reynolds estimate uses the isotropic approximation `sqrt(15) times (2E/3) / sqrt(viscosity times dissipation)`. Component mean-square velocities relative to their isotropic average are {', '.join(f'{x:.3f}' for x in s['component_ratio_to_isotropic'])}; exact isotropy of every realization is not assumed. Resolution diagnostics are measured every 0.5 baseline time units and 0.05 response units; CFL is tracked at every RK stage. [Forced isotropic turbulence](https://arxiv.org/abs/1406.6317) supplies the physical context. The [finite-resolution study](https://doi.org/10.1103/PhysRevFluids.3.064603) motivates explicit refinement rather than relying on a single cutoff diagnostic.

## Independent prediction results

The predicted target is the trapezoidal integral of difference energy divided by its initial value over the 1.5-unit response, with saved spacing 0.05 held fixed. RMSE is calculated over all eight targets within each flow, then averaged across the 32 independent flows. The global-current model has 10 inputs; richer current measurements have 30, including local velocity, all velocity-gradient components, local strain, vortex stretching, velocity derivative and initial difference-energy production/loss/forcing rates. Matched history adds 75 causal lag differences for 105 total inputs. Passive-marker history is not tested.

Training-only scaling and the same candidate linear/RBF ridge grid apply to all families. Validation chooses models without test labels. The source and protocol were registered before the new integrations. The immutable model artifact was frozen at **{freeze['created_utc']}**, before every test flow started; source, calibration and pilot hashes are verified at test launch. This chronology is checked again by the delivered audit.

| Comparator | Mean run RMSE | Error reduction with matched history and paired 95% interval |
| --- | --- | --- |
{chr(10).join(comparators)}

Three wrong-history models are trained with whole-run history derangements. Their test error averages every distinct run as a donor. The bootstrap resamples original run identities in both donor and recipient roles, excluding self donors even when resampling duplicates a run. Ten thousand whole-run bootstrap replicates quantify uncertainty. Support requires an advantage against all five comparators and passing all physical/numerical gates. The individual intervals are not simultaneous confidence intervals; uncertainty is conditional on this calibration set and model-selection procedure.

![Prediction errors and paired uncertainty](results/figures/prediction.png)

The disturbance-strength subgroup results below are descriptive and do not replace the declared joint criterion.

| Impulse ratio | History error reduction versus richer current and 95% interval |
| --- | --- |
{chr(10).join(subrows)}

## Physical mechanisms

The difference-energy equation is `dE_delta/dt = P-D+Q`: baseline strain produces `P`, viscosity dissipates `D`, and differing feedback-force fields supply `Q`. The solver integrates all terms at the same RK stages as velocity. The numerical residual checks whether these mechanisms account for the observed difference-energy response.

| Impulse ratio | Final difference-energy ratio and 95% interval | Integrated production / initial energy | Integrated viscous loss / initial energy | Integrated forcing work / initial energy |
| --- | --- | --- | --- | --- |
{chr(10).join(physical)}

![Disturbance amplification and energy budgets](results/figures/response_budgets.png)

Advection, strain, viscous loss and the specified forcing feedback explain the resolved response within the measured budget residual. Positive production transfers energy from the background flow. Trajectory separation alone is not an energy or causal measurement. The complete current field under this forcing law determines subsequent evolution; any advantage from partial-observation history would reveal useful missing information, not a new memory force. Failure of this feature set does not rule out other history measurements or targets.

## Reproduce the results

Use the [source repository]({CODE}) with Python 3.12. Restore archived observations before auditing their bytes:

```text
python -m pip install -r requirements.txt
python restore_recordings.py
python -m pytest test_solver.py test_design.py -q
python audit.py
```

In a separate copy, `forecast.py evaluate`, `analyze.py` and `report.py` rebuild predictions, tables and figures from recorded observations. Fresh integrations require a clean copy without earlier observations or frozen models; keep the delivered archive intact, then run `pipeline.py pilots`, `pipeline.py calibration` and `pipeline.py confirmation`. The public workflow uses those same scientific functions with archived stage dependencies. Regenerated floating-point bytes can differ across platforms and require a new manifest.

## Proposed experiments and unanswered questions

- Repeat the complete training, validation and test ensemble with new seeds to estimate variability from model selection, and refine every test flow.
- Increase Reynolds number further with larger resolved grids, longer stationary windows and independently measured scale separation. A broad inertial range remains a separate test.
- Compare repeated stationary ensembles across several viscosities with the same intervention design to estimate how prediction benefit changes with Reynolds number.
- Compare richer spatial current measurements, passive-marker geometry histories, noisy observations, other source positions and longer response horizons, using new frozen models and independent runs.
- Replay identical external forcing across copies to isolate fluid dynamics from forcing feedback, and compare other forcing laws.
- Match a physical fluid or laboratory dataset. Droplet evaporation, smoke particle effects, walls and free surfaces need additional models; ice requires separate solid mechanics and thermal physics.

These are future experiments, separate from the completed numerical results above.
'''
    (ROOT/'README.md').write_text(text,encoding='utf-8');(ROOT/'report.html').write_text(render(text),encoding='utf-8')
    save('results/completion.json',dict(completed_utc=now(),independent_baselines=64,independent_test_runs=32,test_targets=256,refinement_run_files=8,refinement_events=p['refinement_events'],numerical_pass=checks['all_passed'],calibration_stationarity_pass=cal['passed'],test_stationarity_pass=s['passed'],history_supported=bool(support),code_url=CODE,streams_url=STREAMS))
    print('Report generated from completed values')

if __name__=='__main__':main()
