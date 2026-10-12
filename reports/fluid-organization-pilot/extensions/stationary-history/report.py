"""Create the report from completed measurements, with future experiments separate."""
import csv,html,re
from common import ROOT,read,protocol,save,now
CODE='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-organization-pilot/extensions/stationary-history'
STREAMS='https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization/extensions/stationary-history'

def inline(s):
 s=html.escape(s);s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s);s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s);return re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
def render(text):
 out=[];lines=text.splitlines();i=0
 while i<len(lines):
  line=lines[i]
  if not line.strip():i+=1;continue
  if line.startswith('```'):
   code=[];i+=1
   while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
   out.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>');i+=1;continue
  if line.startswith('|'):
   rows=[]
   while i<len(lines) and lines[i].startswith('|'):rows.append([c.strip() for c in lines[i].strip('|').split('|')]);i+=1
   out.append('<div class="table"><table><thead><tr>'+''.join('<th>'+inline(c)+'</th>' for c in rows[0])+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+inline(c)+'</td>' for c in row)+'</tr>' for row in rows[2:])+'</tbody></table></div>');continue
  if line.startswith('!['):
   m=re.match(r'!\[([^]]*)\]\(([^)]+)\)',line);out.append('<figure><img alt="'+html.escape(m[1])+'" src="'+html.escape(m[2])+'"></figure>');i+=1;continue
  if line.startswith('#'):
   level=len(line)-len(line.lstrip('#'));out.append(f'<h{level}>'+inline(line[level:].strip())+f'</h{level}>');i+=1;continue
  if line.startswith('- '):
   items=[]
   while i<len(lines) and lines[i].startswith('- '):items.append('<li>'+inline(lines[i][2:])+'</li>');i+=1
   out.append('<ul>'+''.join(items)+'</ul>');continue
  para=[line];i+=1
  while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','|','![','- ','```')):para.append(lines[i]);i+=1
  out.append('<p>'+inline(' '.join(para))+'</p>')
 return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>History prediction in sustained turbulence</title><style>body{font:16px/1.65 system-ui,Segoe UI,sans-serif;color:#24323b;max-width:1120px;margin:auto;padding:32px}h1{line-height:1.2;font-size:2.15rem}h2{margin-top:2.1em;line-height:1.3}a{color:#086580}code{background:#eef3f6;padding:2px 4px}pre{background:#eef3f6;padding:16px;overflow:auto;font-size:13px}.table{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.91em}th,td{border-bottom:1px solid #d7e0e5;text-align:left;padding:8px 10px}th{background:#edf3f6}img{width:100%;height:auto}figure{margin:25px 0}</style><main>'+''.join(out)+'</main></html>'

def main():
 p=protocol();r=read('results/prediction_summary.json');stat=read('results/test_stationarity.json');cal=read('results/calibration_stationarity.json');phys=read('results/physical_summary.json');checks=read('results/numerical_checks.json');ref=read('results/refinement_prediction_summary.json');freeze=read('results/freeze_receipt.json')
 valid=checks['all_passed'] and stat['passed'] and cal['passed'];supported=r['history_supported'] and valid
 conclusion=('Matched history improved prediction beyond every declared current-state and null comparator.' if supported else 'This experiment did not establish a useful prediction advantage from matched history beyond all the declared controls.')
 if not valid:conclusion='The completed ensemble did not pass every numerical or stationarity criterion; its prediction comparison is exploratory. '+conclusion
 primary=r['contrasts']['present'];lo,hi=primary['reduction_interval'];relo,rehi=stat['re_lambda_interval']
 comparison=[]
 for name,label in [('present','Richer current state'),('global_present','Global current state'),('wrong_history','Wrong-run histories'),('training_prior','Training mean by strength'),('event_prior','Training mean by time and strength')]:
  q=r['contrasts'][name];value=q['null_mean_run_rmse'] if name=='wrong_history' else r['mean_run_rmse'][name];a,b=q['reduction_interval'];comparison.append(f"| {label} | {value:.6f} | {q['error_reduction_percent']:+.2f}% [{a:+.2f}, {b:+.2f}] |")
 with (ROOT/'results/prediction_by_strength.csv').open(encoding='utf-8') as fh:secondary=list(csv.DictReader(fh))
 subgroup=[f"| {float(q['strength']):.2f} | {float(q['present_mean_run_rmse']):.6f} | {float(q['history_mean_run_rmse']):.6f} | {float(q['reduction_percent']):+.2f}% [{float(q['ci_low']):+.2f}, {float(q['ci_high']):+.2f}] |" for q in secondary]
 station=[]
 for key,label in [('energy','Energy'),('dissipation','Dissipation'),('integral_length','Integral length')]:
  q=stat[key];a,b=q['drift_interval'];station.append(f"| {label} | {q['early_mean']:.5f} | {q['late_mean']:.5f} | {100*q['relative_drift']:+.2f}% [{100*a:+.2f}, {100*b:+.2f}] |")
 physical=[]
 for q in phys['responses']:physical.append(f"| {q['relative_impulse']:.2f} | {q['mean_final_energy_ratio']:.4f} [{q['final_ratio_ci_low']:.4f}, {q['final_ratio_ci_high']:.4f}] | {q['mean_production_integral']:.4f} | {q['mean_dissipation_integral']:.4f} | {q['mean_forcing_difference_integral']:.4f} |")
 grid=max(q['max_curve_relative_change'] for q in checks['refinements'] if q['label']=='fine');step=max(q['max_curve_relative_change'] for q in checks['refinements'] if q['label']=='half')
 registration=read('results/public_protocol_receipt.json') if (ROOT/'results/public_protocol_receipt.json').exists() else None
 if registration and registration.get('frozen_models_sha256')!=freeze['sha256']:registration=None
 registration_text=f"The [frozen source and models]({registration['url']}) were published before the held-out integrations." if registration else 'The local freeze receipt records the source, calibration recordings and model hashes before the held-out integrations.'
 execution=read('results/execution.json') if (ROOT/'results/execution.json').exists() else None
 execution_text=(f"Acquisition used {len(execution['calibration_local_seeds'])} local and {len(execution['calibration_hosted_seeds'])} hosted calibration fields, and {len(execution['test_hosted_seeds'])} hosted plus {len(execution['test_local_seeds'])} local test fields. The [execution record](results/execution.json) identifies platforms, seed groups and workflow commits. Cached files are counted only once. Original recording bytes and source hashes are retained." if execution else 'The execution environment and acquisition metadata accompany the source and recordings.')
 if execution and execution['local_20_test_pass'] and execution['linux_20_test_pass']:execution_text+=' The 20 solver/design tests passed locally and on the Linux runner; both receipts are retained.'
 screening_text=''
 if (ROOT/'results/long_refinement_comparisons.json').exists():
  pairs=read('results/long_refinement_comparisons.json');screen_rows=[]
  for pair in pairs:
   changes=pair['relative_changes'];screen_rows.append(f"| {pair['nu']:.3f} | {pair['coarse_n']} to {pair['fine_n']} | {pair['coarse_step']} to {pair['fine_step']} | {100*changes['energy']:+.3f}% | {100*changes['dissipation']:+.3f}% | {100*changes['flux_k4']:+.3f}% |")
  screening_text='''Separate long screening integrations compare the same initial Fourier field across grids. The table gives changes in time-window means from 20 to 45 for one seed per pair. These descriptive differences have no independent ensemble confidence interval. Changes in step and grid are combined where both settings differ; chaotic trajectories need not remain pointwise close over these long intervals.

| Viscosity | Grids | Time steps | Energy change | Dissipation change | Forward-flux change |
| --- | --- | --- | --- | --- | --- |
'''+chr(10).join(screen_rows)+'''

The [screening table](results/screening_regimes.csv) retains every attempted case. The completed 48-cubed pilot fell below the later resolution threshold and was excluded from confirmation. An 80-cubed, viscosity-0.02 pilot at step 0.025 was rejected by its initial-step CFL guard; this is a numerical rejection, not evidence of physical instability. The [source provenance](results/screening_provenance.json) also documents the first driver's source-hash and filename bookkeeping errors. Screening cases do not count toward the independent prediction sample.
'''
 text=f'''# History prediction in sustained three dimensional turbulence

**Completed result: {conclusion}** There are {r['independent_test_runs']} independent test runs, separate from 24 training and 8 validation runs. Each flow supplies four intervention times and two impulse strengths, giving 256 test targets while retaining only 32 independent sampling units. Matched-history mean run error is **{r['mean_run_rmse']['history']:.6f}**. Against the richer current-state predictor, the error reduction is **{primary['error_reduction_percent']:+.2f}%**, with a paired 95% interval **[{lo:+.2f}%, {hi:+.2f}%]**. Positive reductions favor history; the full control comparison appears below.

The measured approximate Taylor Reynolds number is **{stat['re_lambda_mean']:.2f} [95% interval {relo:.2f}, {rehi:.2f}]**, compared with about 19.75 in the [earlier transient study](../stress-history/README.md). This is a higher Reynolds number, sustained-flow experiment at a still modest Reynolds number. It does not establish a broad inertial range or high-Reynolds asymptotic behavior. Forcing, initialization and measurements also changed, so differences between studies cannot be attributed to Reynolds number alone.

[Runnable code, tests and restorable recordings]({CODE}) · [Streams and Rocks report]({STREAMS})

## Completed experiment and parameters

| Item | Executed setting |
| --- | --- |
| Fluid equations | Three dimensional, incompressible, constant density Navier–Stokes |
| Domain and boundaries | Periodic cube of side 6, dimensionless model units |
| Main resolution | {p['n']} cubed, strict two thirds Fourier cutoff, double precision |
| Viscosity and power input | Viscosity {p['nu']}; power per unit mass {p['power']} |
| Forcing | Constant-power negative damping in modes 0 less than mode length at most 2 |
| Integration | RK4; baseline step {p['spin_dt']}, response step {p['response_dt']} |
| Independent initial conditions | 64 fresh fields from Gaussian random draws, projected to be solenoidal and restricted to low modes; normalized RMS speed {p['initial_speed']} |
| Interventions | Copies at times 32, 36, 40 and 44; localized divergence-free impulses of 0.25 and 0.75 baseline RMS speed |
| Response horizon | 1.5 model-time units per intervention |
| Histories | Current observations and differences at lags 0.5, 1 and 2 |
| Refinement controls | Four predeclared test runs; half steps and an {p['fine_n']} cubed grid at the same smaller steps |

The disturbance is a fixed Fourier polynomial derived from the curl of a smooth periodic vector potential. Its strength is a global RMS velocity ratio, not a percentage increase in total kinetic energy. Every intervention starts from identical copies of the same baseline field. Both copies use the same power-control law evaluated on their respective states, so their forcing fields can differ.

Large-scale constant-power forcing follows an established [stationary-turbulence DNS approach](https://arxiv.org/abs/1406.6317). Constant energy alone is insufficient evidence for turbulence: [low-Reynolds flows under negative damping can self-organize](https://doi.org/10.1088/1751-8113/48/25/25FT01). Here the stationarity assessment also checks the energy spectrum, mean dissipation, turnover times and sustained nonlinear transfer to higher modes. This is an idealized forced fluid; material properties of water, air, fog, smoke or ice are not calibrated by this experiment.

{execution_text}

## Stationarity and numerical validation

The calibration stationarity check **{'passed' if cal['passed'] else 'failed'}** before model freezing. The independent test check **{'passed' if stat['passed'] else 'failed'}**. The first intervention follows at least **{stat['minimum_spin_turnovers']:.2f}** accumulated integral turnover times. Means in windows 20 to 32 and 32 to 44 are compared across whole independent runs; time samples are not treated as independent replicates.

| Quantity | Earlier window | Later window | Relative change and paired 95% interval |
| --- | --- | --- | --- |
{chr(10).join(station)}

The declared equivalence bound for these three relative changes is ±15%. The mean dissipation/input mismatch is {100*stat['power_balance']['relative_mismatch']:+.2f}%, with its 95% interval required to lie within ±10%. The total variation between normalized mean spectra is {100*stat['spectrum']['total_variation']:.2f}%, with a 95% upper bound of {100*stat['spectrum']['upper_95_percent_bound']:.2f}% against a 10% limit. Mean forward flux through mode length 4 is {stat['flux_mean']:.5f}; its 95% interval is [{stat['flux_interval'][0]:.5f}, {stat['flux_interval'][1]:.5f}], with the lower limit required to exceed 5% of input power. These finite-window criteria support approximate statistical stationarity; they cannot prove stationarity of every higher-order statistic for unlimited time.

Component mean-square velocities relative to their isotropic mean are {', '.join(f'{x:.3f}' for x in stat['component_ratio_to_isotropic'])}. The reported Taylor Reynolds number uses the isotropic approximation `sqrt(15) times (2 E / 3) / sqrt(viscosity times dissipation)`. The integral Reynolds number in the data tables is `sqrt(2 E) times integral length / viscosity`; its speed is the full vector RMS speed. Exact isotropy of every finite realization is not assumed.

![Stationarity diagnostics and spectra](results/figures/stationarity.png)

All declared numerical checks **{'passed' if checks['all_passed'] else 'did not pass'}**. Minimum sampled cutoff times Kolmogorov length is {phys['minimum_keta']:.4f}, against a 1.5 requirement. The cutoff is `(ceil(N/3) - 1) times 2 pi / 6`, using the inscribed axis cutoff rather than the larger corner wavenumber; the Kolmogorov length is `(viscosity cubed / dissipation) to the one-quarter power`. The largest saved normalized difference-energy balance residual is {phys['maximum_difference_budget_residual']:.3g}. Resolution and spectral-tail diagnostics are sampled every 0.5 baseline time units and 0.05 response time units; CFL is tracked at every RK stage. Tests cover an exact decaying solution, a forced equilibrium without a cascade, incompressibility, nonlinear energy conservation, constant power input, identical clones, the full difference-energy budget, matched grids, stationarity rejection controls, and protection against prediction leakage.

The largest relative change in a normalized response curve is {100*step:.5f}% for the half-step comparison and {100*grid:.5f}% for finer grid versus half-step at the same time step. Each control restarts from the same base-grid field two time units before an intervention, then recomputes history and response. These are dependent numerical controls, not additional independent evidence. They test short-time resolution of the prediction calculation; long-run screening comparisons address statistical grid behavior separately. They do not constitute continuum extrapolation of the entire ensemble.

![Resolution checks and unchanged frozen predictors](results/figures/resolution.png)

{screening_text}

## Independent prediction comparison

The target is the integral of `E_delta(t) / E_delta(0)` from time 0 to 1.5 after the impulse, evaluated by trapezoidal quadrature on saved samples spaced 0.05 time units apart. It measures disturbance persistence and amplification in model-time units. The error score first takes RMSE over all eight targets within one run, then averages those 32 run errors. The recording cadence is held fixed in the integration-step and grid comparisons.

The global-current predictor receives eight global measurements, impulse strength and intervention time. The richer-current predictor adds local velocity, all nine velocity-gradient components, local strain and vortex stretching, the local instantaneous velocity derivative, and the initial normalized production, dissipation and forcing-difference rates: 30 inputs in total. The history predictor adds 75 causal lag differences, for 105 inputs. These are Eulerian measurements, not a test of passive-marker history. Full feature names are recorded in the source.

All input families use the same candidate linear/RBF ridge models and validation grid. Scaling uses training data only, and validation selects models by mean run RMSE. All times and strengths from a run remain together in a split. The model freeze is recorded at **{freeze['created_utc']}**. {registration_text} No test targets are used for tuning. The protocol was developed during screening and calibration; the final source/model freeze is prospective with respect to the independent test, not the calibration observations. The earlier design registration is retained separately.

| Comparator | Comparator mean run RMSE | Error reduction with matched history and 95% interval |
| --- | --- | --- |
{chr(10).join(comparison)}

Wrong-history controls use three models trained with complete run histories permuted between independent runs. Their test error averages every different run as a history donor, keeping times and strengths together. The uncertainty calculation resamples run identities in both recipient and donor roles and excludes the same original run even if bootstrap sampling duplicates it. Donor pairs are not independent observations.

Support requires a benefit against every declared comparator and passing numerical and stationarity checks. This joint criterion prevents a weak current-state model, extra model inputs, gradual settling, or an unusually favorable history permutation from establishing the claim by itself. Individual 95% comparison intervals are not simultaneous confidence intervals. Uncertainty is conditional on this calibration ensemble, feature set and model-selection procedure; retraining on new calibration ensembles remains a separate validation question.

![Independent prediction errors and paired uncertainty](results/figures/prediction.png)

The strength-specific comparison below is descriptive. Its intervals do not replace the combined, predeclared control criterion and are not corrected for exploratory subgroup comparisons.

| Impulse ratio | Richer current RMSE | History RMSE | History error reduction and 95% interval |
| --- | --- | --- | --- |
{chr(10).join(subgroup)}

Across the four refinement runs, the largest change in the history-minus-current error contrast is {ref['max_absolute_step_contrast_change']:.6g} for time step and {ref['max_absolute_grid_contrast_change']:.6g} for grid spacing. The complete comparisons are in [the refinement table](results/refinement_predictions.csv).

## Physical response and mechanisms

The resolved difference-energy equation is `dE_delta/dt = P - D + Q`, where `P = -mean(delta_u_i delta_u_j d_j u_base_i)` is production by baseline strain, `D = viscosity times mean(abs(grad delta_u)^2)` is viscous loss, and `Q = mean(delta_u dot (f_disturbed - f_base))` is forcing-difference work. The code integrates all three terms with the same RK stages as the velocity field. In this experiment forcing feedback does not cancel.

| Impulse RMS ratio | Final energy ratio and 95% interval | Integrated production / initial energy | Integrated loss / initial energy | Integrated forcing difference / initial energy |
| --- | --- | --- | --- | --- |
{chr(10).join(physical)}

![Disturbance response and the complete energy budget](results/figures/response_budgets.png)

Nonlinear advection and strain, viscous damping, and the specified forcing feedback account for the recorded response within the checked numerical residuals. A positive difference-energy production term transfers energy from the background flow; it does not create energy without a source. A larger normalized separation of trajectories is not by itself an energy or causal measurement. This study does not identify weakening molecular interactions, and a change in prediction accuracy is not a causal intervention on organization.

The full current fluid field with this forcing law determines its subsequent evolution. Any benefit from a short history of partial measurements would be useful evidence about missing current information, not proof of an additional memory force in Navier–Stokes. Conversely, failure of these declared history features does not exclude benefits from other observations, targets, horizons or independently tested model families.

## Reproduce the completed results

Use the [code repository]({CODE}) with Python 3.12. Restore the archived observations and audit the original bytes before regenerating files:

```text
python restore_recordings.py
python -m pip install -r requirements.txt
python -m pytest test_solver.py test_design.py -q
python audit.py
```

In a separate working copy, `forecast.py evaluate`, `analyze.py`, `report.py`, `manifest.py` and `audit.py` rebuild the statistical results and figures from recordings. Fresh DNS runs are more costly. Preserve the delivered directory, and use another copy with the prior `data/runs`, `data/cloud`, `data/logs`, `data/predictions.npz` and entire `results` directory archived outside that copy. Keep the protocol and `data/screen` historical screening records. This removes stale models, public-registration receipts and execution records from the fresh copy. Run:

```text
python pipeline.py calibration --workers 8
python pipeline.py confirmation --workers 8
```

The calibration phase refuses to freeze models unless its numerical and stationarity checks pass. The confirmation phase verifies frozen source/model hashes before test runs and does not refit them. Run `environment.py` to record a reproduction environment. Full fields are regenerable from seeds; selected full restart fields, time series, observed targets, predictions, source hashes, failed screening records and exact archive checksums are retained. Regenerated floating-point bytes can differ across platforms; a new manifest describes that reproduction rather than the original snapshot.

## Remaining experiments and unanswered questions

- Repeat the whole training, validation and test procedure with new independent ensembles to measure model-selection variability, and refine every test run to extend the present four-run numerical check.
- Test substantially larger Reynolds numbers with a wider resolved cascade and longer observation windows. The present Reynolds number remains modest.
- Compare multiple stationary Reynolds-number regimes with the same forcing, observations and intervention design to isolate Reynolds-number effects.
- Test richer spatial current measurements, passive-marker geometry histories, noisy observations, different source positions and other response horizons, using new frozen models and independent runs.
- Vary the forcing law or replay identical external forcing across copies to isolate the contribution of forcing feedback from fluid dynamics.
- Calibrate a physical fluid or laboratory dataset. Fog evaporation and heat exchange, smoke particle effects, walls and free surfaces need their own models and checks; ice requires separate solid mechanics and thermal physics.

These are proposed follow-ups. They are not completed outcomes of the numerical experiment above.
'''
 (ROOT/'README.md').write_text(text,encoding='utf-8');(ROOT/'report.html').write_text(render(text),encoding='utf-8')
 save('results/completion.json',dict(completed_utc=now(),independent_baselines=sum(len(x) for x in p['splits'].values()),independent_test_runs=len(p['splits']['test']),test_targets=len(p['splits']['test'])*len(p['events'])*len(p['strengths']),main_interventions=sum(len(x) for x in p['splits'].values())*len(p['events']),refinement_run_files=2*len(p['refinement_seeds']),numerical_pass=checks['all_passed'],calibration_stationarity_pass=cal['passed'],test_stationarity_pass=stat['passed'],history_supported=bool(supported),code_url=CODE,streams_url=STREAMS))
 print('Report generated from completed measurements.')
if __name__=='__main__':main()
