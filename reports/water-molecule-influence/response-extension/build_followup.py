"""Generate a concise HTML report and GitHub README from the measured follow-up."""
from pathlib import Path
import base64,hashlib,json

R=Path(__file__).resolve().parent;D=R/'data'
results=json.loads((D/'results.json').read_text())
controls=json.loads((D/'controls.json').read_text())
analysis=json.loads((D/'analysis.json').read_text())
protocol=json.loads((D/'protocol.json').read_text())
rows=results['rows']
times=[.005,.02,.1,.2]

def at(t):return [r for r in rows if abs(r['time_ps']-t)<1e-10]
def triplet(values,format_spec='.3f'):return ', '.join(format(x,format_spec) for x in values)
def table(headers,data):
    return '<table><tr>'+''.join('<th>'+str(x)+'</th>' for x in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in data)+'</table>'

summary_table=table(['Time after impulse','Force/response rank correlations (runs 1–3)','Original force leader’s rank (of 216)'],[
    [f'{t*1000:g} fs',triplet([r['spearman_force_vs_response'] for r in at(t)]),
     ', '.join(str(r['initial_force_leader_response_rank']) for r in at(t))] for t in times])
detail_table=table(['Run','Initial force leader','Response leaders at 5, 20, 100, 200 fs','Response max/min ratio at 200 fs'],[
    [rep+1,at(.005)[rep]['initial_force_leader'],', '.join(str(at(t)[rep]['response_leader']) for t in times),
     f"{at(.2)[rep]['response_max_min_ratio']:.2f}"] for rep in range(3)])
half_kick=max(c['half_kick']['relative_L2_difference'] for c in controls)
half_kick_single=max(c['half_kick']['max_relative_difference'] for c in controls)
half_step=max(c['half_timestep']['relative_L2_difference'] for c in controls)
free=max(c['no_interactions_max_score_ps'] for c in controls)
replay=max(c['same_initial_state_replay_max_COM_difference_nm'] for c in controls)
energy=max(c['max_energy_drift_kJ_mol_per_water'] for c in controls)
early_half_step=max(d['half_timestep_by_time']['relative_L2'][0] for d in analysis['diagnostics'])
late_half_step=max(d['half_timestep_by_time']['relative_L2'][-1] for d in analysis['diagnostics'])
late_half_step_single=max(d['half_timestep_by_time']['maximum_relative'][-1] for d in analysis['diagnostics'])
early_rho=triplet([r['spearman_force_vs_response'] for r in at(.005)])
late_rho=triplet([r['spearman_force_vs_response'] for r in at(.2)])
late_ranks=', '.join(str(r['initial_force_leader_response_rank']) for r in at(.2))
image='data:image/png;base64,'+base64.b64encode((R/'finite_time_results.png').read_bytes()).decode()
html=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Finite-time water influence: follow-up experiment</title>
<style>body{{font:17px/1.6 system-ui,sans-serif;max-width:1000px;margin:48px auto;padding:0 24px;background:#fafbfc;color:#172630}}h1{{font-size:35px;line-height:1.2}}h2{{font-size:23px;margin-top:32px}}img{{width:100%}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:20px 0}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #d4dde2}}th{{background:#e8eff2}}.lead{{font-size:21px}}.box{{background:#edf5f5;border-left:4px solid #087e83;padding:14px 20px}}code,pre{{background:#edf1f3}}pre{{padding:15px;overflow:auto}}a{{color:#087e83}}li{{margin:8px 0}}</style>
<p>FOLLOW-UP TO THE IDENTICAL-WATER-MOLECULE EXPERIMENT</p>
<h1>A snapshot leader can lose its lead as influence spreads</h1>
<p class="lead">The original force ranking predicts the response to a tiny impulse extremely well at 5 femtoseconds. By 200 femtoseconds its rank correlation is {late_rho} in the three tested states. The original highest-ranked molecules have fallen to ranks {late_ranks} out of 216.</p>
<p>This follow-up directly changes molecular velocities and follows the resulting motion. It strengthens the evidence for roles that depend on the current state and observation time. It does not establish a lasting hierarchy, nor does it validate this water model against physical experiments.</p>
<img src="{image}" alt="Rank correlations, leader ranks, response scatter, and response dispersion at four observation times">
{summary_table}
{detail_table}
<p>A femtosecond (fs) is 10⁻¹⁵ seconds; 200 fs is 0.2 ps. All molecule IDs are arbitrary zero-based labels. Correlations compare 216 molecules within each state; these are not 216 independent experimental replicates.</p>
<h2>What changed from the original experiment</h2>
<p>The earlier water score asked how a tiny <em>displacement</em> of one molecule changes other molecules’ forces at the same instant. Here the question is whether that ranking predicts motion after a tiny <em>velocity impulse</em>. The last saved configuration (50 ps production time) from each of the original three TIP3P runs was used. Each was assigned new constrained thermal velocities at 300 K using seeds 20261011–20261013, with overall center-of-mass velocity removed.</p>
<p>Every molecule was tested as a source. For each source, positive and negative impulses were applied separately along all three Cartesian axes. All three atoms of the source received the same velocity change, ±0.01 nm/ps (10 m/s). Each other molecule received the opposite change divided by 215, so the added total momentum was zero. The average kinetic-energy increase of the positive/negative pair is the same for every source; individual signed trials can have different changes because of their initial velocities.</p>
<p>These are six trajectories per source, 216 sources per state and three states: <strong>3,888 primary perturbed trajectories</strong>. Each was followed for 0.2 ps and sampled at 5, 20, 100 and 200 fs. Additional trajectories provided convergence, replay, relabeling and zero-interaction controls. The physical model and periodic box match the original experiment.</p>
<h2>How influence was measured</h2>
<p>For source i, target j, and kick axis b, the paired position difference gives A<sub>ji</sub> = [R<sub>j</sub>(+δv) − R<sub>j</sub>(−δv)]/(2δv), where R is the molecular center of mass. Each A is a 3 × 3 directional response matrix. The response caused directly by the compensating background kick is known and removed:</p>
<div class="box">B<sub>ji</sub>(t) = (N−1)/N × [A<sub>ji</sub>(t) + tI/(N−1)]<br>S<sub>i</sub>(t) = Σ<sub>j ≠ i</sub> ‖B<sub>ji</sub>(t)‖<sub>F</sub>.</div>
<p>The rescaling makes B the central response equivalent to a single-source impulse, using the invariance of internal motion to a uniform velocity boost. Excluding the source prevents its direct ballistic motion from inflating the score. Summing Frobenius norms combines all axes and target molecules. The score’s units are ps (position divided by velocity); it is not an energy fraction, and a sum over many targets can exceed the elapsed time.</p>
<p>At sufficiently short times the equations of motion predict B<sub>ji</sub>(t) ≈ t³K<sub>ji</sub>/(6M), where K is the original translational force Jacobian and M is molecular mass. This explains why the original score should predict the early response. The observed correlations at 5 fs were {early_rho}. As the interacting configuration evolves, a fixed initial Jacobian becomes a less complete predictor.</p>
<h2>Numerical method and controls</h2>
<p>The probes used deterministic NVE dynamics, double precision, 1 fs time steps and a constrained velocity-Verlet/RATTLE integrator with 10⁻⁸ constraint tolerance. Positions and velocities refer to the same time, which is necessary for meaningful comparisons across time steps. The thermostat and automatic center-of-mass remover were absent. OpenMM documents this <a href="https://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html">constrained velocity-Verlet implementation</a>. Input positions were projected to the tighter constraint tolerance; the tiny correction is recorded for each state.</p>
<ul><li><strong>Half impulse:</strong> ten sources per state, spaced across the initial force ranks. The largest aggregate relative L2 score difference was {100*half_kick:.4f}%. The largest individual score difference was {100*half_kick_single:.2f}% at the earliest 5 fs sample; aggregate errors are dominated by the larger later-time responses.</li>
<li><strong>Half time step:</strong> the weakest, middle and strongest force-score sources in each state. The largest aggregate relative L2 difference was {100*half_step:.3f}%. At 200 fs the largest across-state L2 difference was {100*late_half_step:.3f}%, and the largest individual difference was {100*late_half_step_single:.3f}%. At the very earliest 5 fs sample, the maximum L2 difference was {100*early_half_step:.2f}%; five 1-fs steps leave a noticeable amplitude-discretization error even when the ordering is stable.</li>
<li><strong>No interactions:</strong> all forces were removed while rigid constraints were retained. The largest corrected response score was {free:.2e} ps, verifying removal of the direct compensating motion.</li>
<li><strong>Replay:</strong> repeating the same unperturbed state produced a maximum center-of-mass discrepancy of {replay:.2e} nm.</li>
<li><strong>Relabeling:</strong> three checked sources in state 1 were tested after complete position/velocity blocks were reassigned. The response agreed with the predicted reassignment to relative L2 error {controls[0]['permutation']['relative_L2_difference']:.2e}.</li>
<li><strong>Energy:</strong> the largest unperturbed energy change over the sampled interval was {energy:.4g} kJ/mol per molecule. This short-run diagnostic is not a long-time stability guarantee.</li></ul>
<p>An exploratory leapfrog calculation was discarded after its half-step velocity convention made the initial physical states differ between time-step settings. Switching to on-step velocities reduced the checked aggregate discrepancy from approximately 4.6% to 0.06% in a preliminary comparison. All results in this report use the corrected method. The official <a href="https://docs.openmm.org/latest/api-python/generated/openmm.openmm.VerletIntegrator.html">VerletIntegrator documentation</a> identifies the built-in method as leapfrog.</p>
<h2>What the results support—and what remains open</h2>
<p>The data support a specific conclusion: <strong>an instantaneous force-based ranking is useful at very short times, but its leading identity and predictive value can change rapidly as motion develops</strong>. In these three states, the initial force leader was not the 200 fs response leader. The response still differed between molecules, so loss of the original ranking does not imply equal influence at that later time.</p>
<p>The instantaneous potential remains reciprocal. An additional saved diagnostic measures the asymmetry of the finite-time source/target response-norm matrix; forward responses on an evolving configuration need not retain the same symmetry as the instantaneous force Jacobian. This exploratory diagnostic is not evidence for permanent one-way control.</p>
<ul><li>Only three position/velocity microstates were sampled. These are conditional responses from the original short trajectories, not independent long-time equilibrium estimates.</li>
<li>The longer-time response depends on initial velocities as well as positions. This follow-up does not isolate oxygen-center position from molecular orientation or velocity.</li>
<li>Impulse and time-step controls cover selected sources, not every molecule. Full-ranking convergence at smaller time steps remains to be established.</li>
<li>The physical model is still rigid, nonpolarizable TIP3P with the original finite cutoff and PME settings. No alternate model, box-size study, nuclear quantum treatment or chemical reactions were added.</li>
<li>No significance test treats molecules or neighboring time samples as independent. The three-state values are reported individually rather than converted into an overconfident confidence interval.</li></ul>
<p>The next priorities are independent thermal velocity draws at fixed positions, all-source smaller-time-step comparisons, additional independently equilibrated configurations, larger boxes, alternate water models and longer response horizons with renewed linearity checks.</p>
<h2>Reproduce</h2>
<p>From the parent experiment directory, with its pinned dependencies installed:</p>
<pre>python response-extension/impulse_response.py --platform OpenCL
python response-extension/analyze_response.py
python response-extension/build_followup.py</pre>
<p>The OpenCL route requires double-precision support. Use <code>--platform Reference</code> for a portable but slower calculation. To preserve the delivered measurements during a new run, add <code>--out reproduced-response</code>. The plotting/report scripts read the delivered <code>response-extension/data</code> by default.</p>
<p>Saved arrays include initial atomic positions and velocities, baseline center-of-mass trajectories and energies, all source/target response norms at all four times, source scores, and convergence-control arrays. Source reproducibility and saved-data analysis are provided; bitwise agreement across hardware is not guaranteed. The parent <a href="../report.html">original report</a> documents the initial simulation protocol and its limitations.</p>
</html>'''
(R/'report.html').write_text(html,encoding='utf-8')
readme=f'''# Finite-time response: does the snapshot leader stay influential?

This follow-up tests all 216 molecules in each of three saved water states with paired, momentum-balanced velocity impulses. It follows the motion of the other molecules for 5–200 femtoseconds.

**Finding:** the original force ranking predicts the very early response, but loses predictive value over the measured interval. The original force leaders rank **{late_ranks} out of 216** by 200 fs. This supports temporary, observation-time-dependent roles, not a permanent hierarchy.

![Finite-time response results](finite_time_results.png)

| Time | Force/response Spearman correlations, runs 1–3 | Original force leader ranks |
| --- | --- | --- |
'''
for t in times:
    readme+=f"| {t*1000:g} fs | {triplet([r['spearman_force_vs_response'] for r in at(t)])} | {', '.join(str(r['initial_force_leader_response_rank']) for r in at(t))} |\n"
readme+=f'''
## Methods and checks

- 3,888 primary perturbed trajectories: 216 sources × 3 axes × 2 signs × 3 microstates.
- Rigid TIP3P, short NVE trajectories, double precision, 1 fs steps, constrained velocity Verlet with velocities and positions at the same time.
- Source kick: ±0.01 nm/ps, with compensating kicks to keep total added momentum zero. The direct compensating ballistic motion is removed analytically.
- Score: sum of Frobenius norms of other molecules' center-of-mass velocity-to-position responses; units are ps, not energy fraction.
- Half-kick controls: largest aggregate relative L2 difference {100*half_kick:.4f}% on ten selected sources per state; largest individual change {100*half_kick_single:.2f}% at 5 fs.
- Half-timestep controls: largest aggregate relative L2 difference {100*half_step:.3f}% on three selected sources per state. The earliest 5 fs amplitude has up to {100*early_half_step:.2f}% L2 discretization sensitivity.
- With interactions removed, the largest corrected score was {free:.2e} ps.
- Complete-state relabeling agreed to relative L2 error {controls[0]['permutation']['relative_L2_difference']:.2e} for the three checked sources.

Read the [full report](report.html) for definitions, numerical controls, energy checks, limitations and next validation steps. Download the HTML file and open it in a browser; GitHub displays its source.

## Reproduce

From the parent experiment directory, after installing its `requirements.txt`:

```sh
python response-extension/impulse_response.py --platform OpenCL
python response-extension/analyze_response.py
python response-extension/build_followup.py
```

Use `--platform Reference` for portable, slower execution. Use `--out reproduced-response` to keep a new run separate from the delivered data. Plotting and report generation use `response-extension/data` by default.

See [protocol](data/protocol.json), [numeric results](data/results.csv), [controls](data/controls.json), and [analysis diagnostics](data/analysis.json). The saved response arrays contain initial states and all source/target norm matrices, so the reported scores and figures can be independently reanalyzed.

## Limits

Three conditional microstates are not a long-time ensemble. Finite-time influence depends on initial velocities as well as positions. Convergence controls cover selected sources; full-ranking convergence, larger boxes, other water models, more thermal velocity draws and longer horizons remain to be tested. Molecules and time samples are correlated; no independent-sample p-value is claimed.

An initial exploratory leapfrog run was replaced because its time-step-dependent velocity convention did not provide matching physical starting states. The delivered results use corrected on-step velocities. Implementation follows the [OpenMM constrained velocity-Verlet documentation](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html).

[Return to the original experiment](../README.md).
'''
(R/'README.md').write_text(readme,encoding='utf-8')
manifest={p.relative_to(R).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted(R.rglob('*')) if p.is_file() and p.suffix!='.pyc' and p.name!='manifest_sha256.json'}
(R/'manifest_sha256.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Follow-up report, README and manifest saved.')
