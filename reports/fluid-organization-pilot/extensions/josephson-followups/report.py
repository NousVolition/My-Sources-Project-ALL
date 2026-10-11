"""Render the completed numerical evidence without converting proposals into results."""
import json,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from support import ROOT,save
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.dpi':160})
CODE='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-organization-pilot/extensions/josephson-followups'
STREAMS='https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization/extensions/josephson-followups'
def load(name):return json.loads((ROOT/f'results/{name}.json').read_text())
def figsave(fig,name):
    fig.savefig(ROOT/f'figures/{name}.png',bbox_inches='tight');fig.savefig(ROOT/f'figures/{name}.svg',bbox_inches='tight');plt.close(fig)
def ci(x):return f"{x['mean']:.5g} [{x['ci95'][0]:.5g}, {x['ci95'][1]:.5g}]"
def table(headers,rows):return '<table><thead><tr>'+''.join('<th>'+str(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'
def main():
    long=load('long_cycles');ramps=load('ramps');rr=load('ramps_refined');noise=load('noise');cap=load('capacitance');cp=load('circuit_prediction');fp=load('fluid_prediction')
    audit_count=load('verification')['total'] if (ROOT/'results/verification.json').exists() else 0
    fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    for a,ax in zip([0,.5,2],axes):
        for j,b in enumerate([1.005,1.02,1.2]):
            rows=[r for r in long['regimes'] if r['alpha']==a and r['bias']==b];finite=rows[:-1];whole=rows[-1]
            x=[r['window'] for r in finite];y=[r['mean_abs_up_down_gap'] for r in finite]
            line=ax.plot(x,y,'o-',label=f'i={b:g}')[0];ax.fill_between(x,[r['ci_low'] for r in finite],[r['ci_high'] for r in finite],color=line.get_color(),alpha=.12)
            ax.axhline(whole['mean_abs_up_down_gap'],color=line.get_color(),ls=':',alpha=.8)
        ax.set(xscale='log',xlabel='Averaging duration',title=f'Load ratio α = {a:g}',ylabel='Mean |up − down| voltage');ax.grid(alpha=.15)
    axes[0].legend(fontsize=9);fig.suptitle('Finite windows approach orbit-dependent whole-cycle means (dotted)');figsave(fig,'whole_cycles')
    z=np.load(ROOT/'data/ramps.npz');fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for leg in [200,800]:
        rows=[r for r in ramps['regimes']+rr['regimes'] if r['leg_duration']==leg]
        axes[0].loglog([r['current_increment'] for r in rows],[r['mean_run_rms_vs_smooth'] for r in rows],'o-',label=f'Leg duration {leg}')
    axes[0].set(xlabel='Staircase current increment',ylabel='Run RMS voltage difference from smooth ramp',title='Refinement at fixed elapsed time');axes[0].legend();axes[0].grid(alpha=.2)
    from matplotlib.ticker import FixedLocator,FixedFormatter,NullFormatter
    axes[0].xaxis.set_major_locator(FixedLocator([.0015625,.00625,.025,.05]));axes[0].xaxis.set_major_formatter(FixedFormatter(['0.00156','0.00625','0.025','0.05']));axes[0].xaxis.set_minor_formatter(NullFormatter())
    t=z['time_800'];b=z['bias_800'];q=z['smooth_fine_800'];v=2*np.diff(q[:,:,0],axis=0).mean(1)/np.diff(t);mid=(b[:-1]+b[1:])/2
    axes[1].plot(mid[:32],v[:32],label='Smooth up');axes[1].plot(mid[32:],v[32:],label='Smooth down')
    axes[1].set(xlabel='Current i',ylabel='Bin-averaged total voltage',title='A dynamic loop, not a static hysteresis test');axes[1].legend();figsave(fig,'ramps')
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for ax,b in zip(axes,[.99,1.005,1.02]):
        rows=[r for r in noise['regimes'] if r['bias']==b];y=np.array([r['mean_total_voltage'] for r in rows]);lo=np.array([r['ci_low'] for r in rows]);hi=np.array([r['ci_high'] for r in rows])
        ax.errorbar(range(5),y,yerr=[y-lo,hi-y],fmt='o',capsize=3);ax.axhline(rows[0]['matched_identical_noiseless_mean'],ls=':',color='gray',label='Identical / noiseless')
        ax.set(xticks=range(5),xticklabels=['0 / .02','1e−5 / 0','1e−5 / .02','1e−3 / 0','1e−3 / .02'],xlabel='Noise Γ / mismatch δ',ylabel='Mean total voltage, time 100–400',title=f'i = {b:g}');ax.tick_params(axis='x',labelrotation=35)
    axes[0].legend(fontsize=9);fig.suptitle('Dimensionless sensitivity: 24 starting draws and paired noise controls');figsave(fig,'noise_mismatch')
    q=np.load(ROOT/'data/capacitance.npz');fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout='constrained');n=16
    for j,ax in enumerate(axes):
        for name,ls in [('base','-'),('double_settle',':')]:
            v=q[name+'_voltage'][:,:, :].mean(-1)
            ax.plot(q['biases'],v[:n,j],color='#277da8',ls=ls,label='Up, S='+('100' if name=='base' else '200'))
            ax.plot(q['biases'],v[n:,j][::-1],color='#bc4b51',ls=ls,label='Down, S='+('100' if name=='base' else '200'))
        ax.set(xlabel='Current i',ylabel='Mean total voltage',title=f'Capacitance β = {q["betas"][j]:g}')
    axes[0].legend(fontsize=8);fig.suptitle('Capacitive pair on the exactly synchronized branch');figsave(fig,'capacitance')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,d,title in zip(axes,[cp,fp],['Circuit: 16 held-out starts','Fluid: 12 new independent fields']):
        names=['present','history','shuffled_history'];y=np.array([d[k]['mean'] for k in names]);lo=np.array([d[k]['ci95'][0] for k in names]);hi=np.array([d[k]['ci95'][1] for k in names])
        ax.errorbar(range(3),y,yerr=[y-lo,hi-y],fmt='o',capsize=5);ax.set(xticks=range(3),xticklabels=['Present','Present + history','Shuffled history'],ylabel='Mean per-run forecast RMSE',title=title);ax.tick_params(axis='x',labelrotation=15)
    fig.suptitle('Initial pilots: history benefit was unresolved against present measurements');figsave(fig,'held_out_prediction')
    # Structured headline values are also machine-readable.
    summary={'date':'2026-10-10','scope':'Completed dimensionless numerical pilots; no laboratory validation','whole_cycle_regimes':[r for r in long['regimes'] if r['window']=='whole_cycle'],'circuit_history_minus_present':cp['history_minus_present'],'fluid_history_minus_present':fp['history_minus_present'],'fluid_history_minus_shuffled':fp['history_minus_shuffled_history'],'remaining':['Physical calibration and laboratory validation','LRC load and general off-synchrony capacitive dynamics','Joint noise/mismatch/capacitance/ramp studies','More independent test starts and locked external validation','Physical transfer to fluid mechanisms, gas/fog/smoke/ice extensions']}
    save('results/summary.json',summary)
    whole_table=table(['α','i','Whole-cycle mean |up − down|','95% paired bootstrap interval'],[[r['alpha'],r['bias'],f"{r['mean_abs_up_down_gap']:.6g}",f"[{r['ci_low']:.6g}, {r['ci_high']:.6g}]"] for r in summary['whole_cycle_regimes']])
    pred_table=table(['Prediction study','Present measurements','Present + history','Shuffled history','Paired history − present'],[[title,ci(d['present']),ci(d['history']),ci(d['shuffled_history']),ci(d['history_minus_present'])] for title,d in [('Circuit',cp),('Fluid',fp)]])
    validate_table=table(['Check','Observed value','Interpretation'],[
        ['Cycle log-separation return',f"{long['checks']['max_log_separation_return_error']:.3g}",'Below 1e−5 gate'],
        ['Direct time integration versus cycle continuation',f"{long['checks']['max_direct_time_voltage_error']:.3g}",'Six trajectories to time 12,800; below 2e−6'],
        ['Cycle energy budget',f"{long['checks']['max_normalized_cycle_energy_error']:.3g}",'Input = change in potential + dissipation'],
        ['Smooth-ramp half-step voltage change',f"{ramps['max_smooth_step_refinement_voltage_error']:.3g}",'Below 2e−6'],
        ['Noisy paired step-refinement RMS',f"{noise['max_resolution_rms']:.3g}",'Below 0.005; coupled Brownian increments'],
        ['Free-diffusion noise covariance relative error',f"{noise['covariance_relative_error']:.3g}",'Below 0.02; independent 200,000-draw control'],
        ['Noisy potential-budget residual',f"{noise['max_normalized_potential_budget_residual']:.3g}",'Stratonovich drift + thermal work; diagnostic'],
        ['Capacitance half-step voltage change',f"{cap['max_step_voltage_change']:.3g}",'Below 2e−5'],
        ['Capacitance energy budget',f"{max(r['max_normalized_balance_error'] for r in cap['variants'].values()):.3g}",'Includes capacitive stored energy; below 1e−6'],
        ['New-fluid maximum energy budget / divergence',f"{fp['validation']['max_energy_budget_error']:.3g} / {fp['validation']['max_divergence']:.3g}",'Unchanged solver; independent new fields'],
        ['New-fluid largest time/grid energy difference',f"{max(max(r['time_curve_error_over_initial'],r['grid_curve_error_over_initial']) for r in fp['validation']['refinement']):.3g}",'Two new fields; below 0.01 of initial energy']])
    body=f'''<header><p class="eyebrow">STREAMS AND ROCKS · COMPLETED NUMERICAL FOLLOW-UPS · 10 OCTOBER 2026</p><h1>Does the past improve prediction?</h1><p class="lead">Some circuit history effects survive whole-cycle averaging. Extra predictive value from recorded history remains unresolved against present measurements in these small held-out tests.</p><p><a href="{STREAMS}">Streams and Rocks report</a> · <a href="{CODE}">Runnable code, tests and recorded arrays</a> · <a href="protocol.json">Declared protocol</a></p></header>
<section><h2>What was actually completed</h2><p>These are dimensionless simulations, not measurements of a fabricated circuit. The fluid study is a separate Navier–Stokes forecast replication.</p>
{table(['Study','Completed numerical work'],[['Longer observation / whole cycles','324 checkpoint cases, two cycles each; windows 200, 800, 3,200 and 12,800 via verified periodic continuation; six direct integrations to 12,800. Checkpoints reuse the prior 12 phase-pair draws.'],['Sweep protocols','12 new initial phase pairs, two ramp speeds, six staircase increments, and smooth ramps with a half-step control. Three finer increments were added after the original error was nonmonotone.'],['Noise / mismatch','360 condition-trajectories at each of two coupled time steps; 24 starting draws, three currents, five noise/mismatch combinations, plus 72 matched identical/noiseless controls.'],['Capacitance','Three capacitances, six equal-phase starts, up/down currents 0–1.5 at spacing 0.1, three numerical/settling variants: 1,728 plateaus. Exactly synchronized branch only.'],['Circuit forecasts','24 training, eight validation and 16 independent test phase pairs; three observation times and two interventions per run.'],['Fluid forecasts','Old 18/6 training/validation fields; 12 newly simulated test fields, each cloned into three forcing arms. Two new fields also checked at half the time step and a finer grid.']])}
<p>Repeated windows, current values, arms and numerical refinements are paired measurements, not extra independent samples. The report retains the original sweep’s numerical failures in its <a href="../josephson-sweep/README.md">separate earlier record</a>.</p></section>
<section><h2>1. Longer windows reveal different recurring orbits</h2><p>The unchanged overdamped model is M φ̇ = i1 − sin φ, with M = I + α11ᵀ, α = r/R, time τ = 2eI<sub>c</sub>rt/ℏ, and total voltage in units I<sub>c</sub>r. For identical junctions define p = (φ₁+φ₂)/2, q = (φ₁−φ₂)/2, z = log|tan(q/2)| and k = 1+2α. Then ṗ = (i + sin p tanh z)/k and ż = −cos p.</p>
<p>At the tested currents i &gt; 1, p increases monotonically. Integrating with p as the independent variable avoids the bottleneck near threshold. The phase/log-separation returns were checked over two full turns; six direct time integrations independently check the long-window continuation. At α=0, whole-cycle voltage agrees with 2√(i²−1), independently of phase separation. Coupled cases retain up/down differences:</p>{whole_table}
<figure><img src="figures/whole_cycles.png"><figcaption>Shading: 95% paired bootstrap intervals across the 12 inherited starts. Dotted levels: complete-cycle averages. Finite-window approach need not be monotone.</figcaption></figure>
<p>A direct identity explains the ideal-model persistence. Set C(z)=cosh(z)<sup>1/k</sup>. The quantity H(p,z)=sin(p)C(z)+(i/k)∫₀ᶻC(u)du is constant at fixed i. Differentiating gives dH/dτ=0. At fixed p, ∂H/∂z=C(z)(i+sin p tanh z)/k &gt; 0 for i&gt;1; hence returning to p modulo 2π returns to the same z. Each orbit has a whole-cycle mean 4π/T. A prior sweep can select different orbit labels through its current phase state.</p>
<p>This is a property of the specified ideal equations, not evidence for a new force or hidden memory beyond the complete state. It is not attracting-bistability hysteresis. Constants of motion in Josephson arrays are established research context (<a href="https://www.stevenstrogatz.com/articles/constants-of-motion-for-superconducting-josephson-arrays">Watanabe–Strogatz</a>; <a href="https://arxiv.org/abs/0812.4481">Marvel–Strogatz</a>). The large-N array result is not used as proof for this two-junction pair; the identity above is derived directly. The asymptotic conclusion here is restricted to identical, noiseless, overdamped junctions and sampled running currents; it does not resolve threshold i=1.</p></section>
<section><h2>2. Staircases converge to ramps, sometimes nonmonotonically</h2><p>The current follows 0.9 → 1.3 → 0.9. Each leg lasts either 200 or 800 time units. Midpoint staircases use increments 0.05, 0.025, 0.0125, then 0.00625, 0.003125 and 0.0015625. Total elapsed time and the 32 voltage bins per leg are identical. The last three increments are labeled as a follow-up to the nonmonotone original long-leg error.</p><p>At the finest increment, RMS differences from a smooth ramp are {rr['regimes'][2]['mean_run_rms_vs_smooth']:.4g} for the 200-unit leg and {rr['regimes'][5]['mean_run_rms_vs_smooth']:.4g} for the 800-unit leg. A narrower increment is not guaranteed to reduce every coarse-grid error: accumulated phase and the location of the final averaging boundary matter. The 100-unit subcritical burn makes these 12 starts nearly indistinguishable; the very narrow bootstrap intervals reflect that collapse, not broad physical certainty.</p><figure><img src="figures/ramps.png"><figcaption>Original and subsequent refinement measurements are both shown. The loop on the right is measured while current changes; it is not by itself evidence of static hysteresis.</figcaption></figure></section>
<section><h2>3. Dimensionless noise and mismatch alter the response</h2><p>Critical currents are c₁=1−δ and c₂=1+δ. The overdamped stochastic model is dφ=M⁻¹(i1−c sin φ)dτ + √(2Γ) chol(M⁻¹)dW, with all resistors assigned the same dimensionless temperature Γ=k<sub>B</sub>T/E<sub>J</sub>. Thus phase-increment covariance is 2ΓM⁻¹dτ; the common load’s thermal contribution is included. Additive stochastic Heun uses steps 0.01 and 0.02 with exactly summed Brownian increments. Conditions share random drivers for paired sensitivity comparisons; independent-run uncertainty uses 24 groups.</p><p>At i=0.99, the matched identical/noiseless mean voltage is about 1.5×10⁻⁶. A 2% critical-current mismatch yields 0.08426; Γ=0.001 with no mismatch yields 0.01999 [0.01453, 0.02542]. These are finite-window observations, not calibrated laboratory temperatures or proof that noise erases asymptotic history. Voltage is measured by phase increments from time 100 to 400; an instantaneous derivative of a white-noise path is not used.</p><figure><img src="figures/noise_mismatch.png"><figcaption>Error bars: bootstrap 95% intervals over independent starting/noise draws. Small ideal-model phase differences and noise are distinct mechanisms.</figcaption></figure><p>The normalization follows the standard classical RCSJ thermal-noise convention; see the dimensionless equation and noise covariance in <a href="https://link.aps.org/accepted/10.1103/PhysRevB.85.224507">Brenner et al., Appendix A</a>. Its experimental nanowire results are not calibration data for this simulated pair.</p></section>
<section><h2>4. Capacitance introduces stored energy and bistability</h2><p>For equal phases and velocities, the exact invariant branch obeys βp̈+(1+2α)ṗ+sin p=i, where β=2eI<sub>c</sub>r²C/ℏ. The pair energy is E=βṗ²−2cos p, with input power 2iṗ and dissipation 2(1+2α)ṗ². Tests use α=0.5 and β=0.2, 2 and 10. Half-step and doubled-settling sweeps were both completed.</p><p>For β=10 at i=0.9, upward continuation remains essentially at zero voltage while downward continuation remains running at about 0.789. The difference remains when settling doubles. This is consistent with ordinary inertial switching/retrapping in the capacitive model. Smaller finite-window voltage differences above threshold can still depend on averaging phase. Six common-phase starts converge onto the same branch, so their tiny spread is not a robustness test to arbitrary phase mismatch.</p><figure><img src="figures/capacitance.png"><figcaption>Solid: settling 100; dotted: settling 200. Voltage window 100 in both. The half-step curves agree closely and are recorded in the data. General off-synchrony stability and LRC loads were not simulated.</figcaption></figure></section>
<section><h2>5. Held-out history prediction: useful hints, no resolved advantage over the present</h2><p>The circuit target is the difference in the next 20-unit mean voltage between a changed-current copy (±0.1) and an unchanged copy. Present features contain both individual voltages, the old bias and the intervention bias; a simpler comparator uses total voltage. History adds total voltages 1, 2, 5 and 10 units earlier. Complete initial phase pairs stay in one split. The matched null shuffles whole histories within each prior-bias stratum and split, preserving observation/arm slots.</p><p>The fluid target is the integrated normalized energy difference from unchanged forcing over times 0–8. The physical protocol reuses the prior validated periodic 3D incompressible Navier–Stokes solver: grid 18³, ν=0.2, dt=0.02, two-unit burn and 12-unit response, unchanged/sudden/gradual forcing arms. Present features are eight <em>global flow summaries</em> (energy, enstrophy, helicity, forcing alignment, low-mode energy fraction, maximum speed, fourth speed moment and RMS acceleration) plus intervention indicators. History adds their changes over 0.5 and 1 time unit. These are not a complete velocity field or pointwise local measurements, and no new marker-history experiment was run here.</p><p>Both studies use the same RBF kernel-ridge family, training-only scaling, and hyperparameters selected only on validation runs. Models are frozen before test scoring. Old fluid test seeds are unused. Bootstrap intervals resample whole test runs (4,000 resamples); they are exploratory, with no multiplicity adjustment and one shuffled-history realization.</p>{pred_table}
<figure><img src="figures/held_out_prediction.png"><figcaption>Lower error is better. Units differ between the circuit and fluid panels. A paired difference interval containing zero does not resolve improvement.</figcaption></figure><p>Circuit history slightly worsens mean test error; its paired interval includes zero. Fluid history lowers mean error by about 24%, but the paired interval also includes zero. Its advantage over this shuffled-history control is {ci(fp['history_minus_shuffled_history'])}. That comparison alone does not establish benefit beyond instantaneous measurements. The complete phase-state physics comparator predicts the circuit response to a maximum error of {cp['full_state_physics_max_error']:.3g}; extra measured history can only help recover information omitted from partial observables in this deterministic model.</p></section>
<section><h2>Numerical checks and reproducibility</h2>{validate_table}<p><a href="results/verification.json">Saved-data audit</a> · <a href="tests.xml">Automated test record</a> · <a href="manifest.json">File checksums and provenance</a>. All current result tables are under <a href="results/">results/</a>. The linked source repository provides executable Python, frozen models, exact checkpoints, recordings, restoration checksums and reproduction instructions. The original protocol is preserved; the added staircase refinements are explained in <a href="results/ramps_refined.json">their own record</a>.</p></section>
<section><h2>What remains proposed</h2><p><strong>Physical calibration and laboratory validation:</strong> obtain the circuit diagram, I<sub>c1</sub>, I<sub>c2</sub>, shunt resistances, load impedance, capacitance, temperature, acquisition bandwidth, current waveform and time-resolved voltage. Convert to dimensionless units using the supplied calibration helper; estimate instrument noise separately from thermal noise. Use one dataset to fit parameters, then freeze them and predict withheld ramp speeds, sweep directions and repeated starts. Compare voltage traces, switching distributions and observation-window dependence, with uncertainties from independent repetitions. This has not been done because no hardware or measurement dataset was supplied.</p><p><strong>Additional simulations:</strong> validate an explicitly specified LRC topology and its stored-energy budget; test capacitive states away from synchrony; combine noise, mismatch, capacitance and ramps; expand independent test sets without reusing these test outcomes for model tuning; repeat shuffled nulls and check stronger present-state features and history horizons. Larger/noisier fluid ensembles, boundaries, actual marker geometry and other media require their own physical models.</p><p><strong>Transfer to fluids:</strong> the new fluid forecast replication asks the same statistical question, but it does not identify Josephson orbit invariants in water. To claim a shared mechanism, first define and measure the proposed fluid structures, map interventions to measurable physical quantities, compare against established advection/diffusion/stretching and dissipation explanations, and validate on independent flows and measurements. Similar pictures, marker separation, or lower prediction error alone would not establish that mechanism.</p><p>No laboratory, calibrated-temperature, LRC, general-capacitance, fog/ice or circuit-to-fluid causal result is claimed. No new law is inferred from the established models tested here.</p></section>'''
    css='''body{font:16px/1.65 system-ui,sans-serif;color:#172a36;background:#f3f6f7;margin:0}header,section{max-width:1080px;margin:24px auto;padding:32px;background:white;border-radius:12px}h1{font-size:40px;line-height:1.15}h2{font-size:25px;color:#174e67}.lead{font-size:21px}.eyebrow{letter-spacing:.08em;font-size:12px;color:#476d7b}a{color:#086783}table{border-collapse:collapse;width:100%;font-size:14px}td,th{text-align:left;border-bottom:1px solid #dce5e8;padding:9px;vertical-align:top}th{background:#edf4f6}figure{margin:24px 0}img{width:100%;height:auto}figcaption{font-size:14px;color:#4c6270}@media print{body{background:white}section,header{margin:0;break-inside:auto}figure{break-inside:avoid}}'''
    (ROOT/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Josephson follow-ups — Streams and Rocks</title><style>'+css+'</style>'+body+'</html>',encoding='utf-8')
    readme=f'''# Josephson follow-ups and independent fluid forecast replication

Completed 10 October 2026. Dimensionless numerical experiments; no laboratory calibration or established circuit-to-fluid mechanism.

The ideal coupled circuit retains some up/down differences after whole-cycle averaging. The held-out prediction tests do **not** resolve a history advantage over present measurements: circuit history-minus-present error is {ci(cp['history_minus_present'])}; fluid is {ci(fp['history_minus_present'])} (95% paired bootstrap intervals).

| Completed study | Main result |
|---|---|
| Longer / whole-cycle observations | 324 checkpoints; six direct runs to 12,800. At α=0.5, i=1.02, whole-cycle mean absolute up/down gap ≈0.07171. Uncoupled whole-cycle gap is numerical zero. |
| Staircases and continuous ramps | Matched elapsed durations; six increments. Finest RMS discrepancies ≈0.000103 and 0.00000966 for the two speeds. Original nonmonotone refinement retained. |
| Dimensionless noise and mismatch | 24 independent starting/noise draws per regime; coupled step checks pass. No physical temperature fitted. |
| Capacitance | Three β values, up/down sweeps and time-step/settling checks, restricted to the exactly synchronized branch. β=10 supports different locked/running responses. |
| Held-out forecasts | 16 independent circuit test starts and 12 new fluid test fields. Fluid history beats one shuffled-history control, but improvement over present measurements is unresolved. |

[Full illustrated report](report.html) (download and open locally; GitHub shows HTML source) · [Parameter/result tables](results/) · [Protocol](protocol.json) · [Numerical audit](results/verification.json) · [Test record](tests.xml).

[Streams and Rocks report]({STREAMS}) · [Executable code, models and recordings]({CODE}).

![Whole-cycle study](figures/whole_cycles.png)

![Held-out forecast comparison](figures/held_out_prediction.png)

Remaining experiments are explicitly separate: physical calibration and laboratory validation; an LRC load; general off-synchrony capacitive dynamics; combined perturbations; larger locked validation sets; marker-geometry and other-media studies; and any demonstrated physical transfer to fluids. The report explains the data and controls each needs.

## Reproduce from the linked source repository

Download the code folder linked above, including `recordings/`. Restore the checksum-verified arrays before running the audit:

```text
python restore_recordings.py
python -m pip install -r requirements.txt
python -m pytest test_models.py test_inference.py -q
python verify.py
python run_all.py
```

The completed run passed **13 automated tests and {audit_count} saved-data checks**. `run_all.py` reruns all studies and rebuilds the figures/report; it takes several minutes. Fluid trajectories are cached with source and data hashes; changed solver settings require removing or relocating only the affected cached files before recomputation. `report.py` rebuilds the presentation from saved results without resimulating. Reproduction can differ in final floating-point digits and runtime metadata.

Raw circuit arrays and frozen models are in `data/`; new fluid arrays and hash-checked metadata are in `fluid/data/fluid/`. The restoration script verifies every archive part and refuses to overwrite different local files. `provenance.json` records the unchanged solver sources, inherited circuit checkpoints and prior fluid training/validation inputs. The fluid test fields are newly generated.

Array conventions: cycle case axes flatten `(bias, history mode, alpha, initial draw)`; noise arrays are `(time, bias, condition, draw, junction)`; capacitance arrays are `(sweep plateau, sample time, beta, draw, [phase, velocity, input work, dissipation])`. Predictive targets are changed-minus-unchanged responses. All voltages and times are dimensionless.

`calibrate.py` converts supplied component values into these units; no physical values were supplied or fitted in this study. The laboratory procedure and remaining controls are described in the full report.
'''
    (ROOT/'README.md').write_text(readme,encoding='utf-8')
    from confirmation_report import add_confirmation
    add_confirmation()
    print('Created five scientific figures, report, structured summary and repository README.',flush=True)
if __name__=='__main__':main()
