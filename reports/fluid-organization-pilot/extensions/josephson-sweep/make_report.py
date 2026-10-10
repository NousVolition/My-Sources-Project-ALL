"""Generate scientific figures and a report from saved accepted measurements."""
from pathlib import Path
import csv,json,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import from_state,voltage,synchronized_voltage
from analyze import canonical
ROOT=Path(__file__).resolve().parent
CODE='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-organization-pilot/extensions/josephson-sweep'
STREAMS='https://github.com/NousVolition/Nous-Volition/tree/main/studies/fluid-organization/extensions/josephson-sweep'


def main():
    c=json.loads((ROOT/'protocol.json').read_text());s=json.loads((ROOT/'results/summary.json').read_text())
    r=json.loads((ROOT/'results/reference_summary.json').read_text())
    if not r['refinement']['all_passed'] or not r['power_gate_passed']:raise RuntimeError('Refined references must pass')
    figures=ROOT/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titlesize':12,'figure.dpi':120,'savefig.dpi':160})
    colors=['#156c9c','#ca671d','#198777'];names=['Fresh starts','Upward sweep','Downward sweep']
    b=np.array(c['biases']);a=np.array(c['alphas']);n=c['n_independent_initializations']
    v={}
    for setting in c['settings']:
        with np.load(ROOT/f'data/{setting["name"]}.npz') as z:v[setting['name']]=canonical(z['voltage']).sum(-1)
    with np.load(ROOT/'data/reference_windows.npz') as z:p=from_state(z['states'],z['sign'])
    windows={int(w):canonical(voltage(p[0],p[round(w/c['record_stride'])],w)).sum(-1) for w in c['window_checks']}
    def save(fig,name):
        fig.savefig(figures/f'{name}.png',bbox_inches='tight');fig.savefig(figures/f'{name}.svg',bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(13,4.2),sharex=True)
    for k,ax in enumerate(axs):
        for m in range(3):
            ax.plot(b,v['base'][:,m,k].mean(-1),color=colors[m],label=names[m])
            ax.plot(b,v['both'][:,m,k].mean(-1),color=colors[m],ls=':',lw=2)
        dense=np.linspace(0,2,501)
        ax.plot(dense,synchronized_voltage(dense,a[k]),color='#777777',ls='--',lw=1,label='Synchronized branch')
        ax.axvline(1,color='#aaa',lw=.8);ax.set(title=f'Load ratio α = {a[k]:g}',xlabel='Normalized bias i')
        ax.grid(alpha=.15)
    axs[0].set_ylabel('Mean total voltage / (Ic r)');axs[2].legend(fontsize=8)
    fig.suptitle('Dimensionless Josephson voltage sweeps');fig.tight_layout();save(fig,'voltage_sweeps')
    fig,axs=plt.subplots(1,3,figsize=(13,4.2))
    mask=(b>=.98)&(b<=1.02)
    for k,ax in enumerate(axs):
        for m in range(3):
            ax.plot(b[mask],windows[200][mask,m,k].mean(-1),'-o',ms=3,color=colors[m],label=names[m])
            ax.plot(b[mask],windows[800][mask,m,k].mean(-1),'--',color=colors[m])
            low,high=np.quantile(windows[800][mask,m,k],[.05,.95],axis=-1)
            ax.fill_between(b[mask],low,high,color=colors[m],alpha=.10)
        ax.axvline(1,color='#aaa',lw=.8);ax.set(title=f'α = {a[k]:g}',xlabel='Normalized bias i');ax.grid(alpha=.15)
    axs[0].set_ylabel('Mean total voltage / (Ic r)');axs[2].legend(fontsize=8)
    fig.suptitle('Critical window: W=200 solid, W=800 dashed; 5–95% spread shaded')
    fig.tight_layout();save(fig,'critical_window')
    fig,axs=plt.subplots(1,3,figsize=(13,4.2),sharey=True)
    step=np.maximum(abs(v['half_step']-v['base']),abs(v['both']-v['double_settle']))
    with np.load(ROOT/'data/reference_S100.npz') as z:r100=canonical(z['voltage']).sum(-1)
    with np.load(ROOT/'data/reference_S200.npz') as z:r200=canonical(z['voltage']).sum(-1)
    settle=abs(r200-r100);wd=abs(windows[800]-windows[200])
    for k,ax in enumerate(axs):
        for values,label,color in [(step,'Half time step','#385c9d'),(settle,'Double settling','#bd551c'),(wd,'Longer window','#247e63')]:
            ax.semilogy(b,np.maximum(values[:,:,k].max(axis=(1,2)),1e-14),'-o',ms=3,label=label,color=color)
        ax.axhline(c['gates']['max_absolute_voltage_step_difference'],color='#888',ls=':',label='Step acceptance limit')
        ax.set(title=f'α = {a[k]:g}',xlabel='Normalized bias i');ax.grid(alpha=.15)
    axs[0].set_ylabel('Maximum absolute total-voltage change');axs[2].legend(fontsize=8)
    fig.suptitle('Separate numerical step, settling and observation effects');fig.tight_layout();save(fig,'separate_checks')
    rng=np.random.default_rng(c['bootstrap_seed']);rs=rng.integers(0,n,(c['bootstrap_resamples'],n))
    fig,axs=plt.subplots(1,3,figsize=(12,4.1),sharey=True)
    for k,ax in enumerate(axs):
        means=[];lows=[];highs=[]
        for w in windows:
            gap=abs(windows[w][mask,1,k]-windows[w][mask,2,k]).mean(0)
            lo,hi=np.quantile(gap[rs].mean(-1),[.025,.975]);means.append(gap.mean());lows.append(lo);highs.append(hi)
        ax.errorbar(list(windows),means,yerr=[np.array(means)-lows,np.array(highs)-means],fmt='o-',capsize=4,color=colors[k])
        ax.set(title=f'α = {a[k]:g}',xlabel='Voltage averaging window W',xticks=[200,400,800]);ax.grid(alpha=.15)
    axs[0].set_ylabel('Mean absolute up/down voltage gap')
    fig.suptitle('History-associated gaps: critical-current average, 95% paired bootstrap intervals')
    fig.tight_layout();save(fig,'history_windows')
    with np.load(ROOT/'data/precision_diagnostic.npz') as z:
        t=z['time'];ref=z['reference_0'];log=z['log_state_0'];direct=[z[f'direct_0_{dt}'] for dt in [.05,.025,.0125]]
    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    for d,dt in zip(direct,[.05,.025,.0125]):axs[0].plot(t,d.sum(-1),label=f'Direct phases, dt={dt}')
    axs[0].plot(t,ref.sum(-1),'k--',lw=2,label='Log-separation reference')
    axs[0].set(xlabel='Dimensionless time τ',ylabel='Unwrapped total phase',title='A failed direct-phase validation case')
    axs[0].legend(fontsize=8);axs[1].plot(t,log[:,1],color='#156c9c')
    axs[1].axhline(28,color='#c26336',ls=':',label='Severe cancellation: separation about 10⁻¹²')
    axs[1].set(xlabel='Dimensionless time τ',ylabel='z = log |tan(q/2)|',title='Separation modulo 2π can be lost')
    axs[1].legend(fontsize=8);fig.tight_layout();save(fig,'precision_control')
    def ci(d):return f"{d['mean']:.6f} [{d['ci95'][0]:.6f}, {d['ci95'][1]:.6f}]"
    rows=''.join(f"<tr><td>{rr['alpha']:g}</td><td>{old['max_step_change']:.3g}</td><td>{rr['max_settling_change']:.3g}</td><td>{rr['max_window_change']:.3g}</td><td>{ci(rr['critical_mean_abs_up_down_gap_W200'])}</td><td>{ci(rr['critical_mean_abs_up_down_gap_W800'])}</td></tr>" for old,rr in zip(s['regimes'],r['regimes']))
    refinement_rows=''.join(f"<tr><td>{name}</td><td>{html.escape(str(val))}</td></tr>" for name,val in r['refinement']['comparisons'].items())
    gates=''.join(f"<tr><td>{html.escape(k.replace('_',' '))}</td><td>{g.get('value','—')}</td><td>{g.get('limit','—')}</td><td>{'PASS' if g['passed'] else 'FAIL'}</td></tr>" for k,g in s['numerical_gates'].items())
    settings=''.join(f"<tr><td>{x['name']}</td><td>{x['dt']}</td><td>{x['settle']}</td><td>200</td></tr>" for x in c['settings'])
    caption='Dimensionless Josephson sweep. Solid curves use the initial integration settings; dotted curves use half the step and twice the settling time. Critical-window voltage can remain sensitive to history and observation time.'
    report=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Josephson voltage sweep — completed numerical pilot</title>
<style>body{{font:17px/1.65 system-ui,sans-serif;max-width:1140px;margin:44px auto;padding:0 24px;color:#243a45}}h1,h2{{color:#164d65;line-height:1.2}}h1{{font-size:40px}}h2{{margin-top:42px}}.status,.equation{{padding:16px 22px;background:#eef5f7;border-left:4px solid #26768c}}.caution{{background:#fff6e6;border-color:#b8791d}}a{{color:#146787}}figure{{margin:30px 0}}img{{max-width:100%}}figcaption{{font-size:14px;color:#485a63}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #d7e2e6}}th{{background:#eef5f7}}.scroll{{overflow-x:auto}}code{{font-size:14px}}footer{{border-top:1px solid #ccd;margin-top:40px;font-size:14px}}</style></head><body>
<p>STREAMS AND ROCKS · COMPLETED NUMERICAL EXTENSION · 10 OCTOBER 2026</p>
<h1>Josephson voltage sweeps: history, observation time and numerical precision</h1>
<div class="status"><strong>Completed:</strong> {s['counts']['primary_segments']:,} fixed-bias segments across all four original step/settling settings, 3,672 adaptive-reference segments, 612 targeted refinement segments, longer-window replays, 24 independent-start solver checks and nine exact-cycle controls. Eighteen automated tests pass. An original sweep refinement gate failed; additional checks below resolve the affected history without erasing that failure. This is a dimensionless model study, not device measurements or a fluid experiment.</div>
<p><a href="{CODE}">Runnable code, tests and recorded arrays</a> · <a href="results/reference_regimes.csv">Resolved regime table</a> · <a href="results/reference_summary.json">Reference summary and uncertainty</a> · <a href="protocol.json">Original protocol and coordinate amendment</a> · <a href="refinement_protocol.json">Additional refinement protocol</a></p>
<h2>What the completed experiment establishes</h2>
<p>The original step comparison changes total voltage by as much as <strong>{s['numerical_gates']['max_absolute_voltage_step_difference']['value']:.3g}</strong>, exceeding the declared 2×10⁻⁵ limit in part of the α=0.5 downward sweep. Further halving to Δτ=0.0125 and 0.00625, plus independently integrated full-history references, resolves that case. The final finest-step/reference discrepancy there is <strong>{r['refinement']['comparisons']['dt00625']['max_total_voltage_difference']:.3g}</strong>. The results below distinguish that numerical problem from the larger settling and observation effects.</p>
<p>Upward, downward and freshly initialized runs give different finite-observation voltages after refinement. The reference table shows which differences decrease with a longer window and which remain in this pilot. Full-history references cover all original currents and starting draws; further fixed-step refinement is targeted to the history that failed its original gate.</p>
<p>These results test the earlier sweep caption quantitatively. They do not establish an additional physical memory law or asymptotic hysteresis. In this model, the two phases and the prescribed current determine the future; a scalar averaged voltage does not fully specify those phases.</p>
<figure><img src="figures/voltage_sweeps.png" alt="Three load-ratio voltage sweeps comparing fresh, upward and downward protocols"><figcaption>{caption} These are the original Δτ=0.05/S=100 and Δτ=0.025/S=200 comparisons, whose step error was not uniformly converged. Curves are means over 12 starting phase pairs. Gray dashed curves are the analytic synchronized branch, not a prediction for every arbitrary phase arrangement.</figcaption></figure>
<h2>Model, assumptions and physical budget</h2>
<p>The supplied circuit is interpreted as two identical overdamped resistively shunted Josephson junctions in series, each with shunt r, with a shared load R across the pair. The conventional RSJ approximation neglects junction capacitance; see <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7314811/">Couëdo et al. (2020)</a>. More general loaded arrays are treated by <a href="https://arxiv.org/abs/chao-dyn/9408004">Wiesenfeld and Swift (1995)</a>; their LRC model is background, not an experimental validation of this idealized pair. No capacitance, inductance, thermal evolution, noise or junction mismatch is included here.</p>
<div class="equation">α = r/R, &nbsp; i = I<sub>b</sub>/I<sub>c</sub>, &nbsp; τ = (2e I<sub>c</sub> r / ℏ)t<br>
(1+α)φ′<sub>1</sub> + αφ′<sub>2</sub> = i − sin φ<sub>1</sub><br>
αφ′<sub>1</sub> + (1+α)φ′<sub>2</sub> = i − sin φ<sub>2</sub><br>
v<sub>j</sub> = V<sub>j</sub>/(I<sub>c</sub>r) = φ′<sub>j</sub>; &nbsp; v<sub>total</sub> = v<sub>1</sub> + v<sub>2</sub>.</div>
<p>The shared branch carries normalized current αv<sub>total</sub>. α=0 removes that load (R→∞ with r fixed). Voltage averages use unwrapped phase increments over the observation window, retaining both individual junction voltages and their sum.</p>
<div class="equation">U = −cos φ<sub>1</sub> − cos φ<sub>2</sub><br>
i v<sub>total</sub> = dU/dτ + (φ′<sub>1</sub>)² + (φ′<sub>2</sub>)² + α(v<sub>total</sub>)².</div>
<p>The resistors dissipate nonnegative energy. Input work, potential change and accumulated dissipation are checked for each segment. This is a driven dissipative system, so its potential energy alone is not conserved. On the invariant synchronized branch, v̄<sub>total</sub>=2√(i²−1)/(1+2α) for i&gt;1, and the locked branch has zero voltage below threshold. The single-cycle period 2π(1+2α)/√(i²−1) grows near i=1, making a fixed observation interval particularly sensitive there.</p>
<h2>Reproducible comparison</h2>
<p>Three load ratios (0, 0.5, 2), 17 currents from 0 to 2, and 12 independent phase pairs drawn uniformly from [−π,π)² with seed 20261010. Current samples are concentrated between 0.98 and 1.02. Each pair is reused across conditions as a paired control; it is not counted as a new independent replicate for every current.</p>
<table><tr><th>Setting</th><th>Step Δτ</th><th>Settling time S</th><th>Averaging time W</th></tr>{settings}</table>
<p>Up and down are separately initialized monotone sweeps. Each next current starts from the endpoint after S+200 at the preceding current. Fresh runs reset to the original phases at every current. Settling is a prescribed waiting interval, not a claim that every trajectory has reached its asymptotic regime. The four settings use identical current grids and starting draws; changing S also changes the accumulated history of a continuation sweep.</p>
<p>For the finest step with S=200, each stored post-settling state is replayed for W=800. Nested W=200, 400 and 800 averages and four disjoint 200-unit blocks are measured. Longer observations do not feed back into the primary sweep continuation. All quantities are dimensionless; no physical component values or voltages in volts are inferred.</p>
<figure><img src="figures/separate_checks.png" alt="Separate time-step, settling and averaging-window changes"><figcaption>Maximum absolute total-voltage changes over the 12 starts and three protocols at each current. The step curve uses the larger change from the two original fixed-settling comparisons. Settling and observation curves use the full-history adaptive references: S=100 versus 200 with W fixed at 200; W=200 versus 800 from the same settled state. Floors at 10⁻¹⁴ are for plotting only.</figcaption></figure>
<h2>Quantitative regimes and uncertainty</h2>
<div class="scroll"><table><tr><th>α</th><th>Original max step change</th><th>Reference max settling change</th><th>Reference max W change</th><th>Reference critical up/down gap, W=200 [95% CI]</th><th>Reference critical up/down gap, W=800 [95% CI]</th></tr>{rows}</table></div>
<p>The first three columns are maxima across all currents, modes and starting draws; they are sensitivities, not error bars. The gap is the mean absolute up-minus-down voltage difference over the seven sampled currents from 0.98 through 1.02, first computed separately for each starting draw. Intervals use 4,000 paired bootstrap resamples of the 12 independent draws (seed 10012026), preserving currents and modes together. These intervals describe this initial-phase distribution, not device-to-device uncertainty. A gap near a precision floor should not be interpreted as physically nonzero.</p>
<figure><img src="figures/critical_window.png" alt="Critical-current voltages at short and long observation windows"><figcaption>Adaptive-reference mean voltage with W=200 (solid) and W=800 (dashed), both after S=200. Every long-window reference was also replayed with RK4 at Δτ=0.025 from its identical settled state. Shading is the 5–95% range over initial draws at W=800, not a confidence interval.</figcaption></figure>
<figure><img src="figures/history_windows.png" alt="Up-down gap versus averaging window with bootstrap intervals"><figcaption>Paired uncertainty over starting phase draws. Nested windows are dependent measurements. A remaining W=800 gap is a completed finite-time observation; it does not establish whether the gap persists at arbitrarily long times.</figcaption></figure>
<h2>A numerical artifact identified and corrected</h2>
<div class="status caution">The initial direct-phase implementation passed its power budget but failed the voltage refinement checks. Its independent-solver disagreement reached 0.03015 for an individual junction. These superseded numbers are kept in <a href="{CODE}/audit">the audit record</a> and are excluded from the accepted tables.</div>
<p>Phases can approach equality modulo 2π so closely that rounding and subtraction errors severely distort their small separation, and can eventually erase it. That separation can later grow again. Merely halving the time step does not reliably recover information already lost to rounding. An adaptive solver using the same ill-conditioned phase coordinates can also be affected. The illustrated reference reaches z≈29.9, corresponding to a phase separation modulo 2π of roughly 4×10⁻¹³.</p>
<div class="equation">p=(φ<sub>1</sub>+φ<sub>2</sub>)/2, &nbsp; q=(φ<sub>1</sub>−φ<sub>2</sub>)/2, &nbsp; z=log|tan(q/2)|<br>
p′=(i+sin p tanh z)/(1+2α), &nbsp; z′=−cos p.<br>
The sign of q is retained separately. Exact q=0 is handled as z=−∞.</div>
<p>These are algebraically equivalent equations on the invariant initial phase strip, not a new circuit model. They preserve small separations internally even when reconstructed phases round to the same value modulo 2π. Voltage is still measured from unwrapped phases. Tests compare both formulations before severe conditioning develops, validate the transformed velocities against Kirchhoff's equations, and retain the synchronized special case. Acceptance thresholds were fixed before production and were not relaxed after failure.</p>
<figure><img src="figures/precision_control.png" alt="Direct-phase failure and log-separation coordinate diagnostic"><figcaption>An independent validation start at α=2 and i=1.005. Direct-phase integrations at three step sizes diverge from the accurately resolved trajectory. Large positive z corresponds to q approaching π, which is synchronization modulo 2π. The danger scale is illustrative; actual rounding depends on phase magnitude and arithmetic.</figcaption></figure>
<h2>Validation and reproducibility</h2>
<p>The original fixed-step matrix is retained below, including its failed step gate. Passing fresh-start checks and an energy budget did not ensure an entire continuation history was accurate. The original fixed-step means are therefore labeled as the initial comparison; resolved gap estimates use separate full-history references.</p>
<table><tr><th>Check</th><th>Measured</th><th>Acceptance limit</th><th>Status</th></tr>{gates}</table>
<p><strong>Additional completed checks:</strong> DOP853 replays both entire S=100 and S=200 sweeps from the original initial conditions (relative tolerance 10⁻¹¹, absolute tolerance 10⁻¹², maximum step 0.5). The α=0.5, S=200 downward histories are rerun with fixed steps 0.0125 and 0.00625 and with a tighter DOP853 reference (10⁻¹²/10⁻¹³, maximum step 0.25). All four follow-up gates pass. The finest-step pair and finest-step/reference limits remain 2×10⁻⁵; the reference-tolerance limit is 2×10⁻⁶. All long-window reference trajectories are cross-checked with fixed RK4 from exactly the same settled states. These checks do not certify arbitrary longer histories or other circuits.</p>
<table><tr><th>Follow-up comparison</th><th>Measured discrepancy</th></tr>{refinement_rows}</table>
<p>Independent validation uses 24 new initial draws from seed 20261011, covering α=0,0.5,2 and i=0.99,1,1.005,1.2. Fixed RK4 results are compared with DOP853 in the stable coordinates (relative tolerance 10⁻¹¹, absolute tolerance 10⁻¹², maximum step 0.5). Nine exact synchronized-cycle checks cover the three loads and i=1.005,1.2,2. Eighteen automated tests cover circuit balance, passive dissipation, stationary and periodic solutions, coordinate equivalence, small separation, unwrapped voltage, refinement and replay.</p>
<p><a href="results/primary_voltages.csv">Original four-setting voltages</a> · <a href="results/reference_voltages.csv">Resolved reference voltages</a> · <a href="results/reference_means.csv">Reference means and intervals</a> · <a href="results/reference_windows.csv">Reference nested windows</a> · <a href="results/reference_disjoint_blocks.csv">Reference disjoint blocks</a> · <a href="results/independent_validation.json">Independent validation</a> · <a href="results/verification.json">Original numerical gates</a> · <a href="results/refinement.json">Additional refinement gates</a> · <a href="tests.xml">Automated test record</a>. Source, protocol and recorded-array hashes are in the code package; arrays retain initial, settled and endpoint states and sampled histories.</p>
<h2>Physical interpretation and relation to the fluid question</h2>
<p>The current-dependent Josephson phase dynamics, the resistive load, different phase arrangements, critical slowing and finite observation windows explain the responses available in this model. At α=2 and i=1.005 the synchronized period is about 314 time units: W=200 does not even cover one cycle and W=800 covers only about 2.6. Dissipation and coupling are linked by the circuit equations; increasing α does not isolate a universal interaction-strength variable independently of transport and damping.</p>
<p>The model also preserves extraordinarily small phase differences in an exactly identical, noiseless circuit. Whether such sensitivity survives component mismatch and environmental fluctuations is a physical question that this idealized calculation cannot answer.</p>
<p>Previous observations can help infer phase information missing from an instantaneous scalar voltage. That is consistent with the broader organization/history question, but this sweep does not fit or validate a history-based predictor against a complete-state predictor. It supplies no evidence that fluid markers interact like junction phases, nor that the fluid has this circuit's phase states. The initial numerical failure also shows why apparent history effects require a precision check as well as an energy check.</p>
<h2>Proposed experiments — not completed here</h2>
<ul><li>Extend settling and observation through additional orders of magnitude, use full-cycle averages and analyze phase-orbit structure to distinguish persistent branch selection from long transients.</li><li>Vary sweep increments and dwell protocols, and include a genuinely continuous current ramp; this pilot uses 17 discrete current plateaus.</li><li>Add controlled junction mismatch and calibrated thermal noise, then test whether phase-sensitive differences survive those physical perturbations.</li><li>Add capacitance or an LRC load only with the corresponding equations, new energy checks and convergence controls; inertial hysteresis is outside this overdamped pilot.</li><li>Train and freeze voltage-history predictors on separate draws, then test against instantaneous voltage and the complete phase state. No Josephson learning result is claimed here.</li><li>Validate any laboratory comparison using measured component values and uncertainty. Transfer of these mechanisms to the fluid experiments requires separate evidence.</li></ul>
<footer>Completed numerical results belong to this extension. Earlier fluid and recovery/network results are preserved. <a href="{STREAMS}">Streams and Rocks</a> · <a href="{CODE}">Source repository</a>.</footer></body></html>'''
    (ROOT/'report.html').write_text(report,encoding='utf-8')
    readme=f'''JOSEPHSON VOLTAGE SWEEP — COMPLETED NUMERICAL PILOT

Read report.html for results, figures, equations, uncertainty and limitations.
Streams and Rocks: {STREAMS}
Code and recordings: {CODE}

Python 3.11+; install requirements.txt. requirements-lock.txt records this run.
  python -m pip install -r requirements.txt
  python -m pytest test_model.py -q
  python run_all.py

Individual stages: run_sweep.py, validate.py, diagnose_precision.py,
analyze.py, refine.py, analyze_reference.py, make_report.py, verify.py.
The original fixed-step gate remains a reported FAIL; all follow-up refinement
gates must pass for the final report. run_sweep.py caches source/protocol-matched
recordings. If changing the model/protocol, copy the package to a fresh directory
without data/ or use --out NEW_DIRECTORY for numerical exploration. Analysis
scripts read the package data/ directory. Do not reuse old data after a change.

If recordings are delivered as archive parts, run restore_recordings.py first.
All subsequent analysis can use recorded data without rerunning the sweeps.
  python analyze.py
  python analyze_reference.py
  python make_report.py
  python verify.py

The accepted equations use a mean/log-separation coordinate transform of the
original conditional circuit_rhs in the fluid pilot. No new physical memory
term is added. audit/ preserves the initial failed direct-phase checks and
instructions to reproduce them separately. Their raw traces are superseded.

Array layout in base/half_step/double_settle/both.npz:
  voltage: [bias index, mode, alpha, initialization, junction]
  states: [bias index, recorded time, mode, alpha, initialization, (p,z)]
  sign: [mode, alpha, initialization]; initial/settled/endpoint store (p,z).
  bias_order[:,0:3]: actual fresh/up/down current sequences.
Descending mode is reversed only when forming increasing-current tables.
windows.npz: states [recorded time, bias index, mode, alpha, initialization, (p,z)].
Use model.from_state(states,sign) to reconstruct unwrapped phases. Preserve z
when restarting: reconstructing it from rounded output phases loses precision.
Every primary and window trajectory is sampled every 5 dimensionless time units;
integration uses the separately recorded .05/.025 steps. Voltages use exact
integration endpoints, not differencing the coarsely sampled output.

All voltages are normalized by Ic*r. All durations use tau=(2e Ic r/hbar)t.
No dimensional device parameters or laboratory measurements are supplied.
'''
    (ROOT/'README.txt').write_text(readme,encoding='utf-8')
    print('Created five PNG/SVG figures, report.html and README.txt')


if __name__=='__main__':main()
