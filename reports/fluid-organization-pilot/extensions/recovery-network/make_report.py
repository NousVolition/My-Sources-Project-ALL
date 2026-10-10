"""Scientific plots and transparent report generated only from saved results."""
import csv,html,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,P

def read(name):return json.loads((ROOT/'results'/name).read_text(encoding='utf-8'))
def ci(d,digits=4):return f"{d['mean']:.{digits}g} [{d['ci95'][0]:.{digits}g}, {d['ci95'][1]:.{digits}g}]"
def figsave(fig,name):
    fig.savefig(ROOT/'figures'/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(ROOT/'figures'/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def table(headers,rows):
    return '<table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'
def image(name,caption):return f'<figure><img src="figures/{name}.png" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>'

def main():
    (ROOT/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    r=read('recovery_summary.json');rp=read('recovery_prediction.json');n=read('network_forecast.json');nh=read('network_history_prediction.json');nv=read('network_validation.json');par=read('learned_network_parameters.json');rec=read('recurrence.json')
    colors=['#537786','#d46a32','#17918a'];names=['Unchanged forcing','Abrupt shutoff','Gradual shutoff']
    z=np.load(ROOT/'results'/'recovery_curves.npz');t=z['time'];curves=z['energy']
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for a,name in enumerate(names):
        q=curves[:,:,a];mean=q.mean(0);lo,hi=np.quantile(q,[.1,.9],axis=0)
        ax[0].plot(t,mean,label=name,color=colors[a]);ax[0].fill_between(t,lo,hi,color=colors[a],alpha=.12)
        if a:
            vals=[]
            with (ROOT/'results'/'recovery_runs.csv').open() as f:
                for row in csv.DictReader(f):
                    if row['arm']==P['fluid']['arms'][a] and row['sustained_arrival']:vals.append(float(row['sustained_arrival']))
            ax[1].scatter(np.full(len(vals),a)+np.linspace(-.12,.12,len(vals)),vals,color=colors[a],s=18,alpha=.7)
    ax[0].axhline(.1,ls=':',color='black',label='Settling threshold');ax[0].set(xlabel='Time after intervention',ylabel='Kinetic energy / initial energy',title='36 independent starting flows');ax[0].legend(fontsize=9)
    ax[1].set(xticks=[1,2],xticklabels=['Abrupt','Gradual'],ylabel='Sustained settling time',title='Below 10% energy for one time unit')
    figsave(fig,'recovery')
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for j,(target,dd) in enumerate(rp.items()):
        models=['current','history','shuffled_history'];means=[dd[k]['mean'] for k in models]
        errors=np.array([[dd[k]['mean']-dd[k]['ci95'][0] for k in models],[dd[k]['ci95'][1]-dd[k]['mean'] for k in models]])
        ax[j].bar(range(3),means,yerr=errors,color=colors,capsize=4)
        ax[j].set(xticks=range(3),xticklabels=['Present','+ History','+ Shuffled'],ylabel='Mean test-run RMSE',title=['Integrated recovery energy','Change relative to unchanged forcing'][j])
    fig.suptitle('12 independent test flows; 95% run-bootstrap intervals');figsave(fig,'recovery_prediction')
    fig,ax=plt.subplots(2,1,figsize=(11,6),layout='constrained',sharex=True)
    sample=np.load(ROOT/'data'/'network'/'cycle_c1_seed2024.npz')
    for i,color in enumerate(colors):
        ax[0].plot(sample['time'],np.exp(sample['base'][:,i]),color=color,label=f'State component {i+1}')
        ax[1].plot(sample['time'],sample['base'][:,i],color=color)
    ax[0].legend(ncol=3,fontsize=9);ax[0].set(ylabel='Component amplitude',title='Specified deterministic three-saddle cycle')
    ax[1].set(xlabel='Model time',ylabel='Log amplitude',title='Small components remain represented; no artificial floor')
    figsave(fig,'network_cycle')
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    ex=np.load(ROOT/'results'/'network_forecast_example.npz')
    for key,style in [('truth','-'),('learned','--'),('fixed_prior',':')]:ax[0].plot(ex['time'],ex[key][:,0],style,label=key.replace('_',' '))
    ax[0].set(xlabel='Forecast time',ylabel='Component 1',title='Independent disturbed start\nNo teaching signal');ax[0].legend()
    models=['fixed_prior','learned','shuffled_teacher','linear'];means=[n[k]['mean'] for k in models]
    errors=np.array([[n[k]['mean']-n[k]['ci95'][0] for k in models],[n[k]['ci95'][1]-n[k]['mean'] for k in models]])
    ax[1].bar(range(4),means,yerr=errors,capsize=3,color=['#9caaad','#17918a','#d46a32','#6971a3']);ax[1].set_yscale('log')
    ax[1].set(xticks=range(4),xticklabels=['Fixed prior','Learned','Shuffled','Linear'],ylabel='Mean test-run state RMSE',title='30-unit autonomous forecasts\n12 test runs; 95% intervals')
    figsave(fig,'network_learning')
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    d=nh['restricted_wait_to_next_dominant_state'];models=['label_only_markov','age_conditioned_semi_markov','current','history','shuffled_history','learned_dynamics']
    errors=np.array([[d[k]['mean']-d[k]['ci95'][0] for k in models],[d[k]['ci95'][1]-d[k]['mean'] for k in models]])
    ax[0].barh(range(6),[d[k]['mean'] for k in models],xerr=errors,capsize=3,color=['#a8b6bd','#7e929c','#537786','#17918a','#d46a32','#344f73'])
    ax[0].text(.15,5,f"{d['learned_dynamics']['mean']:.4f}",va='center',fontsize=9)
    ax[0].set(yticks=range(6),yticklabels=['Label only','Label + age','Full present, ridge','Present + history','Shuffled history','Learned dynamics'],xlabel='Mean test-run waiting-time RMSE',title='Waiting time capped at 30\nCensored observations retained')
    with (ROOT/'results'/'network_regimes.csv').open() as f:regimes=list(csv.DictReader(f))
    for c in P['network']['coupling_values']:
        rr=[q for q in regimes if q['variant']=='cycle' and float(q['coupling'])==c and q['arm']=='base' and int(q['seed']) in P['network']['test_seeds'][:8]]
        ax[1].scatter([c]*len(rr),[int(q['transitions_between_saddle_visits']) for q in rr],s=30,alpha=.65)
    ax[1].set(xlabel='Off-diagonal competition multiplier',ylabel='Transitions between axial-state visits / 180 units',title='Visits near axial states\nDistance < 0.2')
    figsave(fig,'history_and_coupling')
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for a,name in enumerate(names):
        ds=np.stack([np.load(ROOT/'data'/'recurrence'/f'seed{s}.npz')['distance_to_initial'][:,a] for s in P['fluid']['test_seeds']])
        ax[0].plot(t,ds.mean(0),color=colors[a],label=name)
    ax[0].axhline(.1,color='black',ls=':');ax[0].axhline(.2,color='gray',ls='--');ax[0].legend(fontsize=8)
    ax[0].set(xlabel='Time',ylabel='Relative full-field distance from initial state',title='12 test flows: departure and possible return')
    for label,col in [('time','#17918a'),('grid','#d46a32')]:
        rr=[q for q in r['validation']['refinement'] if q['comparison']==label]
        ax[1].scatter([q['seed'] for q in rr],[q['max_energy_difference_over_initial'] for q in rr],label=label,color=col)
    ax[1].axhline(.01,color='black',ls=':',label='Declared 1% gate');ax[1].set_yscale('log');ax[1].legend()
    ax[1].set(xlabel='Independent test seed',ylabel='Max energy-curve difference / initial energy',title='Time and spatial refinement assessed separately')
    figsave(fig,'recurrence_and_validation')
    recovery_rows=[]
    for arm in P['fluid']['arms']:
        q=r[arm];recovery_rows.append([arm,q['settled'],q['censored'],ci(q['settling_time_among_observed']) if q['settling_time_among_observed'] else 'Not observed',ci(q['integrated_energy_0_8'])])
    preds=[]
    for target,dd in rp.items():
        preds.append([target,ci(dd['current']),ci(dd['history']),ci(dd['shuffled_history']),ci(dd['history_minus_current'])])
    regress=[[k,ci(d[k])] for k in models]
    forecasts=[[k,ci(n[k]),ci(n['paired_disturbance_response'][k])] for k in ('oracle','fixed_prior','learned','shuffled_teacher','linear')]
    body=f'''<h1>Switching, recovery, and learned dynamics</h1>
<p class="lead">Completed numerical extension pilot · 10 October 2026</p>
<div class="status">Two experiments were run: a three-arm forcing intervention in 3D incompressible flow, and a separate three-saddle competition model with autonomous learned forecasts. These results do not identify a heteroclinic network in the fluid.</div>
<p><strong>Main findings:</strong> gradual shutoff delays observed settling by about {r['gradual_minus_abrupt_settling']['mean']:.2f} model-time units on average. Fluid history lowers integrated-energy prediction error against the present-only model, but its paired advantage over shuffled history has an interval spanning zero. The intervention-specific energy difference likewise has no statistically resolved history advantage. In the specified network, frozen coefficient learning forecasts independent trajectories accurately; full-state dynamics outperform the history regression.</p>
<p><a href="README.txt">Run instructions</a> · <a href="protocol.json">Declared protocol</a> · <a href="results/recovery_runs.csv">Fluid results table</a> · <a href="results/network_regimes.csv">Network regimes</a> · <a href="results/tests.xml">Automated tests</a></p>
<h2>1. What was completed</h2>
<p>The fluid ensemble has 36 independent starting fields, each copied into unchanged, abrupt-shutoff, and gradual-shutoff branches (108 primary trajectories). Eight additional three-branch calculations separately check halved time steps and finer grids on four seeds. A separate preliminary screen is excluded from ensemble statistics. The network has 36 primary paired starts, 40 paired coupling/topology controls and 12 paired numerical refinements. Each prediction study uses 18 training, 6 validation and 12 entirely separate test initializations. Multiple branches or times from one initialization never count as independent replicates.</p>
<h2>2. Fluid model and intervention</h2>
<p>We integrate ∂u/∂t = −P[(u·∇)u] + ν∇²u + A(t)f₀ on a periodic cube of side 6, with ∇·u=0, ν=0.2, 18³ grid points and RK4 step 0.02. Strict two-thirds Fourier dealiasing uses the existing pilot operators. f₀ is a fixed solenoidal Beltrami field with RMS 0.06. A=1 for the control, A=0 after abrupt shutoff, and A=max(1−t/2,0) for gradual shutoff. Both interventions reach zero forcing, but their cumulative input work differs and is recorded.</p>
<p>Random low-mode initial fields have RMS speed drawn from 0.4–0.8 and evolve for 2 time units before copying. This is a reproducible transient baseline, <strong>not a claim of statistical stationarity</strong>. The post-intervention observation lasts 12 units. All lengths, times and coefficients are model units; these are not new water/air material calibrations. This extension observes field organization; it introduces no forces between passive markers.</p>
<p>Settling means kinetic energy at or below 10% of its value at shutoff, continuously for one observed time unit. Samples are 0.1 units apart, so crossing times have that resolution. First arrivals and sustained arrivals are saved separately. Runs without a complete residence window are right-censored. The table reports mean and 95% run-bootstrap intervals; censored runs are counted and excluded only from the explicitly conditional observed-settling-time mean.</p>
{table(['Arm','Settled /36','Censored','Observed settling time [95% CI]','∫₀⁸ E/E₀ dt [95% CI]'],recovery_rows)}
<p>Gradual minus abrupt settling time, among jointly observed pairs: {ci(r['gradual_minus_abrupt_settling'])}. This difference includes the finite ramp duration and extra forcing work; it is not itself a memory measurement.</p>
{image('recovery','Mean energy curves; shaded ranges are the 10th–90th percentiles across starting fields, not confidence intervals. The adjacent points show each observed settling time.')}
<p>Overshoot is max(E/E₀−1,0), saved per trajectory. The abrupt unforced branch must lose total energy; any substantial energy rise would fail the physics check. Threshold sensitivity covers energy fractions 0.05, 0.1, 0.2 and residence times 0.5, 1, 2 in <a href="results/recovery_summary.json">the full summary</a>.</p>
<h2>3. Does pre-intervention organization predict recovery?</h2>
<p>Present observations are energy, enstrophy, helicity, forcing alignment, low-wave-number energy fraction, maximum speed, fourth speed moment and RHS norm, plus intervention identity. The history adds changes over the previous 0.5 and 1 units. Every model has the same quadratic present-state basis; history terms are additional linear inputs. Shuffling complete history blocks between runs at the same intervention supplies a comparator with identical input count. Regularization is chosen using validation runs only.</p>
<p>The first target is integrated normalized energy over 0–8 units. The second is its paired difference from unchanged forcing. These defined finite-horizon quantities retain all runs regardless of settling. They are related recovery measures, not interchangeable with settling-time prediction or future marker deformation.</p>
{table(['Target','Present RMSE','+History RMSE','Shuffled-history RMSE','Paired history − present RMSE'],preds)}
{image('recovery_prediction','Mean per-run prediction error on 12 held-out initializations; uncertainty resamples complete runs. Negative paired history-minus-present error indicates improvement.')}
<p>For integrated recovery energy, the paired history-minus-shuffled RMSE difference is {ci(rp['integrated_normalized_energy_0_8']['history_minus_shuffled'])}. For the intervention-specific target, history-minus-present is {ci(rp['integrated_difference_from_unchanged_0_8']['history_minus_current'])}. Both intervals span zero. The modest first-target improvement against present-only regression therefore remains tentative after the capacity-matched shuffled control; the second target does not establish the main disturbance-response hypothesis.</p>
<p>The current predictor observes only eight summaries of a large velocity field. A history benefit here would show utility under incomplete observation, not an additional constitutive memory law. Nonlinear predictor families, more complete present-state inputs, more independent flows and measurement noise remain necessary before claiming generality. The test outcomes are exploratory across the reported targets; no multiple-testing correction is claimed.</p>
<h2>4. Returns and numerical consistency</h2>
<p>A descriptive replay of all 12 test starts measures the relative spatial L2 distance from the complete velocity field at shutoff. After departure beyond radius 0.2, a return requires distance ≤0.1 for at least 0.2 units. There were <strong>{rec['total_returns']} observed returns</strong> across these 36 branches. Radii 0.05, 0.1 and 0.2 are also recorded. This measures return near that one reference field, not all past states, and does not demonstrate an attracting basin. Replay energy curves match the original recordings.</p>
{image('recurrence_and_validation','Departure from the initial full velocity field and independent time/grid controls. The recurrence diagnostic was added descriptively after the initial network analysis and did not select prediction models.')}
<p>Maximum normalized post-intervention energy-budget residual: {r['validation']['max_energy_budget_residual']:.3g}; maximum divergence: {r['validation']['max_divergence']:.3g}; maximum advective CFL: {r['validation']['max_cfl']:.3g}. Exact cloning error: {r['validation']['clone_error_max']:.3g}. The budget includes stage-time forcing work and viscous loss. Refinement compares (18³,0.02) with (18³,0.01), then (18³,0.01) with (26³,0.01). Passing means maximum energy-curve difference below 1% of initial energy and settling-time difference ≤0.11. These four-seed checks do not certify other Reynolds numbers or wall boundaries.</p>
<h2>5. Specified heteroclinic network</h2>
<p>The separate model is xᵢ′=xᵢ[1−ΣⱼCᵢⱼxⱼ], with Cᵢᵢ=1, Cᵢ,ᵢ₊₁=1.6c, Cᵢ,ᵢ₊₂=0.6c and cyclic indices. At c=1, the three axial equilibria have eigenvalues −1, −0.6 and +0.4. Invariant coordinate planes contain the connections 1→2→3→1, with cycle contraction product (0.6/0.4)³=3.375. An independent boundary-orbit test checks a representative connection. This is the standard three-component competition construction, not the supplied nine-saddle graph or the cited learning paper's exact template.</p>
<p>Initial positive components are uniform on [0.05,0.8], followed by a 20-unit preparation and a 180-unit observation. The disturbed copy receives +0.2 in the log amplitude of the next cyclic component at the intervention, equivalent to a multiplicative amplitude kick. Integration uses log coordinates and DOP853, rtol 10⁻⁹, atol 10⁻¹¹, maximum step 0.2. There is no added noise and no artificial positive floor. Refinements halve the maximum step and tighten tolerances tenfold.</p>
{image('network_cycle','An independent test trajectory approaches successive saddle states. Smaller components and increasing residence times remain visible in log coordinates.')}
<p>Off-diagonal multipliers c=0,0.8,1,1.2 test uncoupled growth, stable coexistence and two cyclic regimes. We also remove one inhibitory matrix entry and reverse one row's two off-diagonal weights, preserving its total strength. The physical competition matrix and the phase-space transition graph are distinct. The regime CSV reports transition counts and completed dominance intervals, with first/last censoring. In a stable coexistence regime, arbitrarily small oscillations or ties can change the largest-component label; this is not evidence of switching between saddle states.</p>
<h2>6. Frozen learning and autonomous prediction</h2>
<p>Using only baseline training trajectories, nonnegative least squares estimates the nine competition coefficients from sampled log derivatives. Parameters are frozen before independent test trajectories start. The forecast is initialized once from each test state and runs 30 units without teacher forcing or access to later observations. Both unchanged and kicked starts are evaluated. Comparators are the exact known equations (oracle), a declared wrong fixed prior (a=1.3,b=0.8), a derivative/state-pair shuffled identification null, and a learned unconstrained linear vector field with nine coefficients. Equal parameter count does not imply equal representational power; the linear model is a deliberately restrictive comparator.</p>
<p>Maximum absolute coefficient error: {par['max_abs_coefficient_error']:.4g}. This deliberately favorable identification experiment has noiseless complete state measurements and the correct functional model class. It validates a workflow; it does not establish learning from unknown fluid dynamics or robustness to observation noise.</p>
{table(['Predictor','State RMSE [95% CI]','Paired kick-response RMSE [95% CI]'],forecasts)}
{image('network_learning','Frozen learned dynamics predict new starts without a continuing teaching signal. A correct functional form and complete clean observations make this a favorable control.')}
<h2>7. History, the next state and its waiting time</h2>
<p>The next destination in this simple cycle is already determined by the current label. A transition table learned from training data achieves {nh['next_destination']['label_only_accuracy']:.1%} accuracy on {nh['next_destination']['n_observed_test_transitions']} observed test transitions. This is a topology control, not a history discovery.</p>
<p>The harder target is time until the next largest-component change, capped at 30 units. All {nh['n_test_observations']} test observations are retained, including {nh['censored_test_observations']} capped cases. This restricted waiting time is not the uncensored mean dwell time. Observation times are fixed in the protocol. Present measurements include all three amplitudes, their logs, instantaneous log rates and intervention identity. History includes preceding changes, elapsed dominance age and the previous completed dominance duration. A label-only constant-hazard model and a censoring-aware age-conditioned empirical survival model provide coarse baselines. The exact and learned ODE forecasts use only the current full state.</p>
{table(['Waiting-time predictor','Mean run RMSE [95% CI]'],regress)}
{image('history_and_coupling','History helps a restricted regression model, while a frozen full-state dynamical model supplies the stronger comparator. Transition counts under different coupling are descriptive, not automatic evidence for heteroclinic behavior.')}
<p>History-minus-present regression RMSE: {ci(d['history_minus_current'])}. The known full-state dynamics remain sufficient in this specified deterministic model. An improvement over a finite regression basis is a statement about representation and forecasting, not evidence that two identical complete states have different futures.</p>
<h2>8. What remains open</h2>
<ul><li>Stationary or turbulent fluid baselines, new forcing shapes and endpoints, walls, obstacles, and broader independent regimes.</li><li>Settling-time prediction with survival methods, beyond the completed finite-horizon energy prediction.</li><li>More complete current-state predictors, matched nonlinear alternatives, sensitivity of fitted predictions to grid refinement, measurement noise and larger test ensembles.</li><li>A branching or nine-saddle network; multiple possible next states; varying noise with a declared stochastic model and time-step controls.</li><li>The cited online teacher-coupled adaptation rule, unknown model classes, and out-of-regime learned forecasting.</li><li>Independent identification of fluid invariant states and connecting orbits. No such identification is claimed here.</li><li>Josephson sweeps, phase changes, and ice remain outside these two completed pilots.</li></ul>
<h2>Sources and provenance</h2>
<p>The forcing-intervention comparison is an original fluid pilot inspired by the distinction between state kicks and parameter changes in <a href="https://man-aravind.github.io/assets/pdf/heteroclinic_relaxation.pdf">Aravind &amp; Meyer-Ortmanns (2023)</a>; it does not reproduce that paper's neural interpretation or nested network. The independent-learning controls are motivated by <a href="https://www.frontiersin.org/journals/applied-mathematics-and-statistics/articles/10.3389/fams.2019.00063/full">Voit &amp; Meyer-Ortmanns (2019)</a>, with a different, explicitly stated offline method. The competition model follows the class studied by <a href="https://doi.org/10.1137/0129022">May &amp; Leonard (1975)</a>.</p>
<p>Fluid Fourier operators are unchanged from the original pilot, whose upstream numerics were pinned to My-Sources-Project-ALL commit 1b12558739736d472b0eb75845cf721371fcbb71. See <a href="manifest.json">file hashes</a>, <a href="results/verification.json">verification gates</a>, and <a href="requirements-lock.txt">recorded dependencies</a>. This extension leaves the original pilot results intact.</p>'''
    css='body{font:17px/1.65 system-ui,sans-serif;max-width:1080px;margin:40px auto;padding:0 24px;color:#20343e}h1,h2{color:#175770}h1{font-size:36px}.lead{font-size:21px}.status{background:#ecf5f3;border-left:5px solid #17918a;padding:16px}table{border-collapse:collapse;width:100%;font-size:14px;display:block;overflow:auto}td,th{padding:10px;border-bottom:1px solid #ccd6db;text-align:left}th{background:#eef3f5}img{width:100%;height:auto}figure{margin:30px 0}figcaption{font-size:14px;color:#4a606a}a{color:#126588}'
    (ROOT/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Recovery and network pilot</title><style>'+css+'</style>'+body+'</html>',encoding='utf-8')
    print('Six scientific figures and report generated',flush=True)

if __name__=='__main__':main()
