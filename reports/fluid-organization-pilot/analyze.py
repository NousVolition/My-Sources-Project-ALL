"""Rebuild tables, static scientific plots and the self-contained HTML report."""
import base64
import csv
import html
import json
import platform
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import ROOT,PilotFlow,wrap,sha,save_json,L
from prediction import neighbors,edges

RESULTS=ROOT/'results';FIG=RESULTS/'figures';DATA=ROOT/'data'
COLORS=['#116b83','#df8435','#8c4e9f','#23825c','#b54d50','#5b6da8']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.prop_cycle':plt.cycler(color=COLORS),'savefig.facecolor':'white'})


def read(path):return json.loads(Path(path).read_text())


def csv_write(path,rows):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def figsave(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close(fig)


def summary(path):
    m=read(path);a=m['diagnostics'][0];b=m['diagnostics'][-1];c=m['config'];e0=a['perturbation_energy'];z=np.load(path.with_suffix('.npz'))
    ix=np.argmin(abs(z['time']-c['kick_time']));x=z['positions_base'][ix];nn=neighbors(x);e=edges(x,nn);ef=edges(z['positions_perturbed'][-1],nn)
    newnn=neighbors(z['positions_perturbed'][-1]);ret=np.mean([len(set(u)&set(v))/8 for u,v in zip(nn,newnn)])
    return {'name':m['name'],'seed':m['seed'],'n':c['n'],'dt':c['dt'],'nu':c['nu'],'speed':c['speed'],'kick_rms':c['amplitude'],
            'forcing':c['forcing'],'nonlinear':c['nonlinear'],'initial':c['initial'],'Re_box':c['speed']*L/c['nu'],
            'energy_gain':b['perturbation_energy']/e0 if e0 else 0.,
            'baseline_energy_retained':b['baseline_energy']/a['baseline_energy'],
            'disturbed_energy_retained':b['disturbed_energy']/a['disturbed_energy'],
            'production_over_Edelta0':b['integrated_production']/e0 if e0 else 0.,
            'dissipation_over_Edelta0':b['integrated_perturbation_dissipation']/e0 if e0 else 0.,
            'source_radius_final':b['radius_from_source'],'source_radius_rate':(b['radius_from_source']-a['radius_from_source'])/(c['end']-c['kick_time']),
            'participation_fraction_final':b['participation_fraction'],
            'marker_pair_rms':b['marker_difference_rms'],'neighbor_length_ratio':float(np.median(np.linalg.norm(ef,axis=-1)/np.linalg.norm(e,axis=-1))),
            'neighbor_retention':ret,'vorticity_pair_correlation':b['vorticity_pair_correlation'],
            'vorticity_initial_correlation':b['vorticity_pattern_correlation_initial'],
            'vorticity_rms_change':b['vorticity_difference_rms'],'max_high_band_fraction':max(r['high_band_energy_fraction'] for r in m['diagnostics']),
            'max_budget_residual':max(abs(r[k]) for r in m['diagnostics'] for k in ['baseline_energy_residual','disturbed_energy_residual','perturbation_budget_residual']),
            'max_divergence':max(r['divergence_max'] for r in m['diagnostics']),'max_cfl':m['max_cfl'],'seconds':m['seconds']}


def aggregate(rows,p):
    out=[];rng=np.random.default_rng(117)
    for label in p['regimes']:
        rr=[r for r in rows if r['name'].startswith(label+'_')]
        if not rr:continue
        row={'regime':label,'runs':len(rr),'nu':rr[0]['nu'],'speed':rr[0]['speed'],'kick_rms':rr[0]['kick_rms'],'Re_box':rr[0]['Re_box']}
        for k in ['energy_gain','source_radius_rate','participation_fraction_final','marker_pair_rms','neighbor_length_ratio','neighbor_retention','vorticity_pair_correlation','production_over_Edelta0','dissipation_over_Edelta0']:
            vals=np.array([r[k] for r in rr]);bs=vals[rng.integers(0,len(vals),(4000,len(vals)))].mean(1)
            row[k+'_mean']=float(vals.mean());row[k+'_ci_low'],row[k+'_ci_high']=map(float,np.quantile(bs,[.025,.975]))
        out.append(row)
    return out


def compare_numerics(p):
    rows=[]
    for seed in p['test_seeds']:
        base=summary(DATA/'ensemble'/f'seed{seed}.json');z=np.load(DATA/'ensemble'/f'seed{seed}.npz')
        for kind in ['grid','dt']:
            r=summary(DATA/'refinement'/f'{kind}_{seed}.json');zz=np.load(DATA/'refinement'/f'{kind}_{seed}.npz')
            rows.append({'regime':'reference','seed':seed,'comparison':kind,'gain_relative_change':abs(r['energy_gain']/base['energy_gain']-1),
                         'marker_rms_change':float(np.sqrt(np.mean(np.sum(wrap(zz['positions_perturbed'][-1]-z['positions_perturbed'][-1])**2,axis=1)))),
                         'fine_budget_residual':r['max_budget_residual']})
    for name in ['fast_water','water','air']:
        base=summary(DATA/'regimes'/f'{name}_7000.json');z=np.load(DATA/'regimes'/f'{name}_7000.npz')
        for kind in ['grid','dt']:
            r=summary(DATA/'refinement'/f'{name}_{kind}_7000.json');zz=np.load(DATA/'refinement'/f'{name}_{kind}_7000.npz')
            rows.append({'regime':name,'seed':7000,'comparison':kind,'gain_relative_change':abs(r['energy_gain']/base['energy_gain']-1),
                         'marker_rms_change':float(np.sqrt(np.mean(np.sum(wrap(zz['positions_perturbed'][-1]-z['positions_perturbed'][-1])**2,axis=1)))),
                         'fine_budget_residual':r['max_budget_residual']})
    return rows


def plots(p,reg,numerics,pred,sens):
    fig,axs=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for label in ['water','air','fast_water','stokes_null']:
        ms=[read(DATA/'regimes'/f'{label}_{s}.json') for s in p['regime_seeds']]
        t=np.array([r['t'] for r in ms[0]['diagnostics']])-p['base']['kick_time']
        for ax,key,norm in [(axs[0,0],'perturbation_energy',True),(axs[0,1],'radius_from_source',False),(axs[1,0],'participation_fraction',False),(axs[1,1],'marker_difference_rms',False)]:
            vals=np.array([[r[key] for r in m['diagnostics']] for m in ms])
            if norm:vals=vals/vals[:,0,None]
            line=ax.plot(t,vals.mean(0),label=label.replace('_',' '))[0]
            ax.fill_between(t,vals.min(0),vals.max(0),color=line.get_color(),alpha=.12)
            ax.set_xlabel('Time since fluid kick (model units)');ax.grid(alpha=.2)
    for ax,title,y in zip(axs.flat,['Disturbance energy','Distance from the source','Spatial extent of disturbance','Corresponding marker displacement'],['Eδ(t) / Eδ(0)','Energy-weighted RMS radius','Participation volume / box volume','RMS displacement']):ax.set(title=title,ylabel=y)
    axs[0,0].axhline(1,color='.4',ls=':',lw=1);axs[0,0].legend(fontsize=9)
    fig.suptitle('Paired 3D simulations | lines: four-seed mean; bands: observed seed range',fontsize=13)
    figsave(fig,'disturbance_response')
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    m=read(DATA/'regimes/fast_water_7000.json');d=m['diagnostics'];t=np.array([r['t'] for r in d])-.4;e=d[0]['perturbation_energy']
    axs[0].plot(t,[r['perturbation_energy']/e for r in d],label='Disturbance energy')
    axs[0].plot(t,[1+r['integrated_production']/e-r['integrated_perturbation_dissipation']/e for r in d],'--',label='1 + production − dissipation')
    axs[0].set(title='Where growth comes from',ylabel='Normalized disturbance energy');axs[0].legend(fontsize=8)
    axs[1].plot(t,[r['baseline_energy']/d[0]['baseline_energy'] for r in d],label='Baseline')
    axs[1].plot(t,[r['disturbed_energy']/d[0]['disturbed_energy'] for r in d],label='Disturbed',ls='--')
    axs[1].set(title='Both total energies decay',ylabel='Total kinetic energy / initial');axs[1].legend(fontsize=8)
    axs[2].plot(t,[r['enstrophy_stretching_rate'] for r in d],label='Vortex stretching contribution')
    axs[2].set(title='3D vortex stretching is retained',ylabel='Enstrophy production rate')
    for ax in axs:ax.set_xlabel('Time since kick');ax.grid(alpha=.2)
    fig.suptitle('Mechanism example | fast-water parameters, seed 7000');figsave(fig,'energy_budget')
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    labels=['current','history','shuffled_history','frozen_rate','training_mean']
    for ax,target,title in zip(axs,['future_deformation','disturbance_response'],['Overall future deformation','Additional deformation caused by kick']):
        means=[pred[target][k]['mean_run_rmse'] for k in labels]
        ax.bar(np.arange(len(labels)),means,color=[COLORS[0],COLORS[1],'#aab7be','#bac9bc','#d4d8dd'])
        for i,k in enumerate(labels):ax.scatter(np.full(8,i),pred[target][k]['per_seed_rmse'],color='#263542',s=15,alpha=.75)
        ax.set_xticks(range(len(labels)),['Current','+ history','Shuffled\nhistory','Frozen\npair rate','Train\nmean'],fontsize=9)
        ax.set(title=title,ylabel='Held-out RMSE (lower is better)');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Prediction from incomplete observations | eight independent test runs');figsave(fig,'history_prediction')
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for label in ['grid','dt']:
        rr=[r for r in numerics if r['comparison']==label and r['regime']=='reference']
        axs[0].semilogy(range(1,len(rr)+1),[max(r['gain_relative_change'],1e-15) for r in rr],'o-',label=label)
    axs[0].set(title='Perturbation-energy sensitivity',xlabel='Held-out run',ylabel='Relative change in final gain');axs[0].legend()
    for j,target in enumerate(['future_deformation','disturbance_response']):
        delta=pred[target]['paired_history_vs_current'];ci=delta['ci95']
        vals=[delta]+[sens[l][target]['paired'] for l in ['grid','dt']]
        for k,v in enumerate(vals):
            y=j*4+k;mean=v['mean_B_minus_A'];lo,hi=v['ci95']
            axs[1].errorbar(mean,y,xerr=[[mean-lo],[hi-mean]],fmt='o',color=COLORS[j],capsize=3)
    axs[1].axvline(0,color='.4',ls=':');axs[1].set_yticks([0,1,2,4,5,6],['Overall: base','Overall: finer grid','Overall: half step','Response: base','Response: finer grid','Response: half step'])
    axs[1].set(title='History − current prediction error',xlabel='Mean RMSE difference; 95% run-bootstrap interval')
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'numerical_controls')
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    pdata=[read(DATA/'particles'/f'air_particles_{s}.json') for s in p['regime_seeds']]
    for label in ['fog_gravity','fog_zero_g','smoke_gravity','smoke_zero_g']:
        rr=[[r for r in m['rows'] if r['species']==label] for m in pdata];t=[r['t'] for r in rr[0]]
        for ax,key in zip(axs,['paired_displacement_rms','distance_from_passive_rms']):
            vals=np.array([[r[key] for r in rows] for rows in rr]);line=ax.plot(t,vals.mean(0),label=label.replace('_',' '),ls='--' if 'zero' in label else '-')[0]
            ax.fill_between(t,vals.min(0),vals.max(0),color=line.get_color(),alpha=.1)
    axs[0].set(title='Particle response to gas disturbance',ylabel='RMS paired displacement')
    axs[1].set(title='Inertia / gravity versus passive tracers',ylabel='RMS distance from passive paths')
    for ax in axs:ax.set_xlabel('Time since gas kick');ax.grid(alpha=.2)
    axs[1].legend(fontsize=8);fig.suptitle('Dilute air-carried particles | four gas seeds; bands show seed range');figsave(fig,'particles')
    z=np.load(RESULTS/'phase_locking.npz');ph=np.linspace(-np.pi,2*np.pi,500)
    fig=plt.figure(figsize=(11,7.2),layout='constrained');gs=fig.add_gridspec(2,3)
    for j,mu in enumerate([0.,.5,1.5]):
        ax=fig.add_subplot(gs[0,j]);ax.plot(ph,mu-np.sin(ph));ax.axhline(0,color='.5',lw=.8)
        if mu<1:
            ax.scatter([np.arcsin(mu)],[0],color=COLORS[0],s=50,label='Stable phase')
            ax.scatter([np.pi-np.arcsin(mu)],[0],facecolors='white',edgecolors=COLORS[0],s=50,label='Unstable phase')
        ax.set(title=f'μ = {mu:g}: '+('phase lock' if mu<1 else 'phase drift'),xlabel='Phase difference φ (rad)',ylabel='dφ/dτ = μ − sin φ');ax.grid(alpha=.2)
    ax=fig.add_subplot(gs[1,:2])
    for mu in [0.,.5,.9,1.,1.2,1.5]:
        arr=z[f'pair{mu}'];sep=abs(np.angle(np.exp(1j*(arr[:,1]-arr[:,0]))))
        ax.semilogy(z['pair_time'],np.maximum(sep,1e-14),label=f'μ={mu:g}')
    ax.set(title='Identical state copies, then one phase kick of 0.05 rad',xlabel='Time after phase kick (τ)',ylabel='Circular phase separation');ax.legend(ncol=3,fontsize=8);ax.grid(alpha=.2)
    ax=fig.add_subplot(gs[1,2]);A=np.linspace(0,2,301);ax.plot(A,np.sqrt(np.maximum(1-A*A,0)),label='Analytic')
    pp=read(RESULTS/'phase_locking.json');ax.scatter([r['A'] for r in pp['coupling_controls']],[r['observed_window_drift'] for r in pp['coupling_controls']],color=COLORS[1],label='Numerical')
    ax.axvline(1,color='.4',ls=':');ax.set(title='Remove / restore coupling',xlabel='Coupling A; fixed detuning Δ=1',ylabel='Mean phase drift');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Phase-locking comparison from the supplied excerpt | a separate oscillator model');figsave(fig,'phase_locking')


def table(headers,rows):
    return '<table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'


def img(name,caption):
    data=base64.b64encode((FIG/(name+'.png')).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>'


def report(p,rows,reg,num,pred,sens,particle_checks):
    total=len(rows);maxbudget=max(r['max_budget_residual'] for r in rows);maxdiv=max(r['max_divergence'] for r in rows)
    gain=pred['future_deformation']['relative_rmse_reduction_percent'];response=pred['disturbance_response']
    text=['''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Fluid organization under stress — completed pilot</title><style>
    body{font:16px/1.6 system-ui,Segoe UI,sans-serif;color:#21323c;background:#f2f5f6;margin:0}main{max-width:1080px;margin:32px auto;background:white;padding:40px 48px}h1{font-size:38px;line-height:1.15;color:#104a60}h2{margin-top:2.4em;color:#116b83}h3{margin-top:1.6em}p{max-width:940px}figure{margin:30px 0}img{max-width:100%;height:auto}figcaption{font-size:13px;color:#52636c}table{border-collapse:collapse;width:100%;font-size:13px;margin:20px 0}td,th{padding:9px 10px;border-bottom:1px solid #dce4e8;text-align:left}th{background:#edf4f6}code,pre{background:#eef3f5;padding:2px 5px}pre{padding:16px;white-space:pre-wrap}.note{border-left:4px solid #da8e43;padding:10px 20px;background:#fff7ee}.lead{font-size:20px}.eq{padding:15px 20px;background:#edf4f6;font-family:Cambria,serif;font-size:19px}a{color:#116b83}li{margin-bottom:8px}@media(max-width:700px){main{margin:0;padding:22px}h1{font-size:30px}table{display:block;overflow-x:auto}}</style><main>
    <p>COMPLETED NUMERICAL PILOT · 9 OCTOBER 2026</p><h1>Fluid organization under stress</h1>''']
    text.append(f'<p class="lead">History improves prediction of overall marker deformation in this smooth-flow pilot. It does not improve prediction of the extra deformation caused by the disturbance.</p><p>On eight independent held-out runs, adding past neighboring-marker geometry reduces the mean run error for overall future deformation by <b>{gain:.1f}%</b>. For the paired disturbance response, the error changes from {response["current"]["mean_run_rmse"]:.6f} to {response["history"]["mean_run_rmse"]:.6f}. The latter does not support a useful history benefit. These are different targets, and their conclusions must remain separate.</p>')
    text.append(f'<p><b>{total} completed fluid configurations</b> include screening, 26 independent prediction seeds (12 train / 6 validation / 8 test), ten regimes with four additional seeds each, and numerical refinements. Configurations sharing a seed are dependent controls. Seven additional air-particle tracking configurations and a separate phase-locking experiment are complete. No results in this report are planned or fabricated values.</p>')
    text.append('<div class="note"><b>Scope:</b> this is a short, smooth, three-dimensional incompressible-flow pilot with periodic boundaries. It includes vortex stretching. It does not establish fully developed turbulence, molecular interactions, phase change, wall effects, ice fracture, a new fluid law, or a Navier–Stokes regularity theorem.</div>')
    text.append('<h2>1. What was reused and what changed</h2><p>The source project is <a href="https://github.com/NousVolition/My-Sources-Project-ALL/tree/1b12558739736d472b0eb75845cf721371fcbb71">NousVolition/My-Sources-Project-ALL, commit 1b125587…</a>. Its matched-stretch Fourier solver is preserved in <code>vendor/numerics.py</code>. This pilot reuses its projection, transforms, curl and Parseval inner product, and adapts the low-mode random initialization from the history study. New code supplies paired fluid disturbances, RK4 integration, strict two-thirds dealiasing, periodic cubic marker interpolation, energy budgets and the tests below.</p><p>The previous <a href="https://github.com/NousVolition/My-Sources-Project-ALL/blob/1b12558739736d472b0eb75845cf721371fcbb71/reports/matched-stretch/history-prediction/README.md">history study</a> reported a 1.025% gain for a different target (delayed interpolated-flow FTLE), with unresolved grid dependence. Its sharp-start marker-stress archives also failed some sensitivity checks. Those values are background, not rerun results or independent validation here. The new target is finite-neighborhood deformation, so the 71% and 1% gains are not comparable effect sizes.</p>')
    text.append('<h2>2. Equations and paired intervention</h2><div class="eq">∂u/∂t + (u·∇)u = −∇p + ν∇²u + f, &nbsp; ∇·u=0.<br>dX/dt = u(X,t).<br>Eδ = ½〈|uₚ−uᵦ|²〉; &nbsp; dEδ/dt = −〈δuᵢ δuⱼ ∂ⱼuᵦᵢ〉 − ν〈|∇δu|²〉.</div><p>Angle brackets denote box-volume averages; energies are per unit mass. A side-6 periodic cube starts from a smooth, exactly periodic low-mode field of prescribed RMS speed. It evolves for 0.4 model time units. The fluid and all marker positions are then copied exactly. Only the second fluid receives a solenoidal curl impulse with global RMS speed 0.05 (other amplitudes: 0.025 and 0.1). Both copies then evolve to time 1.0 under identical forcing and viscosity. The impulse is localized, smooth and band limited to the same Fourier modes on every grid; it has small nonzero spatial tails. Markers exert no force.</p><p>The first term in the difference-energy budget is transfer from the background velocity gradient; the second is viscous loss. Identical forcing cancels directly from this budget, but changes the background flow. The impulse also has a cross term with the original flow: Eδ(0) is not generally equal to the change in total kinetic energy. Production is measured from the projected RHS and checked against the explicit gradient contraction. Periodic pressure redistributes momentum without net energy work.</p><p>The construction and spectral method follow the upstream solver and the <a href="https://jipolanco.github.io/PencilFFTs.jl/dev/generated/navier_stokes/">PencilFFTs Navier–Stokes implementation</a>. The pilot uses classical RK4 with Δt=0.01, 26³ points, and 256 passive markers; controls use Δt=0.005 and 38³ points. The low-mode start and impulse represent the same physical Fourier polynomial on each grid.</p>')
    text.append('<h2>3. Water, air, and parameter regimes</h2><p>For an illustrative dimensional interpretation, take one model length unit as 1 mm and one velocity unit as 0.1 m/s. One time unit is 0.01 s, and the box is 6 mm across. Constant kinematic viscosities of 10⁻⁶ m²/s (water) and 1.5×10⁻⁵ m²/s (air) map to model ν=0.01 and 0.15. They are representative assigned inputs, not fitted thermodynamic measurements. The air speed is safely in the low-Mach regime. Matched dimensionless conditions would give the same incompressible velocity solution for either material; density would rescale energy per volume. This comparison isolates momentum diffusivity, not a universal ranking of the phases.</p>')
    tr=[]
    for r in reg:
        gaintext='identical copies' if r['regime']=='no_kick' else f'{r["energy_gain_mean"]:.3f} [{r["energy_gain_ci_low"]:.3f}, {r["energy_gain_ci_high"]:.3f}]'
        tr.append([r['regime'],r['nu'],r['speed'],r['kick_rms'],f'{r["Re_box"]:.0f}',gaintext,f'{r["source_radius_rate_mean"]:.3f}',f'{r["marker_pair_rms_mean"]:.4f}'])
    text.append(table(['Regime','ν','RMS speed','Kick RMS','Re (box)','Final Eδ/Eδ₀; 95% CI','RMS radius rate','Paired marker RMS'],tr))
    text.append('<p>Intervals resample the four independent seeds as whole runs; n=4 makes them preliminary. Parameter variants use the same seeds. “Fast water” doubles the background speed while keeping the absolute impulse fixed. “Forced” adds a fixed common solenoidal body force of amplitude 0.25 during the entire run. “Organized” uses an ABC Beltrami field. “Stokes null” removes nonlinear advection: it is a mathematical mechanism ablation, not a physically accurate creeping-flow approximation at the listed Reynolds number. A zero-kick energy gain has an undefined denominator and is displayed as identical copies; its numeric sentinel in CSV is zero.</p>')
    text.append(img('disturbance_response','Normalized energy, spatial spreading and marker motion answer different questions. Range bands show variability across four seeds, not confidence bands.'))
    text.append('<p>The source-centered RMS radius includes transport of the disturbance as well as broadening. The participation fraction, (Σe)²/(N³Σe²), measures its spatial extent without choosing a threshold. Their growth is not a finite signal-front speed: viscous diffusion and incompressible pressure do not provide a compact causal front. Periodic wrapping eventually limits both measures. Changes in vorticity, initial-pattern correlation, pair correlation, neighborhood retention and median edge-length ratio are saved in the full CSV. Vorticity correlation is an Eulerian pattern proxy, not an objective material-vortex or Lagrangian-coherent-structure classification.</p>')
    text.append(img('energy_budget','Energy-transfer example: perturbation growth coexists with decay of both total fluid energies. Vortex stretching is measured independently of particle separation.'))
    text.append('<h2>4. Does past organization improve prediction?</h2><p>At the kick, select each marker’s eight nearest neighbors using periodic distance and hold those labels fixed. The target is their mean log distance ratio between the kick and final time, divided by the 0.6 horizon. This is a finite-neighborhood statistic, not an infinitesimal FTLE. The second target subtracts the same statistic in the undisturbed copy. It isolates additional deformation under the specified intervention.</p><p>The current-only model receives 44 measurements: baseline velocity and full gradient, the applied velocity/gradient increment, post-kick strain eigenvalues, current neighborhood geometry and instantaneous neighbor-separation rates. A fixed quadratic expansion gives 1,034 features. The history model retains that expansion and adds 57 past-geometry features from 0.1, 0.2 and 0.3 time units before the kick. Features contain no seed, marker ID, future samples or future-derived neighbors. Histories are selected using the present neighbor labels only.</p><p>Training-only standardization and ridge fitting use 12 runs; six separate validation runs select regularization. Neither operation uses test outcomes. Eight untouched test runs supply the reported scores. The primary comparison and numerical controls use all eight runs, with 4,000 paired whole-run bootstrap samples. Markers in a flow are not treated as independent experimental replicates. Model coefficients are not interpreted as causal effects.</p>')
    hrows=[]
    for target,title in [('future_deformation','Overall deformation'),('disturbance_response','Extra disturbance response')]:
        r=pred[target];d=r['paired_history_vs_current']
        hrows.append([title,f'{r["current"]["mean_run_rmse"]:.6f}',f'{r["history"]["mean_run_rmse"]:.6f}',f'{d["mean_B_minus_A"]:.6f}',f'[{d["ci95"][0]:.6f}, {d["ci95"][1]:.6f}]'])
    text.append(table(['Held-out target','Current RMSE','+ history RMSE','History − current','95% paired interval'],hrows))
    text.append(img('history_prediction','Each dot is an entire held-out run. Shuffling histories within a run removes the large overall-deformation benefit.'))
    text.append('<p>History helps forecast overall deformation on this chosen distribution, including against the quadratic present-only and frozen-pair-rate controls. It fails to improve the specific disturbance-response target. Shuffled training labels provide a negative control. Reversing lag-feature columns and refitting gives the same answer, as ridge regression should under a column permutation; this is an invariance check, not evidence that chronological order is irrelevant. The prediction result is exploratory and does not test all model classes, marker densities or time horizons.</p><p>Past geometry can encode omitted spatial gradients, temporal trends and accumulated deformation. A full present fluid state already determines future evolution while a unique solution exists. Predictive value beyond this deliberately incomplete local observation therefore does not establish additional physical memory, a new force, or passive-marker interactions.</p>')
    text.append('<h2>5. Numerical validation and uncertainty</h2>')
    text.append('<div class="note"><b>Stress limit found:</b> the fast-water run at seed 7000 changes marker paths by 0.0348 RMS from 26³ to 38³, exceeding the declared 0.005 screen. Its final perturbation-energy gain changes by 1.73%, within the 2% energy screen. Thus the short-time energy-growth example is less sensitive than the detailed marker paths; fast-water trajectory organization is not validated on the original 26³ grid. The eight primary prediction runs have grid marker RMS changes no larger than 0.00178. The targeted follow-up refinement below investigates this failure.</div>')
    text.append(f'<p>All analytic and software checks pass (see <code>results/tests.xml</code>). Across the {total} fluid configurations, the largest normalized energy-budget residual is {maxbudget:.3g} and the largest sampled divergence magnitude is {maxdiv:.3g}. Tests include exact Beltrami decay, exact linear Stokes mode decay, incompressible projection, nonlinear energy and enstrophy identities, the perturbation-production identity, identical copies, consistent marker relabeling, periodic transport, interpolation refinement, drag relaxation, disjoint seed groups and phase-locking controls.</p>')
    text.append(table(['Refinement','Overall target RMS change','Response target RMS change','History-current overall RMSE','History-current response RMSE'],[[k,f'{sens[k]["target_rms_change"][0]:.6g}',f'{sens[k]["target_rms_change"][1]:.6g}',f'{sens[k]["future_deformation"]["paired"]["mean_B_minus_A"]:.6g}',f'{sens[k]["disturbance_response"]["paired"]["mean_B_minus_A"]:.6g}'] for k in ['grid','dt']]))
    text.append(img('numerical_controls','Refinement changes the measurements but preserves the distinct conclusions for overall deformation and disturbance-specific response.'))
    text.append(table(['Regime / control','Relative final-gain change','Marker RMS change'],[[r['regime']+' / '+r['comparison'],f'{r["gain_relative_change"]:.3g}',f'{r["marker_rms_change"]:.3g}'] for r in num if r['regime']!='reference']))
    text.append('<p>Numerical sensitivity is reported separately from run-to-run uncertainty. The finer-grid and half-step prediction checks reuse the trained/validated models and evaluate refined measurements for the same eight test seeds; they are dependent robustness checks, not extra independent validation. Cubic marker interpolation is more accurate than the earlier trilinear scheme but is not exactly divergence-free between nodes. Two grids and two time steps establish limited sensitivity evidence, not a rigorous continuum error bound. Only one seed per water/air/stress regime has extra refinement; broader high-stress conclusions remain limited.</p>')
    if (RESULTS/'stress_extension.json').exists():
        stress=read(RESULTS/'stress_extension.json')
        text.append('<h3>Targeted continuation of the stress case</h3><p>After the marker-path failure, the same fast-water seed was repeated on 50³ and 62³ grids, keeping Δt=0.005 fixed. These are dependent sensitivity runs. The original 26³ failure remains part of the result.</p>')
        text.append(table(['Grid comparison','Final marker RMS change','Relative energy-gain change','Finer-grid gain','Marker screen'],[[f'{r["from_n"]}³ → {r["to_n"]}³',f'{r["marker_rms_change"]:.6g}',f'{r["gain_relative_change"]:.6g}',f'{r["fine_gain"]:.6f}','PASS' if r['marker_screen_passed'] else 'FAIL'] for r in stress['comparisons']]))
        last=stress['comparisons'][-1]
        text.append('<p>'+('The finest comparison passes the declared marker screen for this seed. This supports the refined short-time trajectory result, not the original coarse-grid paths or every high-stress realization.' if last['marker_screen_passed'] else 'The finest comparison still fails the marker screen. The stress trajectory result remains unresolved, and further grid refinement is required before interpreting detailed organization.')+'</p>')
    text.append('<h2>6. Fog and smoke as particles in air</h2><div class="eq">dX/dt = V; &nbsp; dV/dt = [u(X,t)−V]/τₚ + g.</div><p>Four independently seeded air flows are replayed with two dilute heavy spherical particle classes, each with gravity present and removed. Fog uses τₚ=0.02 model time, liquid density 1,000 kg/m³ and diameter about 8.05 μm. Smoke uses τₚ=0.002, assigned solid density 1,800 kg/m³ and diameter about 1.90 μm. Particle velocities are copied continuously across the gas impulse; they do not jump to the perturbed gas velocity. Both classes have inertia. The carrier remains unchanged because loading is taken negligible.</p><p>This is the heavy-particle, small-Reynolds-number drag limit of the framework associated with <a href="https://www.dam.brown.edu/people/Maxey_Publications.html">Maxey &amp; Riley (1983)</a>. It omits added mass, the hydrodynamic history force, collisions, Brownian motion, particle feedback and temperature evolution. “History” in the marker-prediction experiment is not the Basset history force. The idealized smoke spheres do not represent combustion chemistry or real aggregate shapes. Fog remains liquid droplets suspended in gas; it is not another thermodynamic phase. Fixed droplets assume approximately saturated, isothermal air; evaporation and condensation were not simulated.</p>')
    text.append(img('particles','Gravity removal isolates a modeled environmental influence. Particle clustering would not imply carrier-fluid compression.'))
    text.append(table(['Particle sensitivity control','Final position RMS change','Largest particle Reynolds number'],[[r['comparison'],f'{r["position_rms_change"]:.3g}',f'{r["max_particle_Re"]:.4f}'] for r in particle_checks]))
    text.append('<p>Particle RK4 substeps satisfy Δt/τₚ≤0.25; the gas velocity is linearly interpolated between cubic spatial interpolants. A halved particle step, halved fluid step and finer fluid grid are checked on seed 7000. These are controlled dilute-particle examples, not validated atmospheric fog or smoke forecasts.</p>')
    text.append('<h2>7. Added comparison: phase locking and loss of synchrony</h2><p>The supplied Figure 4.5.1 excerpt motivates a separate phase-oscillator control. With φ = stimulus phase − oscillator phase, use the Adler equation below. The excerpt’s curve shapes correspond to the normalized form; its original source and edition are not visible in the attachment.</p><div class="eq">dφ/dt = Δ − A sinφ; &nbsp; Δ = Ω−ω.<br>For A&gt;0: μ=Δ/A, τ=At, &nbsp; dφ/dτ=μ−sinφ.</div><p>For |μ|&lt;1 there is a stable phase φ*=arcsin μ modulo 2π and an unstable phase π−arcsin μ. Equal frequencies (μ=0) allow in-phase locking; nonzero detuning locks with a constant offset. At |μ|=1 the fixed points merge in a saddle-node threshold. Robust locking uses the strict inequality |Δ|&lt;A; equality is a nonhyperbolic boundary with critical slowing. For |μ|&gt;1, the phase continually slips, with mean drift sign(μ)√(μ²−1) in normalized time. The instantaneous slip rate is nonuniform.</p>')
    text.append(img('phase_locking','Recomputed phase portraits, paired phase-kick recovery, and coupling removal/restoration. This is not a fitted fluid model.'))
    phase=read(RESULTS/'phase_locking.json')
    text.append(table(['μ','Stable phase (rad)','Recovery rate','Analytic mean slip','Numerical complete-cycle slip'],[[r['mu'],'threshold/drift' if r['stable_phase_rad'] is None else f'{r["stable_phase_rad"]:.4f}',f'{r["linear_recovery_rate"]:.4f}',f'{r["predicted_mean_slip_rate"]:.4f}',f'{r["complete_cycle_slip_rate"]:.4f}'] for r in phase['records']]))
    text.append('<p>Eight starting phases per parameter value were integrated to normalized time 400, with tighter-tolerance verification. Identical copies at time 20 receive a 0.05-radian kick in one copy and are followed for 20 more units. Coupling controls hold Δ=1 fixed while A=0, 0.5, 1 and 2 is removed/restored. This reproduces established synchronization behavior, as also described in <a href="https://arxiv.org/abs/2105.02564">Injection locking and synchronization in Josephson photonics devices</a>; it is not a new prediction.</p><p>The useful addition is a precise example in which coupling can maintain organization and a threshold can destroy it. It does not identify “coupling A” with viscosity or a marker-neighbor bond. To apply phase locking to the fluid, one must first demonstrate a reproducible oscillatory flow mode, define its phase, measure forcing detuning and strength, and test an entrainment interval against an uncoupled control. None of those fluid phase claims is made by this pilot. With full φ, Δ and A observed, this first-order oscillator also has no extra state memory that its past could add.</p><details><summary>Supplied reference excerpt</summary><img src="phase_locking_reference.png" alt="User-supplied Figure 4.5.1 excerpt showing phase locking and phase drift"><p>Supplied by the user; retained unchanged. The report’s own plots are generated from completed numerical controls.</p></details>')
    text.append('<h2>8. What the competing hypotheses establish here</h2>')
    text.append(table(['Hypothesis','Pilot result','Interpretation'],[
        ['H1: weakening interactions causes loss of organization','Not established as a general fluid mechanism','No microscopic coupling-strength parameter is present. Viscosity transfers momentum while dissipating kinetic energy; changing it cannot be labeled simply weakening or strengthening interaction. Phase locking supplies a separate example where reducing oscillator coupling loses synchrony.'],
        ['H2: coupling spreads and reorganizes a disturbance','Compatible with measured NS transfer and spreading','Advection, pressure coupling, background strain and viscous transport explain the paired results. Growth occurs when integrated production exceeds dissipation; total unforced kinetic energy still decreases. This is established continuum physics.'],
        ['H3: history adds predictive information','Supported for overall finite-neighborhood deformation; not for the disturbance-specific increment','The benefit concerns incomplete present observations on a short smooth-flow distribution. It survives the checked refinements but does not establish intrinsic fluid memory or causal marker arrangements.']]))
    text.append('<h2>9. Remaining unanswered questions and proposed extensions</h2><ol><li><b>Disturbance-specific prediction:</b> can history help at other forecast horizons, impulse locations, frequencies or strengths? Use new train/validation/test seed families and keep overall deformation separate from the intervention effect.</li><li><b>Stronger present-state comparators:</b> add present neighborhood velocity/gradient fields, Hessians, pressure information and matched-capacity nonlinear models. Determine whether history is only a proxy for missing current spatial information.</li><li><b>Stress and physical boundaries:</b> extend duration and Reynolds number only after multi-grid convergence, dissipation-scale resolution and interpolation checks. Add no-slip walls, inflow/outflow, moving boundaries and periodic-domain-size controls with an appropriate solver. These boundary comparisons are not completed.</li><li><b>Coherent structures and phase locking:</b> track objective material vortices or deformation barriers; identify a genuine oscillator before measuring fluid entrainment, phase slips or a locking interval. Pattern correlation alone is insufficient.</li><li><b>Fog and smoke:</b> include saturation/temperature transport, evaporation/condensation and latent heat for fog; particle size distributions, aggregate drag, Brownian motion, settling boundaries and finite loading where relevant. Water vapor requires a gas thermodynamic model if phase change or compressibility matters.</li><li><b>Ice — separate extension, not completed:</b> use solid displacement, stress and temperature, for example ρü=∇·σ with σ=C:ε and ρcṪ=∇·(k∇T)+Q, with appropriate anisotropy, boundary constraints and fracture/plasticity when needed. Couple a Stefan/enthalpy phase-change model only if melting is studied. Validate elastic-wave and heat-diffusion controls before paired localized thermal or mechanical impulses. Fluid Navier–Stokes alone is not an ice model.</li><li><b>Generality and uncertainty:</b> repeat on independent flow families and more seeds, extend marker densities and neighborhood scales, include measurement noise and test alternate definitions of organization. Current confidence intervals describe only the declared seed distribution.</li></ol>')
    text.append('<h2>10. Reproduce and inspect</h2><pre>python -m pip install -r requirements.txt\npython -m pytest test_pilot.py -q\npython run_experiments.py --stage all\npython particles.py\npython prediction.py\npython phase_locking.py\npython analyze.py\npython verify.py</pre><p>Run these commands from the <code>fluid_pilot</code> folder. Screening gates must pass before the larger suite. Existing fluid results are reused only when source hashes, configuration and data hashes match. A changed solver requires a fresh destination (for example <code>--data new_data</code>). Plots and tables rebuild from the saved results. The delivered package contains code, protocol, test results, raw paired marker trajectories, initial/final spectral fluid states, all diagnostics, prediction outputs, CSV tables and both PNG/SVG figures.</p>')
    text.append('<p><a href="results/parameter_regimes.csv">Regime table</a> · <a href="results/all_runs.csv">All fluid configurations</a> · <a href="results/numerical_controls.csv">Numerical controls</a> · <a href="results/prediction.json">Prediction metrics</a> · <a href="results/prediction_sensitivity.json">Prediction refinement</a> · <a href="results/phase_locking.json">Phase-locking results</a> · <a href="protocol.json">Protocol</a> · <a href="results/verification.json">Verification</a></p></main></html>')
    from dynamics import supplement
    from pendulum import supplement as pendulum_supplement
    from relaxation import supplement as relaxation_supplement
    from weak_nonlinear import supplement as weak_supplement
    from parametric import supplement as parametric_supplement
    from bifurcations import supplement as bifurcation_supplement
    output='\n'.join(text).replace('<h2>8. What the competing hypotheses establish here</h2>',supplement()+pendulum_supplement()+relaxation_supplement()+weak_supplement()+parametric_supplement()+bifurcation_supplement()+'<h2>8. What the competing hypotheses establish here</h2>')
    output=output.replace('python phase_locking.py\npython analyze.py','python phase_locking.py\npython dynamics.py\npython analyze.py')
    output=output.replace('python dynamics.py\npython analyze.py','python dynamics.py\npython pendulum.py\npython extend_stress.py\npython analyze.py')
    output=output.replace('python pendulum.py\npython extend_stress.py','python pendulum.py\npython relaxation.py\npython extend_stress.py')
    output=output.replace('python relaxation.py\npython extend_stress.py','python relaxation.py\npython weak_nonlinear.py\npython extend_stress.py')
    output=output.replace('python weak_nonlinear.py\npython extend_stress.py','python weak_nonlinear.py\npython parametric.py\npython extend_stress.py')
    output=output.replace('python parametric.py\npython extend_stress.py','python parametric.py\npython bifurcations.py\npython extend_stress.py')
    (ROOT/'report.html').write_text(output,encoding='utf-8')


def main():
    RESULTS.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True);p=read(ROOT/'protocol.json')
    paths=[q for q in DATA.rglob('*.json') if q.parent.name!='particles' and 'diagnostics' in read(q)]
    rows=[summary(q) for q in paths];rows.sort(key=lambda r:r['name'])
    reg=aggregate([r for r in rows if r['seed'] in p['regime_seeds'] and r['n']==26 and r['dt']==.01],p)
    num=compare_numerics(p);pred=read(RESULTS/'prediction.json');sens=read(RESULTS/'prediction_sensitivity.json')
    checks=[];ref=np.load(DATA/'particles/air_particles_7000.npz')['positions'][-1]
    for suffix in ['half_particle_step','half_fluid_step','fine_grid']:
        z=np.load(DATA/'particles'/f'air_particles_7000_{suffix}.npz');m=read(DATA/'particles'/f'air_particles_7000_{suffix}.json')
        checks.append({'comparison':suffix,'position_rms_change':float(np.sqrt(np.mean(np.sum(wrap(z['positions'][-1]-ref)**2,axis=-1)))),'max_particle_Re':max(m['max_particle_Re_by_species'])})
    csv_write(RESULTS/'all_runs.csv',rows);csv_write(RESULTS/'parameter_regimes.csv',reg);csv_write(RESULTS/'numerical_controls.csv',num);csv_write(RESULTS/'particle_controls.csv',checks)
    particle_rows=[]
    for seed in p['regime_seeds']:particle_rows+=read(DATA/'particles'/f'air_particles_{seed}.json')['rows']
    csv_write(RESULTS/'particle_results.csv',particle_rows)
    plots(p,reg,num,pred,sens);report(p,rows,reg,num,pred,sens,checks)
    save_json(RESULTS/'environment.json',{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform()})
    print(json.dumps({'fluid_configurations':len(rows),'max_budget_residual':max(r['max_budget_residual'] for r in rows),'regime_rows':len(reg),'report':str(ROOT/'report.html')},indent=2))


if __name__=='__main__':main()
