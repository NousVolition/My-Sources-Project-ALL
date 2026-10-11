"""Readable scientific report and standalone figures from completed results."""
import base64,html,json,platform
from pathlib import Path
import numpy as np
import scipy,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import ROOT,P,dump,digest
from analyze import interval,load_run,branch_index


def read(name):return json.loads((ROOT/'results'/name).read_text())
def num(x):return f'{x:.6g}'
def ci(r):return f"{num(r['mean'])} [{num(r['ci95'][0])}, {num(r['ci95'][1])}]"
def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.svg',bbox_inches='tight');plt.close(fig)
def error(ax,x,r,**kw):
    lo,hi=r['ci95'];ax.errorbar(x,r['mean'],yerr=[[max(0,r['mean']-lo)],[max(0,hi-r['mean'])]],capsize=4,**kw)


def main():
    (ROOT/'figures').mkdir(exist_ok=True)
    c=read('causal-results.json');pred=read('prediction-results.json');audit=read('numerical-audit.json')
    matched=read('matched-sham-results.json');delay=read('delayed-probe-results.json');aux=read('followup-numerical-audit.json');verify=read('verification.json')
    primary=c['primary'];rows=[r for r in c['per_flow'] if r['horizon']==P['primary_horizon']]
    specific=matched['horizons']['0.4']['matched_minus_scrambled']
    specificity_error=float(np.mean([abs(r['matched_minus_scrambled_grid_difference']) for r in matched['numerical']]))
    specific_supported=bool(specific['ci95'][0]>0 and aux['gates']['matched_grid'] and aux['gates']['matched_half_step'] and specificity_error<.1*abs(specific['mean']))
    values=np.array(primary['per_flow']);retained=np.array([r['retained_history_contrast'] for r in rows])
    rng=np.random.default_rng(9140);ix=rng.integers(0,len(values),(4000,len(values)))
    dz=values.mean()/values.std(ddof=1);bootdz=values[ix].mean(1)/values[ix].std(1,ddof=1)
    relative=100*values.mean()/retained.mean();bootrelative=100*values[ix].mean(1)/retained[ix].mean(1)
    effects=dict(paired_standardized_effect_dz=float(dz),paired_standardized_effect_dz_ci95=np.quantile(bootdz,[.025,.975]).tolist(),
                 relative_attenuation_percent=float(relative),relative_attenuation_percent_ci95=np.quantile(bootrelative,[.025,.975]).tolist(),
                 numerical_error_fraction_of_mean=float(c['finest_grid_effect_change_mean_absolute']/max(abs(primary['mean']),1e-30)))
    dump(ROOT/'results/effect-sizes.json',effects)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    from solver import Solver
    sample,meta=load_run(ROOT/'data',P['test_seeds'][0],keys=['start_fields'])
    solver=Solver(32,P['nu']);selected=[branch_index(meta,arm,h,0) for h in range(2) for arm in ['retained','scramble0']]
    omega=solver.real(solver.curl(solver.unpack(sample['start_fields'][selected])))
    vorticity=np.sqrt(np.sum(omega**2,axis=1))[:,:,:,16]
    fig,axs=plt.subplots(2,2,figsize=(8,7),layout='constrained')
    for i,ax in enumerate(axs.flat):
        im=ax.imshow(vorticity[i].T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],vmin=0,vmax=vorticity.max(),cmap='viridis')
        ax.set(title=('Structure retained' if i%2==0 else 'Phase scramble 0')+(' · + history' if i<2 else ' · − history'),xlabel='x',ylabel='y')
    fig.colorbar(im,ax=axs,label='Vorticity magnitude at z = π',shrink=.85)
    fig.suptitle('Single fresh flow, immediately after intervention\nSame coarse coefficients and scalar power at every Fourier mode',fontsize=12)
    save(fig,'structure-example')
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    x=np.arange(len(rows));axs[0].plot(x,[r['retained_history_contrast'] for r in rows],'o-',label='Structure retained',color='#245a81')
    axs[0].plot(x,[r['scrambled_history_contrast'] for r in rows],'s-',label='Mean of 4 phase scrambles',color='#ac493b')
    axs[0].set(xticks=x,xticklabels=[str(r['seed'])[-2:] for r in rows],xlabel='Fresh test-flow seed suffix',ylabel='RMS history contrast in probe-induced FTLE',title='Identical probe, matched coarse field and modal power');axs[0].legend(fontsize=9)
    for i,t in enumerate(P['horizons']):
        r=c['horizons'][str(t)]['primary_attenuation'];error(axs[1],t,r,fmt='o',color='#245a81')
        axs[1].scatter(np.full(len(values),t)+np.linspace(-.012,.012,len(values)),r['per_flow'],s=13,alpha=.45,color='#245a81')
    axs[1].axhline(0,color='grey',lw=1);axs[1].set(xlabel='Forward observation horizon (dimensionless)',ylabel='Retained minus scrambled contrast',title='Positive means attenuation; bars: whole-flow 95% CI')
    save(fig,'erasure')
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    controls=['current','capacity_current','shuffled_history','across_flow_history'];labels=['Current quadratic','Equal-count current','Within-flow shuffled','Across-flow shuffled']
    for ax,target in zip(axs,P['prediction']['targets']):
        record=pred[f'{target}:T0.4']
        for i,control in enumerate(controls):error(ax,i,record['history_minus_'+control],fmt='o',color='#245a81')
        ax.axhline(0,color='grey',lw=1);ax.set(xticks=range(4),xticklabels=labels,ylabel='History model RMSE minus control RMSE',title=target.replace('native ','').replace('forward ','')+'\nNegative favors history')
        ax.tick_params(axis='x',rotation=22)
    save(fig,'prediction')
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for t in P['horizons']:error(axs[0],t,matched['horizons'][str(t)]['matched_minus_scrambled'],fmt='o',color='#447054')
    axs[0].axhline(0,color='grey',lw=1);axs[0].set(xlabel='Forward horizon',ylabel='Matched translation minus scrambled contrast',title='Equal-size displacement control (12 flows)')
    for wait in [0.,.4,.8]:
        r=interval([r['primary_attenuation'] for r in rows if r['seed'] in P['long_horizon_seeds']]) if wait==0 else delay['horizons'][str(wait)]['0.4']['primary_attenuation']
        error(axs[1],wait,r,fmt='o',color='#9a651f')
    axs[1].axhline(0,color='grey',lw=1);axs[1].set(xlabel='Unforced waiting time before probe',ylabel='Retained minus scrambled contrast, horizon 0.4',title='Later disturbance: same 4 flows (exploratory)')
    save(fig,'controls-and-waiting')
    fig,axs=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
    for ax,key,limit,title in zip(axs,['ftle_rmse','response_ftle_rmse','max_probe_response_field_relative_L2'],
        ['finest_pair_ftle_rmse','finest_pair_response_ftle_rmse','finest_pair_probe_response_field_relative_L2'],['FTLE grid difference','Probe-response FTLE difference','Full response-field relative L2']):
        for j,(coarse,fine) in enumerate([(24,32),(32,48)]):
            v=[r[key] for r in audit['grid'] if r['coarse_n']==coarse];ax.scatter(np.full(len(v),j)+np.linspace(-.07,.07,len(v)),v,s=20)
        ax.axhline(P['numerical_gates'][limit],color='#ac493b',ls='--',label='Finest-pair gate');ax.set(xticks=[0,1],xticklabels=['24 → 32','32 → 48'],yscale='log',title=title);ax.legend(fontsize=8)
    save(fig,'numerical-convergence')
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained');duration={};strength={};viscosity={}
    for t in P['horizons']+P['long_horizons']:
        v=[r['primary_attenuation'] for r in c['sensitivity'] if r['config']['end']==1.6 and r['horizon']==t]
        duration[str(t)]=interval(v);error(axs[0],t,duration[str(t)],fmt='o',color='#245a81')
    axs[0].axhline(0,color='grey',lw=1);axs[0].set(xlabel='Forward observation horizon',ylabel='Retained minus scrambled contrast',title='Longer observation · same 4 flows')
    for amplitude in [.01,.02]:
        v=[r['primary_attenuation']/amplitude for r in c['sensitivity'] if r['config']['amplitude']==amplitude and r['config']['nu']==.12 and r['config']['end']==.8 and r['horizon']==.4] if amplitude==.01 else [r['primary_attenuation']/amplitude for r in rows if r['seed'] in P['amplitude_control_seeds']]
        strength[str(amplitude)]=interval(v);error(axs[1],amplitude,strength[str(amplitude)],fmt='o',color='#245a81')
    axs[1].set(xticks=[.01,.02],xlabel='Probe RMS velocity',ylabel='Contrast divided by probe RMS',title='Probe strength · same 3 flows')
    for nu in [.08,.12,.2]:
        v=[r['primary_attenuation'] for r in c['sensitivity'] if r['config']['nu']==nu and r['horizon']==.4] if nu!=.12 else [r['primary_attenuation'] for r in rows if r['seed'] in P['viscosity_control_seeds']]
        viscosity[str(nu)]=interval(v);error(axs[2],nu,viscosity[str(nu)],fmt='o',color='#245a81')
    axs[2].set(xticks=[.08,.12,.2],xlabel='Positive viscosity ν',ylabel='Retained minus scrambled contrast',title='Viscosity · same 3 flows')
    fig.suptitle('Secondary sensitivity controls; unadjusted and incompletely refined',fontsize=12)
    save(fig,'sensitivity');dump(ROOT/'results/sensitivity-summary.json',dict(duration=duration,probe_normalized=strength,viscosity=viscosity))
    # Same four seeds for waiting comparisons, avoiding an ensemble-size contrast.
    waits=[0.,.4,.8];paired={}
    seeds=P['long_horizon_seeds']
    for w in waits:
        v=[r['primary_attenuation'] for r in rows if r['seed'] in seeds] if w==0 else [r['primary_attenuation'] for r in delay['per_flow'] if r['horizon']==.4 and r['config']['n']==32 and r['config']['dt']==.01 and r['config']['wait']==w]
        paired[str(w)]=interval(v)
    dump(ROOT/'results/same-flow-waiting-summary.json',paired)
    # Environment is actual execution context, not a requirement that other hosts match it.
    env=dict(python=platform.python_version(),platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__)
    try:
        import cupy
        env.update(cupy=cupy.__version__,gpu=cupy.cuda.runtime.getDeviceProperties(0)['name'].decode(),cuda_runtime=cupy.cuda.runtime.runtimeGetVersion())
    except Exception:pass
    dump(ROOT/'results/environment.json',env)
    attenuation='supported' if c['attenuation_supported'] else 'not established'
    direction='attenuated' if primary['mean']>0 else 'increased'
    statement=f"Phase scrambling {direction} the mean history-conditioned probe-response contrast in this ensemble. The primary operational attenuation criterion is {attenuation}. Additional attenuation beyond a strength-matched structured displacement is {'supported in this exploratory comparison' if specific_supported else 'not established'}."
    dump(ROOT/'results/specificity-summary.json',dict(matched_minus_scrambled=specific,mean_absolute_fine_grid_change=specificity_error,additional_attenuation_supported=specific_supported,status='Exploratory, unadjusted; translation also relocates the structure relative to the fixed probe.'))
    history_summary=[]
    for target in P['prediction']['targets']:
        r=pred[f'{target}:T0.4'];history_summary.append(f"{target}: history-versus-current RMSE difference {ci(r['history_minus_current'])}; relative RMSE benefit {r['relative_rmse_benefit_percent']:.2f}%; equal-capacity difference {ci(r['history_minus_capacity_current'])}.")
    count=audit['count']+len(aux['verified_auxiliary_recordings'])+6
    method="""The fluid obeys ∂u/∂t + (u·∇)u = −∇p + ν∇²u + f, ∇·u = 0 in a periodic box [0,2π)³. Pressure projection, a rotational nonlinearity and strict 2/3 Fourier truncation enforce the discrete equations. Float64/complex128 synchronous RK4 advances the field, passive coordinates dx/dt = u(x,t), and tangent matrices dF/dt = ∇u(x,t)F. Velocity and gradients at markers are evaluated directly from the retained Fourier polynomial. FTLE = log σmax(F)/T is a local forward deformation rate, not a molecular bond or leader measure. All times and fields are dimensionless. ν = 0.12 and initial RMS velocity = 0.6; this is a smooth, modest-Reynolds-number regime, not a turbulence universality survey."""
    preparation="""Forty fresh analytic random starting flows are split into 20 training, 8 validation and 12 held-out test seeds; two additional numerical pilot seeds are excluded from inference. Two preparations receive opposite divergence-free high-band forcing until t = 0.2, followed by unforced evolution to t = 0.4. Causal arms then share the complete low-mode field (|k|∞ ≤ 2) and exactly the same scalar power at every high Fourier mode. Each high mode can retain a different vector polarization and phase structure. The fixed probe is a divergence-free velocity impulse of RMS 0.02. Paired unprobed branches remove background evolution. Ten common observation centers and their finite neighbors are passive measurements; they do not exert force."""
    controls_text="""Four seeded high-mode phase scrambles target spatial phase organization while preserving incompressibility, coarse coefficients, energy, enstrophy and every scalar modal power at erasure. A small translation is a sham. A second translation preserves the high-band pattern and matches the mean initial L2 displacement caused by scrambling; it also changes position relative to the probe and background, so it is an imperfect specificity control. A common complete-state replacement makes both preparation states identical; identical subsequent responses are expected deterministically and provide a control, not evidence for a new memory law. Delayed-probe runs wait 0.4 or 0.8 before applying the disturbance. Longer-response runs observe out to 1.6, a different question."""
    controls_text+=' The same velocity impulse need not do the same work: its energy increment is the velocity–impulse inner product plus half the impulse norm squared. A separate post-result work audit is saved in impulse-energy-audit.json. Phase-dependent work is a possible response pathway; the study does not isolate it. Equal-work disturbances and alternate writer/probe positions remain unrun.'
    inference="""The primary outcome at forward horizon 0.4 is, per independent test flow, the RMS over common centers of the difference between the two preparations' probe-minus-unprobed FTLE responses, retained minus the mean of four scrambles. Positive values mean attenuation by scrambling. Means and 95% percentile intervals resample whole flow seeds 4,000 times; centers, preparations, grids and scramble replicates are not additional independent samples. There is one prespecified primary contrast and many unadjusted exploratory follow-ups. Twelve test seeds give limited precision and the bootstrap is not a guarantee of nominal small-sample coverage. Paired effect size dz uses the between-flow standard deviation of the signed contrast. Protocol files were written locally before scientific effects were examined; there was no external registration or blind third-party audit."""
    prediction_text="""Current observations contain exact local velocity, the full local gradient, strain/rotation summaries, finite-neighbor geometry and instantaneous separation rates. Past features add distance and angle changes at lookbacks 0.1, 0.2 and 0.3 and cannot access future fields or tangent matrices. Quadratic ridge current features (560 inputs) are compared with the same features plus 48 history inputs. Controls add 48 current cubic features or shuffle histories within or across flows. Standardization uses training flows; regularization is selected only on validation flows. Thirty models were frozen before opening test observations. Frozen models are also evaluated on finer-grid test data. Every model observes an incomplete local present state; better prediction with past data would support useful partial-observation history, not fundamental non-Markovian Navier–Stokes dynamics."""
    numerics=f"""An analytic divergence-free shear benchmark checks arbitrary-location velocities/gradients, tracer paths, tangent evolution, energy dissipation and fourth-order time convergence. CPU/GPU pilot comparisons agree near roundoff. The former interpolation issue is removed by direct spectral evaluation. Pilot numerical errors triggered a documented promotion from 24³ to 32³ before production effects were examined. All twelve test flows have 24³/32³/48³ controls, all have half-step 32³ controls, and two have half-step 48³ controls. Production numerical gates: {audit['gates']}. Additional-control gates: {aux['gates']}. Mean absolute change in the primary contrast from 32³ to 48³ is {num(c['finest_grid_effect_change_mean_absolute'])}, compared with mean effect {num(primary['mean'])}; the 10%-of-effect gate is {c['numerical_effect_gate']}. Integrity/analytic/relabeling verification passes: {verify['all_passed']} ({verify['analytic_and_analysis_tests']} unit tests). Gates at the primary horizon do not prove global continuum convergence or validate all extended durations. Numerical intervals are assessed separately from sampling intervals."""
    limitations="""Changing a present velocity field can change its future under ordinary Navier–Stokes. This experiment targets how deliberately written history is represented in the current field and modifies response under limited matching. It does not compare identical complete states with different physical histories. Phase scrambling leaves vector polarization relations within each Fourier mode, need not erase every carrier, is an external instantaneous intervention, and can create other structures. The experiment cannot uniquely identify a carrier or transfer Josephson hysteresis, heteroclinic memory, molecular isotope effects, learning, leadership, evolutionary changes in water or new Navier–Stokes physics. No isotope masses or actual molecules are simulated. Fresh flows come from one ensemble family and use one implementation, although CPU/GPU and analytic controls are included. Larger domains, forced turbulent regimes, other writer/probe locations and independent solvers remain untested. The four-flow waiting/long-duration and three-flow viscosity/probe-amplitude controls are exploratory; two-flow delayed resolution checks do not cover every waiting branch."""
    reproducibility="""Install requirements.txt with Python 3.12 (CPU) or add requirements-gpu.txt for a compatible CUDA GPU. To reanalyse the exact recordings, download the code/results ZIP and raw-recordings.json plus every numbered raw part from the GitHub release into a fresh folder, then run `python restore_recordings.py --parts DOWNLOAD_FOLDER --out STUDY_FOLDER`, `python analyze.py`, `python followup_analysis.py`, `python template_analysis.py`, `python impulse_energy.py`, `python verify.py`, and `python report.py` from the study folder. Default archived data use GPU filenames; reanalysis requires no GPU calculation. To rerun the entire design, run `python reproduce.py --destination NEW_EMPTY_PATH --backend cpu` (or gpu). The CPU rerun is implemented but the complete experiment was actually run on GPU; CPU/GPU equivalence was tested on two short pilots. No complete CPU rerun is claimed. The split tar.gz archive copies the original compressed NPZ containers byte for byte, preserving every NPY member and all float64/complex128 values. Checksums, original container hashes, execution sources, protocols and exact fitting-source snapshots are retained. After restoration, verify all source/data hashes before interpreting results. Published raw checkpoints permit every response and feature calculation to be repeated without rerunning the fluid."""
    numerics+=' The smaller-step controls keep the endpoint fixed. Separate longer-window runs extend forward observation to 1.6; their earlier trajectory/response prefixes agree exactly with the short runs for all four checked flows. Delayed probes test waiting before the disturbance separately.'
    sections=[('Finding',statement+f" Primary retained-minus-scrambled contrast: {ci(primary)}. Relative attenuation: {relative:.2f}%; paired dz: {dz:.3f}."),
        ('What was completed',f'{count} simulation configurations including 6 numerical pilots; 40 fresh production flows, with 12 independent held-out flows for the primary result. Ninety-one base/refinement/sensitivity configurations, 26 matched-displacement configurations and 12 delayed-probe configurations. Full fields, marker trajectories, tangent matrices, numerical diagnostics and frozen held-out predictions are saved.'),
        ('History prediction','At horizon 0.4, history improves prediction of ordinary forward deformation by 7.36% over the present-observation model, with an interval favoring history and agreement on the finer grid. It does not improve prediction of the probe-induced change: its error is 0.77% higher. These are different questions; the positive erasure result cannot be substituted for a predictive-response benefit. '+' '.join(history_summary)),('Fluid and measurement',method),('Writing and matching histories',preparation),('Erasure and controls',controls_text),
        ('Inference',inference),('Frozen prediction comparison',prediction_text),('Numerical verification',numerics),('Limits and unanswered questions',limitations),('Reproduce and inspect',reproducibility)]
    link='https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-memory-erasure'
    release='https://github.com/NousVolition/My-Sources-Project-ALL/releases/tag/fluid-memory-erasure-2026-10-10'
    md='# Fluid memory: fresh flows, structure erasure and delayed probes\n\n'+''.join('## '+title+'\n\n'+text+'\n\n' for title,text in sections)
    md+='## Quantitative controls\n\n| Forward horizon | Retained − scrambled (95% CI) | Matched translation − scrambled (95% CI) |\n|---|---|---|\n'
    for t in P['horizons']:md+=f"| {t} | {ci(c['horizons'][str(t)]['primary_attenuation'])} | {ci(matched['horizons'][str(t)]['matched_minus_scrambled'])} |\n"
    md+='\n| Wait before probe | Primary contrast, same four flows (95% CI) |\n|---|---|\n'
    for w,r in paired.items():md+=f'| {w} | {ci(r)} |\n'
    for name in ['erasure','prediction','controls-and-waiting','numerical-convergence','sensitivity','structure-example']:md+=f'\n![{name}](figures/{name}.png)\n'
    md+=f'\n[GitHub study]({link}) · [Raw recordings and release]({release})\n\nNumerical verification context: [NASA grid convergence tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html). Three-grid/time checks are performed here; Richardson extrapolation is not asserted for this nonlinear trajectory ensemble.\n'
    (ROOT/'README.md').write_text(md,encoding='utf-8')
    images=''
    for name in ['erasure','prediction','controls-and-waiting','numerical-convergence','sensitivity','structure-example']:
        b64=base64.b64encode((ROOT/'figures'/f'{name}.png').read_bytes()).decode();images+=f'<figure><img alt="{name}" src="data:image/png;base64,{b64}"></figure>'
    content=''.join(f'<section><h2>{html.escape(title)}</h2><p>{html.escape(text)}</p></section>' for title,text in sections)
    tables=f'<h2>Matched displacement and waiting</h2><table><tr><th>Horizon</th><th>Matched translation − scrambled (95% CI)</th></tr>'+''.join(f"<tr><td>{t}</td><td>{ci(matched['horizons'][str(t)]['matched_minus_scrambled'])}</td></tr>" for t in P['horizons'])+'</table>'
    tables+='<table><tr><th>Waiting time</th><th>Same four flows: retained − scrambled (95% CI)</th></tr>'+''.join(f'<tr><td>{w}</td><td>{ci(r)}</td></tr>' for w,r in paired.items())+'</table>'
    doc=f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Fluid memory erasure results</title><style>body{{font:17px/1.6 system-ui,sans-serif;color:#172c3a;max-width:1120px;margin:36px auto;padding:0 24px;background:#fafbf9}}h1{{font-size:36px;line-height:1.2}}h2{{font-size:23px;margin-top:32px}}p{{max-width:1000px}}img{{width:100%;height:auto}}figure{{margin:28px 0;background:white;padding:12px}}table{{border-collapse:collapse;width:100%;margin:20px 0}}td,th{{padding:12px;border-bottom:1px solid #b9c6c9;text-align:left}}a{{color:#245a81}}.intro{{padding:20px;background:#e9f0f3;border-left:5px solid #245a81}}</style><body><h1>Fluid memory: fresh flows and structure erasure</h1><p>Completed Navier–Stokes experiment · dimensionless · 12 held-out flows · all uncertainty is across independent starting flows</p><div class="intro"><strong>{html.escape(statement)}</strong><p>Primary contrast: {ci(primary)}. Numerical screens: {audit["all_passed"]}. Effect-to-numerical-error gate: {c["numerical_effect_gate"]}.</p></div>{images}{tables}{content}<p><a href="{link}">GitHub code and results</a> · <a href="{release}">Raw recordings</a> · <a href="https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html">NASA numerical verification context</a></p></body></html>'
    (ROOT/'report.html').write_text(doc,encoding='utf-8')
    dump(ROOT/'results/report-provenance.json',dict(report_source_sha256=digest(ROOT/'report.py'),inputs={p.name:digest(p) for p in (ROOT/'results').glob('*.json') if p.name!='report-provenance.json'},figure_files=[p.name for p in (ROOT/'figures').glob('*')]))
    print(statement+'\nPrimary '+ci(primary)+'\n'+'\n'.join(history_summary))


if __name__=='__main__':main()
