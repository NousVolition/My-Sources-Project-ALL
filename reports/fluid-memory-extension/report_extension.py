"""Readable results and publication figures; refuses incomplete studies."""
import json,html,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from final_experiment import ROOT,P,dump,digest
from analysis import effective_regimes

def load(name):return json.loads((ROOT/'results'/name).read_text())
def fmt(x):return f"{x['mean']:.6g} [{x['ci95'][0]:.6g}, {x['ci95'][1]:.6g}]"
def report():
    summary=load('extension-summary.json');main=load('main-analysis.json');common=load('common-analysis.json');off=load('unforced-analysis.json');controls=load('controls-analysis.json')
    figures=ROOT/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'figure.dpi':150,'savefig.dpi':180,'axes.spines.top':False,'axes.spines.right':False})
    keys=list(main['groups']);fig,ax=plt.subplots(figsize=(10,7));colors=['#336699','#a05000','#277344']
    for j,method in enumerate(P['probe_methods']):
        y=np.arange(len(keys))+(j-1)*.2
        for level,marker,delta in [('coarse','o',-.045),('fine','s',.045)]:
            vals=[main['groups'][k][level][method]['attenuation'] for k in keys];means=np.array([v['mean'] for v in vals]);bounds=np.array([v['ci95'] for v in vals]).T
            ax.errorbar(means,y+delta,xerr=np.abs(bounds-means),fmt=marker,ms=4,capsize=2,color=colors[j],alpha=.7 if level=='coarse' else 1,label=f'{method.replace("_"," ")} ({level})')
    ax.axvline(0,color='0.5',lw=.7);ax.set_yticks(np.arange(len(keys)),[k.replace('forced_nu','nu=').replace('_',' ') for k in keys]);ax.set_xlabel('Retained minus mean scramble history-response RMS (FTLE rate)');ax.legend(fontsize=7);ax.set_title('Whole-flow means and 95% bootstrap intervals; paired resolutions');fig.tight_layout();fig.savefig(figures/'equal-work-layout-effects.png');plt.close(fig)
    keys=list(common['groups']);fig,ax=plt.subplots(figsize=(9,4.7))
    for level,shift,color in [('coarse',-.09,'#8599b1'),('fine',.09,'#184c78')]:
        vals=[common['groups'][k][level]['attenuation'] for k in keys];means=np.array([v['mean'] for v in vals]);bounds=np.array([v['ci95'] for v in vals]).T
        ax.errorbar(means,np.arange(len(keys))+shift,xerr=np.abs(bounds-means),fmt='o',capsize=3,color=color,label=level)
    ax.axvline(0,color='0.5',lw=.7);ax.set_yticks(np.arange(len(keys)),[k.replace('_',' ') for k in keys]);ax.set_xlabel('Retained minus scramble contrast, one identical equal-work impulse');ax.legend();fig.tight_layout();fig.savefig(figures/'common-pulse-effects.png');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(11,3.7))
    for ax,(regime,cfg) in zip(axes,effective_regimes().items()):
        nums=[v for v in main['numerical'] if v['config']['regime']==regime and v['config']['n']==cfg['grids'][-1] and v['config']['layout']=='coincident' and v['config']['dt']==cfg['dt']]
        for v in nums:
            spec=np.array(v['field_stats']['energy_spectrum']);k=np.arange(len(spec));mask=(k>0)&(spec>1e-20);ax.loglog(k[mask],spec[mask],alpha=.65)
        re=[v['field_stats']['Re_lambda'] for v in nums];ax.set_title(f'{regime}\nReλ = {min(re):.1f}–{max(re):.1f}');ax.set_xlabel('Rounded wave-number shell');ax.set_ylabel('Kinetic energy per shell')
    fig.suptitle('Prepared fields: finite-time forced regimes, no stationary inertial-range claim');fig.tight_layout();fig.savefig(figures/'regime-spectra.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4))
    vals=[v for v in common['grid']];names=[f"{v['regime']} / {v['layout']} / {v['seed']}" for v in vals]
    x=np.array([v['pulse_relative'] for v in vals]);ax.plot(np.arange(len(vals)),x,'o',ms=4);ax.set_ylabel('Relative difference of constructed impulse');ax.set_xlabel('Paired flow/layout comparisons');ax.set_title('Grid sensitivity of the common-pulse construction');fig.tight_layout();fig.savefig(figures/'common-pulse-grid-sensitivity.png');plt.close(fig)
    lines=['# Fluid memory extension: equal-work probes and wider flow conditions','',
        'Request ID: `2300ecd1-d208-4695-9b39-62b93d705892`.','',
        'This completed follow-up studies history-conditioned response within ordinary periodic, incompressible, three-dimensional Navier–Stokes. Passive observation points measure local deformation. They are not molecules, and they do not influence the fluid.','',
        f"All declared numerical checks passed: **{summary['all_numerical_checks_passed']}**. Numerical failures and unqualified contrasts are reported below rather than interpreted as physical effects.",'',
        '## What was actually run','',
        'Eight fresh smooth-flow seeds and four fresh forced-flow seeds at each of two viscosities; the forced seeds are paired across viscosities. Three writer/probe layouts, three individual probe methods, four independent phase scrambles, a displacement-matched intact translation and an identical-full-state reset. Every condition has two grids; designated seeds have half steps. The production matrix also contains one common equal-work/equal-strength impulse, forcing-off readouts, translation invariance checks, and finer moderate-regime audits. Numerical pilot seeds are excluded from inference.','',
        'The CPU solver independently implements padded conservative flux divergence, full complex Fourier modes, pressure projection and adaptive DOP853. It restarts the exact saved initial states for twelve comparisons, including tighter tolerances and full response windows. This is a second implementation by the same author, not replication by an independent laboratory. It does not repeat every spinup or every production trajectory.','',
        '## Operational memory measure','',
        'Two opposite localized preparation forces create different histories. At readout, the complete low Fourier modes and high modal powers are matched, but the full vector/phase fields remain different. For each history, response is probed FTLE minus unprobed FTLE at the same ten passive centers. The history contrast is the RMS difference across these centers. Attenuation is retained contrast minus the mean of four phase-scrambled contrasts. These definitions test response dependence on retained state structure; they do not add memory beyond the present complete Navier–Stokes state.','',
        'A signed physical-effect claim requires confidence intervals excluding zero on both grids with the same sign, average absolute per-flow grid change below 10% of the fine-grid mean effect, and the numerical controls. Positive attenuation alone does not establish a unique history carrier: the intact-translation comparison tests whether a large spatial rearrangement can produce the same attenuation.','',
        '## Results: individual probes','',
        '| Regime / layout | Probe | Fine-grid attenuation [95% CI] | Resolution-qualified signed effect |','|---|---|---:|---|']
    for key,v in main['groups'].items():
        for method in P['probe_methods']:lines.append(f"| {key} | {method} | {fmt(v['fine'][method]['attenuation'])} | {main['claims'][key][method]['physical_claim_passed']} |")
    lines+=['','![Equal-work and layout results](figures/equal-work-layout-effects.png)','','## Results: one identical impulse for all histories and arms','','| Regime / layout | Fine-grid attenuation [95% CI] | Translation-minus-scramble specificity [95% CI] | Resolution-qualified signed effect |','|---|---:|---:|---|']
    for key,v in common['groups'].items():lines.append(f"| {key} | {fmt(v['fine']['attenuation'])} | {fmt(v['fine']['specificity'])} | {common['claims'][key]['physical_claim_passed']} |")
    lines+=['','![Common-pulse results](figures/common-pulse-effects.png)','','## Numerical checks and implementation correction','','| Check group | Failed checks |','|---|---|']
    for key,v in summary['checks'].items():lines.append(f"| {key} | {', '.join(k for k,passed in v.items() if not passed) or 'None'} |")
    lines+=['',
        'The initial implementation failed lower-viscosity grid tests. A diagnosed reset-control bug used slightly non-integer FFT wave numbers; at exactly parallel integer wave vectors, normalizing an almost-zero cross product amplified rounding noise instead of choosing the defined alternate vector. Corrected wave vectors are exact integers and the strict cutoff is computed with integer arithmetic. Original sources, failed pilots and failed qualification are preserved. New corrected spinups and all corrected pilot comparisons were run; thresholds were never relaxed. Two regression tests specifically cover this degeneracy, alongside analytic shear, independent RHS, equal-work, modal matching and direct-observation checks.','',
        'The forcing-off comparisons use the same complete readout states and impulses, then remove sustaining forcing. They involve only two independent flows per viscosity and are exploratory. Their intervals should not be interpreted as strong evidence of general forcing independence.','',
        '## Forcing-off comparisons','',
        '| Regime | Probe family / method | Off-minus-on attenuation [95% CI], n=2 |','|---|---|---:|']
    for regime,v in off['regimes'].items():
        for family,groups in v['groups'].items():
            for method,value in groups['paired_off_minus_on'].items():lines.append(f'| {regime} | {family} / {method} | {fmt(value)} |')
    lines+=['','## Probe matching and remaining limits','',
        'Fixed-RMS probes use the same shape and strength but generally perform different work on different states. Equal-work/fixed-shape probes vary amplitude and can be much larger than the nominal weak probe. Individual equal-work/equal-RMS probes vary direction by state. The stricter common-pulse experiment uses one identical impulse across all histories and interventions within a comparison, chosen before future outcomes from a localized divergence-free basis and the nullspace of work constraints. Its geometry can differ across flows, layouts and grids. Complete coefficients, angles, ranks, singular values and measured spatial extent are recorded. A localized basis alone does not guarantee that the constructed combination retains the original probe’s localization.','',
        '![Common-pulse construction sensitivity](figures/common-pulse-grid-sensitivity.png)','',
        '![Prepared-flow spectra](figures/regime-spectra.png)','',
        'These are moderate, finite-time forced regimes. A fixed eight-unit spinup is not proof of stationary, fully developed or high-Reynolds-number turbulence. Positive viscosity, energy budgets, incompressibility, CFL, spectral tails, endpoint Kolmogorov resolution, step/grid convergence and solver agreement are screened. Passing a screen does not prove DNS accuracy in every observable. The box size and statistical population remain limited. Bootstrap intervals resample independent flow seeds; markers, layouts and resolutions are not extra samples. Secondary intervals are unadjusted, with only four forced flows. Previous results informed this follow-up, so it is not independent confirmation of the original selection.','',
        'No persistent leader, learning, causal hierarchy, heteroclinic network, isotope effect, molecular water memory, evolutionary change or new Navier–Stokes physics follows from these tests. Complete full-state matching removes deterministic history dependence by construction; agreement of the reset histories is a consistency control, not evidence of a physical erasure mechanism.','',
        '## Reproduction and provenance','',
        'Install dependencies from `requirements.txt`. Run `python reproduce_extension.py --analysis-only` on restored recordings, or `python reproduce_extension.py --run` to perform the full declared study. The latter requires a supported NVIDIA GPU and substantial CPU time/disk space. Raw NPZ containers are preserved without reencoding in the release; verify and restore them with `restore_recordings.py`. The manifest records each container and each NPY member SHA256. Sources, protocols, numerical amendments, execution counts, raw arrays, independent results and statistical per-flow values are included.','',
        'Method references: [SciPy DOP853 documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html), [NASA grid-convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).','']
    text='\n'.join(lines);(ROOT/'README.md').write_text(text,encoding='utf8',newline='\n')
    # Self-contained report content remains simple; README is the authoritative table.
    parts=[]
    for line in lines:
        if line.startswith('!['):
            alt,src=line[2:].split('](',1);parts.append(f'<img alt="{html.escape(alt)}" src="{html.escape(src[:-1])}">')
        elif line.startswith('#'):level=min(len(line)-len(line.lstrip('#')),3);parts.append(f'<h{level}>{html.escape(line[level:].strip())}</h{level}>')
        elif line:parts.append('<p>'+html.escape(line)+'</p>')
    page='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Fluid memory extension</title><style>body{max-width:1100px;margin:36px auto;padding:0 22px;font:16px/1.55 system-ui;color:#183047;background:#fbfcfd}img{max-width:100%;background:white}h1,h2{line-height:1.2}p{overflow-wrap:anywhere}</style><body>'+''.join(parts)+'</body></html>'
    (ROOT/'report.html').write_text(page,encoding='utf8',newline='\n')
    dump(ROOT/'results'/'report-provenance.json',dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,readme_sha256=digest(ROOT/'README.md'),report_sha256=digest(ROOT/'report.html'),sources={f:digest(ROOT/f) for f in ['report_extension.py','followup_analysis.py']},figures={p.name:digest(p) for p in figures.glob('*.png')}))
    print('Created report and four figures',flush=True)
if __name__=='__main__':report()
