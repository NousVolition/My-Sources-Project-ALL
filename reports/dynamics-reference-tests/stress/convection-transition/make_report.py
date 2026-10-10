"""Generate the readable report and figures from recorded results only."""
from pathlib import Path
import csv
import html
import json
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from convection import ROOT, RAC, Q

DATA=ROOT/'data'
FIG=ROOT/'figures'


def next_steps_content():
    """Render planned work separately from the completed numerical evidence."""
    plan=json.loads((ROOT/'next-study-plan.json').read_text(encoding='utf-8'))
    assert plan['completed_new_simulations']==0
    esc=lambda value: html.escape(value,quote=True)
    md=f"## {plan['title']}\n\n**{plan['status']} — no new fluid or pond runs accompany this plan.**\n\n{plan['summary']}\n\n"
    page=f'<section id="next-steps"><h2>{esc(plan["title"])}</h2><p class="tag">{esc(plan["status"])}</p><p>No new fluid or pond runs accompany this plan.</p><p>{esc(plan["summary"])}</p>'
    md+='| Next test | Primary outcome | Go/no-go decision |\n| --- | --- | --- |\n'
    page+='<div style="overflow-x:auto"><table><tr><th>Next test</th><th>Primary outcome</th><th>Go/no-go decision</th></tr>'
    for row in plan['stages']:
        md+=f"| {row['stage']} | {row['measure']} | {row['gate']} |\n"
        page+=f'<tr><td>{esc(row["stage"])}</td><td>{esc(row["measure"])}</td><td>{esc(row["gate"])}</td></tr>'
    md+='\n'
    page+='</table></div>'
    for row in plan['stages']:
        md+=f"### {row['stage']}\n\n{row['change']}\n\n"
        page+=f'<h3>{esc(row["stage"])}</h3><p>{esc(row["change"])}</p>'
    for title,key in [('Where the energy could come from','energy_explanation'),('What the present pond equations predict','pond_identity')]:
        md+=f"### {title}\n\n{plan[key]}\n\n"
        page+=f'<h3>{title}</h3><p>{esc(plan[key])}</p>'
    md+='### Controls and acceptance criteria\n\n'
    page+='<h3>Controls and acceptance criteria</h3><ul>'
    for control in plan['controls']:
        md+=f'- {control}\n'
        page+=f'<li>{esc(control)}</li>'
    page+='</ul>'
    md+='\n'
    for title,key in [('How long to run','duration'),('Decision after these tests','decision')]:
        md+=f"### {title}\n\n{plan[key]}\n\n"
        page+=f'<h3>{title}</h3><p>{esc(plan[key])}</p>'
    refs=[dict(label='Machine-readable next-study plan',url='next-study-plan.json')]+plan['references']
    md+=' · '.join(f'[{r["label"]}]({r["url"]})' for r in refs)+'\n\n'
    page+='<p>'+' · '.join(f'<a href="{esc(r["url"])}">{esc(r["label"])}</a>' for r in refs)+'</p></section>'
    return md,page


def main(reuse_figures=False):
    FIG.mkdir(exist_ok=True)
    s=json.loads((DATA/'summary.json').read_text())
    checks=json.loads((DATA/'checks.json').read_text())
    inertia=json.loads((DATA/'inertia.json').read_text())
    with (DATA/'cases.csv').open(newline='') as f:
        cases=list(csv.DictReader(f))
    by={r['label']:r for r in cases}
    tests=list(ET.parse(DATA/'independent-tests.xml').getroot().iter('testcase'))
    failures=sum(t.find('failure') is not None or t.find('error') is not None for t in tests)
    groups=s['primary_groups']
    next_md,next_html=next_steps_content()
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        if reuse_figures:
            for extension in ['png','svg']:
                assert (FIG/(name+'.'+extension)).is_file(), 'Missing figure to reuse'
            plt.close(fig)
            return
        for extension in ['png','svg']:
            fig.savefig(FIG/(name+'.'+extension),dpi=160,bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    ratios=np.linspace(.78,1.22,250)
    axes[0].plot(ratios,Q*(np.sqrt(ratios)-1),color='#758791',label='Linear theory')
    groups=s['primary_groups']
    axes[0].scatter([g['ratio'] for g in groups],[g['mean_rate'] for g in groups],color='#087f8c',label='Full PDE, refined step',zorder=3)
    axes[0].axhline(0,color='#333',lw=.8)
    axes[0].axvline(1,color='#333',ls=':',lw=.8)
    axes[0].set(xlabel='Heating / theoretical critical heating',ylabel='Disturbance amplitude growth rate',title='A physical change reverses the sign')
    axes[0].legend(fontsize=8)
    for ratio,color in [(.95,'#087f8c'),(1.05,'#ba581b')]:
        with np.load(DATA/f'refined_primary_r{ratio:g}_s410.npz') as a:
            r=a['records'];axes[1].semilogy(r[:,0],r[:,1]/r[0,1],color=color,label=f'{ratio:g} × threshold')
    axes[1].axhline(1,color='#888',lw=.7)
    axes[1].set(xlabel='Time / thermal diffusion time',ylabel='Disturbance / initial disturbance',title='Matched initial velocity disturbance')
    axes[1].legend(fontsize=8)
    for label,text,color in [('refined_primary_r1.05_s410','Accepted step 0.005','#087f8c'),('coarse_negative_control','Rejected step 0.2','#a63829')]:
        with np.load(DATA/(label+'.npz')) as a:
            r=a['records'];axes[2].semilogy(r[:,0],r[:,1]/r[0,1],color=color,label=text)
    axes[2].set(xlabel='Time / thermal diffusion time',ylabel='Disturbance / initial disturbance',title='A coarse step hides real growth')
    axes[2].legend(fontsize=8)
    fig.suptitle('Ideal heated fluid layer · synthetic Boussinesq computations',fontsize=13)
    save(fig,'transition')

    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    with np.load(DATA/'inertia_r-4_pr1.npz') as a:
        pick=a['t']<=.35
        axes[0].plot(a['t'][pick],a['full'][pick,1],label='Full inertia',color='#087f8c')
        axes[0].plot(a['t'][pick],a['stokes'][pick,1],label='Inertia removed',color='#ba581b',ls='--')
    axes[0].axhline(0,color='#777',lw=.8)
    axes[0].set(xlabel='Time / thermal diffusion time',ylabel='Temperature-mode amplitude / initial',title='Stable stratification can still oscillate')
    axes[0].legend()
    stable=[r for r in inertia if r['ratio']==-4]
    axes[1].loglog([r['prandtl'] for r in stable],[100*r['theta_relative_curve_error'] for r in stable],'o-',color='#087f8c')
    axes[1].set(xlabel='Prandtl number',ylabel='Stokes temperature-curve error (%)',title='Removing inertia improves in its own limit')
    fig.suptitle('Separate linear modal comparison · no additional full PDE runs',fontsize=13)
    save(fig,'inertia')

    below=next(g for g in groups if g['ratio']==.95)
    above=next(g for g in groups if g['ratio']==1.05)
    fine=[r for r in cases if float(r['dt'])<=.005]
    budget=max(max(float(r['relative_kinetic_budget']),float(r['relative_thermal_budget'])) for r in fine)
    failed_coarse=[r['label'] for r in cases if r['energy_budgets_passed']=='False']
    rows='\n'.join(f"| {g['ratio']:.2f} | {g['rayleigh']:.3f} | {g['mean_rate']:+.8f} | {g['analytic_rate']:+.8f} | {g['amplification_mean']:.6g} |" for g in groups)
    md=f'''# Heating can switch a fluid disturbance from decay to growth

**Completed synthetic fluid benchmark, with full inertia.** Increasing the maintained bottom-to-top temperature difference changes the outcome in this ideal layer. At 95% of the theoretical critical heating, the initial velocity mode falls to **{below['amplification_mean']:.4f} times** its starting amplitude after six thermal diffusion times; at 105%, it reaches **{above['amplification_mean']:.3f} times** its starting amplitude. These are numerical model results, not water or pond measurements.

![Measured sign change and numerical control](figures/transition.png)

## What was completed

{s['pde_runs']} full nonlinear 2D Boussinesq runs, including controls and selected repeats; 12 separate linear modal full-inertia/Stokes comparisons; **{s['checks_passed']}/{s['checks_total']} declared final checks** and **{len(tests)-failures}/{len(tests)} independent tests** passed. Final matrix elapsed time: {s['elapsed_seconds']:.2f} seconds on the recorded machine; this excludes the rejected initial attempt, test execution and report generation. The complete 79-run matrix was not independently repeated. Selected N=32, Ra/Ra_c=1.05 trajectories were repeated exactly at two steps; the refined comparison includes all recorded arrays.

| Heating / critical heating | Rayleigh number | Refined measured growth rate | Analytic infinitesimal rate | Final / initial mode amplitude |
| --- | ---: | ---: | ---: | ---: |
{rows}

Growth is the slope of log amplitude of the critical vertical-velocity Fourier mode over t=1…6. Energy would grow at twice the amplitude rate for a pure growing mode; those rates are not interchangeable. Each refined row averages four independently seeded mixtures with matched seed across heating conditions. The four horizontal phases nearly collapse onto the same slow linear mode after the transient. Their t intervals in the raw summary are therefore extremely narrow, mainly roundoff-scale, and **are not useful bounds on numerical error or real-world uncertainty**. Step/grid changes, analytic discrepancies and modeling limitations are reported separately.

## Physical mechanism and the actual equations

Rising fluid carries temperature perturbations; buoyancy changes its vertical acceleration. With sufficient heating from below, that feedback beats diffusion and viscous loss. This is classical Rayleigh–Bénard convection. [Fabre's university lecture notes](https://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md) describe this mechanism and the stress-free neutral curve; [Ó Náraigh's course notes](https://maths.ucd.ie/~onaraigh/acm40740/acm_40740_jan2016_v1.pdf) provide further theory. These sources supply theoretical context, not experimental data analyzed here.

The model is u_t + u·grad(u) = −grad(p) + Pr Lap(u) + Pr Ra theta e_z; theta_t + u·grad(theta) = Lap(theta) + w; div(u)=0. The conductive state is motionless with temperature decreasing linearly upward. Here Ra=g alpha DeltaT H^3/(nu kappa), Pr=nu/kappa=1. Length is scaled by H, time by H^2/kappa, velocity by kappa/H, and temperature perturbation by each run's imposed DeltaT. Changing Ra represents changing DeltaT while holding geometry, coefficients and gravity fixed. No actual water properties, temperature in degrees or pond dimensions were calibrated.

The horizontal period is 2 sqrt(2) H. At z=0,H the layer is impermeable, stress-free and held at fixed temperatures: w=0, d_z u=0, theta=0. These ideal walls differ from a no-slip pond bed or a free water surface. Uniform horizontal velocity is constrained to zero; no claim of decay of an unconstrained uniform translation is made.

The same dimensionless velocity perturbation is used at every Ra, so the initial dimensional velocity matches when kappa and H are fixed. The primary thermal perturbation is the same *fraction* of DeltaT, not the same absolute temperature at different DeltaT. Two additional controls divide initial theta by Ra/Ra_c to keep its absolute size relative to the critical DeltaT fixed; they recover the same growth-rate conclusion. Halving the overall seed amplitude also preserves the linear-regime result. Four weak higher spatial modes accompany the critical roll, with independently seeded horizontal phases. Zero disturbance remains exactly zero even above the threshold in these deterministic equations.

For horizontal wavenumber k and vertical mode sin(pi z), let q=k^2+pi^2 and a=k^2/q. The independently evaluated modal matrix on (W,Theta) is [[−Pr q, Pr Ra a],[1,−q]]. Its determinant changes sign at Ra=q^3/k^2. Minimization gives k_c=pi/sqrt(2), Ra_c=27 pi^4/4={RAC:.9f}. At Pr=1, the leading eigenvalue is −q+sqrt(Ra a). The full PDE runs retain nonlinear advection and inertia; they are compared with this infinitesimal prediction, rather than generated by that formula. This study brackets the observed sign change and checks the predicted neutral point; it does not infer an exact experimental threshold from five sample settings.

## Resolution, budgets and rejected calculations

Fourier differentiation uses N=16,24,32 on the horizontally periodic domain and a reflected vertical extension of length 2H. Thus the physical half-domain has N/2 intervals in depth; N is not the number of independent vertical physical cells. Strict 2/3 filtering removes quadratic aliases. Horizontal velocity is even under vertical reflection; vertical velocity and temperature are odd. Both reflection parity and the divergence-free projection are enforced at every stage. Time advancement uses exact half diffusion steps around RK4 for advection and coupling. The split scheme is generally second order; the equal-diffusivity linear Pr=1 benchmark has commuting diffusion/coupling operators and approaches fourth order. That special result is not claimed for arbitrary nonlinear runs or Pr.

The accepted primary step is 0.005, with 0.01 and selected 0.0025 comparisons. Maximum 0.01-to-0.005 growth-rate change is {s['max_temporal_rate_change']:.3g}; maximum refined N=24-to-32 change is {s['max_spatial_rate_change']:.3g}. Separate larger-amplitude runs at theta amplitude 0.02 compare full final fields using Fourier interpolation: N=16→24 difference {s['nonlinear_field_grid_differences'][0]:.3g}; N=24→32 difference {s['nonlinear_field_grid_differences'][1]:.3g}. These smoothly resolved roll cases do not establish resolution for turbulence or different starting fields.

Kinetic energy K=mean(|u|^2)/2 obeys K'=Pr Ra mean(w theta)−Pr mean(|grad u|^2); temperature variance V=mean(theta^2)/2 obeys V'=mean(w theta)−mean(|grad theta|^2). The code integrates production with RK stages and diffusion loss through its exact substeps. It records their separate accumulated residuals, divergence and wall errors. Energy need not remain constant: the maintained thermal background supplies the growing disturbance. These budgets are for the two-dimensional nondimensional model, not total thermodynamic energy of an apparatus.

The original step 0.01 passed growth-rate checks but missed the unchanged 0.1% budget limit near neutrality: the worst relative residual was about 0.155%. Those [initial failed budget checks](data/initial-budget-checks.json) are retained. The complete finer ensemble and finer-grid controls bring the maximum refined budget residual to **{100*budget:.5g}%**. All {len(failed_coarse)} coarse runs failing the budget screen remain marked false in [cases.csv](data/cases.csv); passing final refinement does not relabel them as accepted. At step 0.2, a deliberate negative control turns the true positive growth rate into {s['coarse_control']['fitted_rate']:+.6f}: a bounded, decaying numerical trajectory can conceal a genuine physical instability.

An earlier implementation relied on analytical reflection symmetry without enforcing it each stage. Forbidden wall modes amplified roundoff and produced invalid/nonfinite results. That interrupted attempt and original source are preserved in [rejected-boundary-run](data/rejected-boundary-run/README.md), excluded from every physical conclusion. A regression test now injects a forbidden mode and verifies removal. The initial test XML also records a separate near-neutral cancellation tolerance error: the repaired comparison scales roundoff by the cancelling terms, not their nearly zero difference. These histories are not successful fluid runs.

## Does removing inertia lose real oscillations?

![Independent full-inertia and Stokes modal comparison](figures/inertia.png)

This is a separate linear modal calculation using matrix exponentials, cross-checked against an independent adaptive integrator. Removing inertia imposes W=(Ra a/q)Theta and gives Theta'=(Ra a/q−q)Theta. It preserves the linear onset threshold here, but can change the timing and erase oscillations. Matching the threshold alone does not validate a reduction.

For stable stratification Ra=−4 Ra_c and Pr=1, the full eigenvalues are −{Q:.6f} ± {2*Q:.6f} i. The thermal disturbance oscillates while its modal envelope decays. The Stokes reduction predicts monotonic decay. Initial velocity and temperature are matched using the algebraic Stokes velocity at time zero, so the difference is not a mismatched initial velocity. Its temperature-curve error over the recorded interval falls from {100*inertia[0]['theta_relative_curve_error']:.2f}% at Pr=1 to {100*inertia[3]['theta_relative_curve_error']:.3f}% at Pr=1000. This tests selected modes; it is not a proof for arbitrary initial data. [Wang (2004)](https://doi.org/10.1002/cpa.3047) treats the infinite-Prandtl limit rigorously and identifies its initial-layer character.

## What this means for the user's fluid question

**Yes, changing a physical condition can switch a small fluid disturbance from shrinking to growing.** This is a verified example in a thermally driven Boussinesq layer. The energy source is the maintained temperature difference. A footstep could seed a disturbance in a suitable physical setting, but this study has no footstep input and does not measure how much motion would reach a pond or suspended organisms.

The existing unforced Navier–Stokes study, the fitted dynamic-q equation and the pond-vibration model remain distinct. Their equations and results have not been changed. This calculation does not establish that the existing fluid runs undergo the same transition, that a positive cubic-feedback law is supported by their data, or that microbes select a branch. In a horizontally translation-invariant layer the roll phase is a continuum; opposite signed roll amplitudes can be translations of the same pattern. We therefore do not claim exactly two isolated physical states or prove a global pitchfork basin from these runs. Long-time saturation, hysteresis, nonlinear branch continuation, real no-slip/free-surface boundaries and physical calibration remain untested.

## Reproduce and inspect

Use the pinned dependencies in the parent dynamics package. Run from this directory:

```text
python -m pytest test_convection.py -q --junitxml=data/independent-tests.xml
python run_study.py
python make_report.py
python verify.py
```

Then regenerate the parent package manifest with its `make_report.py` and run its `verify_package.py`. The historical rejected artifacts are preserved, not recreated by a successful rerun. Keep copies before reproducing because current outputs are overwritten.

[Declared plan and refinement rationale](plan.json) · [All case metrics and individual pass/fail](data/cases.csv) · [Summary and environment](data/summary.json) · [Checks](data/checks.json) · [Independent test record](data/independent-tests.xml) · [Delivery audit](data/delivery-audit.json) · [Recorded column names](data/record_columns.json) · [Modal comparison metrics](data/inertia.json).

Every PDE NPZ contains the complete recorded diagnostic time series plus initial, midpoint and final fields (u,w,theta) on the reflected grid. It does **not** contain every internal integration stage. Modal NPZ files contain their time grid, full-state matrix-exponential solution and reduced solution. Sources and the plan are fingerprinted in the summary. Parent manifests cover all delivered files, including clearly rejected historical evidence.
'''
    md=md.replace('## Reproduce and inspect',next_md+'## Reproduce and inspect',1)
    (ROOT/'README.md').write_text(md,encoding='utf-8')
    table=''.join(f'<tr><td>{g["ratio"]:.2f}</td><td>{g["mean_rate"]:+.6f}</td><td>{g["amplification_mean"]:.5g}×</td></tr>' for g in groups)
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>When fluid disturbances grow</title><style>body{{margin:0;background:#f2f6f7;color:#173342;font:17px/1.6 system-ui}}main{{max-width:1100px;margin:auto;padding:32px 22px}}section{{background:white;padding:25px;margin:22px 0;border:1px solid #d1dfe5;border-radius:10px}}h1{{font-size:38px;line-height:1.2}}h2{{line-height:1.25}}img{{width:100%;height:auto}}a{{color:#076c80}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #d8e2e6}}.tag{{color:#1c6873;font-weight:bold}}.warn{{border-left:5px solid #bf6925}}small{{color:#526d79}}</style></head><body><main>
<p><a href="../report.html">Expanded dynamics study</a> · <a href="README.md">Methods, reproduction and limits</a></p><h1>Can heating turn decay into growth?</h1><p class="tag">Yes—in this bounded fluid model.</p><p>At 5% below the theoretical heating threshold, the disturbance fell to {100*below['amplification_mean']:.2f}% of its starting size. At 5% above, it grew to {above['amplification_mean']:.2f} times its starting size over the same six model time units.</p>
<section><h2>A physical change reverses the outcome</h2><img src="figures/transition.png" alt="Measured growth rates cross zero at the theoretical convection threshold; a coarse numerical step falsely predicts decay"><p>These are synthetic Boussinesq simulations of a heated layer with ideal stress-free plates. They reproduce a classical instability, rather than identify a new law or a transition in the project's existing unforced fluid runs.</p><table><tr><th>Heating / critical heating</th><th>Measured amplitude growth rate</th><th>Final / initial amplitude</th></tr>{table}</table><p>The critical Rayleigh number is {RAC:.6f} for these boundary conditions. Real ponds need different boundaries, physical properties and forcing measurements.</p></section>
<section><h2>Validation changed which runs we accepted</h2><p>{s['pde_runs']} PDE runs and 12 separate modal comparisons completed. Final checks: <strong>{s['checks_passed']}/{s['checks_total']}</strong>; independent tests: <strong>{len(tests)-failures}/{len(tests)}</strong>. Grids, time steps, seed amplitude, absolute thermal perturbation, zero disturbance and buoyancy removal were tested.</p><p>The initial step missed the energy-budget limit despite getting the growth sign right. Halving it brought the maximum refined residual to {100*budget:.5g}%. The coarser failures remain recorded. Four seeded starts agree closely because they approach the same linear mode; their tiny statistical intervals do not quantify real-world uncertainty.</p><p>An earlier wall-enforcement error produced false growth. Those runs were rejected, the boundary projection was corrected, and a regression test was added before rerunning the full matrix. <a href="data/rejected-boundary-run/README.md">Failure record</a>.</p></section>
<section><h2>Oscillation is a separate question</h2><img src="figures/inertia.png" alt="A stable mode oscillates with inertia and decays monotonically without it; reduction error falls as Prandtl number increases"><p>A stably stratified fluid mode can oscillate while decaying. Removing inertia erased those oscillations in the Pr=1 comparison. The approximation improved at large Prandtl number. Agreement about stability alone did not guarantee agreement about motion.</p></section>
<section class="warn"><h2>What this establishes—and what remains open</h2><p>Maintained heating can supply energy that amplifies a small fluid disturbance. This is a concrete example supporting the mechanism you asked about. It does not establish equivalent feedback in the existing dynamic-q fit or show that footsteps trigger this instability in a pond. No microbial response, global branch selection, long-time turbulent state or Navier–Stokes smoothness result was tested.</p><p>Read <a href="README.md">the complete study</a> for equations, controls, uncertainty, source links and raw-data definitions. <a href="https://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md">Classical convection theory</a> · <a href="https://doi.org/10.1002/cpa.3047">Inertia-free limiting theory</a>.</p></section></main></body></html>'''
    page=page.replace('<h1>Can heating', '<p><a href="#next-steps">New: next tests for the original flow and pond</a></p><h1>Can heating',1)
    page=page.replace('</main></body></html>',next_html+'</main></body></html>',1)
    (ROOT/'report.html').write_text(page,encoding='utf-8')
    print(json.dumps(dict(report=str(ROOT/'report.html'),independent_tests=len(tests),failed_tests=failures,refined_budget_max=budget)))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse-figures',action='store_true',help='Preserve existing figure bytes for a report-only update.')
    main(reuse_figures=parser.parse_args().reuse_figures)
