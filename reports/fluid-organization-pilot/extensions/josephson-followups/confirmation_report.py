"""Append the larger independent fluid result, keeping the original pilot separate."""
import json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from support import ROOT
def add_confirmation():
    result=ROOT/'results/fluid_confirmation.json'
    if not result.exists():return
    s=json.loads(result.read_text());d=s['history_minus_present'];rel=s['relative_mean_rmse_reduction'];z=np.load(ROOT/'data/fluid_confirmation.npz')
    def ci(x):return f"{x['mean']:.6g} [{x['ci95'][0]:.6g}, {x['ci95'][1]:.6g}]"
    passed=all(s['gates'].values());resolved=s['primary_improvement_resolved'] and passed
    decision='History improved prediction in this fixed fluid regime.' if resolved else 'This larger test did not resolve a history advantage over present measurements.'
    if not passed:decision='The numerical validation failed; no confirmed prediction advantage is claimed.'
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    names=['present','history','shuffled_history'];means=np.array([s[k]['mean'] for k in names]);low=np.array([s[k]['ci95'][0] for k in names]);high=np.array([s[k]['ci95'][1] for k in names])
    axes[0].errorbar(range(3),means,yerr=[means-low,high-means],fmt='o',capsize=5,color='#156b88')
    axes[0].set(xticks=range(3),xticklabels=['Present','Present + history','Shuffled history'],ylabel='Mean per-field prediction RMSE',title='95% intervals resample whole fields');axes[0].tick_params(axis='x',labelrotation=12)
    current=np.sqrt(np.mean((z['present']-z['truth']).reshape(-1,2)**2,axis=1));history=np.sqrt(np.mean((z['history']-z['truth']).reshape(-1,2)**2,axis=1))
    lo=min(current.min(),history.min())*.8;hi=max(current.max(),history.max())*1.2
    axes[1].loglog(current,history,'o',alpha=.6,ms=4);axes[1].plot([lo,hi],[lo,hi],'--',color='gray');axes[1].set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Present-only RMSE for one field',ylabel='History-augmented RMSE for the same field',title='Below the line: history predicts better')
    fig.suptitle('120 independent fluid fields · frozen models · same forcing interventions')
    for ext in ('png','svg'):fig.savefig(ROOT/f'figures/fluid_confirmation.{ext}',bbox_inches='tight',dpi=160)
    plt.close(fig)
    budget=s['numerical'];ref=max(max(r['time_error_over_initial_energy'],r['grid_error_over_initial_energy']) for r in budget['refinements'])
    change=-100*rel['value'];change_ci=[-100*rel['ci95'][1],-100*rel['ci95'][0]]
    section=f'''<section id="fluid-confirmation"><p class="eyebrow">LARGER INDEPENDENT TEST · COMPLETED AFTER THE 12-FIELD PILOT</p><h2>{decision}</h2>
<p><strong>120 entirely new fluid starting fields</strong> (seeds 5000–5119), each cloned into unchanged, abrupt-shutoff and gradual-shutoff forcing: 360 response trajectories. Eight additional three-arm numerical controls use half the time step and a finer grid on four declared seeds. The physical parameters are the same smooth periodic-flow regime as the pilot.</p>
<p>The present-only, history and shuffled-history models were frozen before generating the new fields; their exact file hashes and the protocol were recorded. There was no refitting, re-selection, feature change or removal of poorly predicted runs. The original 18 training and six validation fields are unchanged; the old 12-field pilot is excluded from this test.</p>
<table><tr><th>Frozen predictor</th><th>Mean per-field RMSE [95% interval]</th></tr><tr><td>Present summaries</td><td>{ci(s['present'])}</td></tr><tr><td>Present + history</td><td>{ci(s['history'])}</td></tr><tr><td>Shuffled history</td><td>{ci(s['shuffled_history'])}</td></tr></table>
<p>The prespecified primary statistic, history RMSE minus present-only RMSE, is <strong>{ci(d)}</strong>. Its interval uses 20,000 paired bootstrap resamples of whole fields, preserving the two interventions together. The success rule required the upper endpoint to be below zero. Relative mean-error change with history is <strong>{change:+.2f}%</strong> [95% interval {change_ci[0]:+.2f}%, {change_ci[1]:+.2f}%]; positive values mean increased error. History predicts better in {s['fields_improved']} of 120 fields; this count is a secondary diagnostic.</p>
<figure><img src="figures/fluid_confirmation.png"><figcaption>Independent confirmation. The left panel shows group-bootstrap uncertainty; the right shows all 120 paired errors. The diagonal is equal performance.</figcaption></figure>
<p>Maximum energy-budget residual: {budget['max_energy_budget']:.3g}; burn-in residual: {budget['max_burn_budget']:.3g}; divergence: {budget['max_divergence']:.3g}; CFL: {budget['max_cfl']:.3g}; clone discrepancy: {budget['max_clone_error']:.3g}. Largest tested time/grid energy-curve difference is {ref:.3g} of initial energy. All declared numerical gates {'passed' if passed else 'did not pass'}.</p>
<p><strong>Meaning and limits:</strong> this tests whether recorded history adds useful predictive information to the specified eight current <em>global flow summaries</em> and the intervention. It does not compare against a complete velocity field, exhaust every possible present-state predictor, establish a new memory force, demonstrate turbulent or laboratory performance, or transfer Josephson physics to fluids. The result is conditional on the frozen predictor, original training sample and this distribution of flows. A different flow regime or new model selection requires another untouched validation set.</p>
<p><a href="fluid_confirmation_protocol.json">Pre-run confirmation protocol</a> · <a href="results/fluid_confirmation.csv">All 120 paired results</a> · <a href="results/fluid_confirmation.json">Statistics and numerical checks</a> · <a href="https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/fluid-organization-pilot/extensions/josephson-followups">Frozen models, hashes and recordings</a>.</p></section>'''
    path=ROOT/'report.html';text=path.read_text(encoding='utf-8')
    text=text.replace('<section><h2>1. Longer windows',section+'<section><h2>1. Longer windows',1)
    text=text.replace('Extra predictive value from recorded history remains unresolved against present measurements in these small held-out tests.','The initial pilots were inconclusive; the separate 120-field fluid confirmation below tests the frozen models on new flows.')
    text=text.replace('5. Held-out history prediction: useful hints, no resolved advantage over the present','5. Initial held-out pilots: useful hints, no resolved advantage over the present')
    path.write_text(text,encoding='utf-8')
    path=ROOT/'README.md';text=path.read_text(encoding='utf-8');heading='\n## Larger independent fluid confirmation\n\n'
    block=f'''{heading}**{decision}** The new test used **120 independent starting fields**, each cloned into three forcing arms, with the prediction models frozen beforehand.

- Present-only mean RMSE: **{s['present']['mean']:.6g}**.
- Present + history mean RMSE: **{s['history']['mean']:.6g}**.
- Paired difference: **{ci(d)}** (95% interval, 20,000 whole-field bootstrap resamples).
- Relative mean-error change with history: **{change:+.2f}%** [{change_ci[0]:+.2f}%, {change_ci[1]:+.2f}%]; positive means increased error.
- Numerical controls: **{'all passed' if passed else 'failures retained'}**; four new fields checked at half step and finer grid.

[Declared confirmation protocol](fluid_confirmation_protocol.json) · [All 120 results](results/fluid_confirmation.csv) · [Statistics and controls](results/fluid_confirmation.json).

![Independent fluid confirmation](figures/fluid_confirmation.png)

This is evidence about prediction from the chosen global flow summaries in this simulated regime. It is not a new physical memory law, a complete-state comparison, a laboratory result, or proof of transfer from the Josephson circuit.

## Earlier pilots and circuit sensitivity studies

'''
    title,sep,rest=text.partition('\n');text=title+'\n'+block+rest
    text=text.replace('The held-out prediction tests do **not** resolve','The initial small held-out prediction tests did **not** resolve')
    # Explain the extra reproduction step and never overwrite confirmation models.
    text+='\nRun `python fluid_confirmation.py --workers 4` to reproduce the larger test, or add `--analyze` to score the saved fields. Frozen models are reused and hash-checked; no new fitting is performed.\n'
    path.write_text(text,encoding='utf-8')
    summary=json.loads((ROOT/'results/summary.json').read_text());summary['fluid_confirmation']=s
    (ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':add_confirmation()
