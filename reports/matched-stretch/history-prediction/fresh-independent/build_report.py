"""Render the frozen comparison and its limits; never relabel pilot scores as tests."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent


def main():
    p=json.loads((HERE/'protocol.json').read_text());r=json.loads((HERE/'results.json').read_text())
    primary=r['cases']['primary'];models=primary['models'];pair=primary['paired']['B_history_vs_A_full']
    ap=pair['spatial_AP'];rmse=pair['relative_RMSE_reduction']
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for ax,metric,scale,label,reference in [(axes[0],'spatial_AP',100,'Spatial AP change (percentage points)',.5),
                                          (axes[1],'relative_RMSE_reduction',100,'Relative RMSE reduction (%)',1.)]:
        labels=['B_history','A_matched_current','B_reversed','B_irrelevant','B_shuffled_91','B_mirrored_history']
        for i,name in enumerate(labels):
            v=primary['paired'][name+'_vs_A_full'][metric];lo,hi=np.array(v['ci95'])*scale;mean=v['mean']*scale
            ax.errorbar(mean,i,xerr=[[mean-lo],[hi-mean]],fmt='o',color='C0')
        ax.axvline(0,color='0.5');ax.axvline(reference,ls='--',color='C1');ax.set_yticks(range(len(labels)),labels)
        ax.set_xlabel(label);ax.set_title('Paired whole-run 95% intervals')
    fig.savefig(HERE/'comparisons.png',dpi=180);fig.savefig(HERE/'comparisons.svg');plt.close(fig)
    fig,ax=plt.subplots(figsize=(5,4),layout='constrained')
    for name in ['A_full','B_history']:
        rows=models[name]['calibration'];ax.plot([v['probability'] for v in rows],[v['observed'] for v in rows],'o-',label=name)
    ax.plot([0,1],[0,1],color='0.5',ls='--');ax.set(xlabel='Predicted event probability',ylabel='Observed event fraction',title='Held-out calibration');ax.legend()
    fig.savefig(HERE/'calibration.png',dpi=180);plt.close(fig)
    with np.load(HERE/'primary-predictions.npz',allow_pickle=False) as z:
        use=(z['group']==p['test_seeds'][0])&np.isclose(z['anchor'],.5);xyz=z['position'][use];event=z['event'][use]
        fig=plt.figure(figsize=(11,5),layout='constrained')
        for i,name in enumerate(['A_full','B_history']):
            ax=fig.add_subplot(1,2,i+1,projection='3d');prob=z[name+'_probability'][use];idx=np.argsort(-prob,kind='stable')[:int(np.ceil(.1*len(prob)))]
            ax.scatter(*xyz.T,c=event,cmap='coolwarm',s=12,vmin=0,vmax=1,alpha=.5)
            ax.scatter(*xyz[idx].T,facecolors='none',edgecolors='black',s=50)
            ax.set(xlabel='x',ylabel='y',zlabel='z',title=name+'; circles = top 10% predicted')
        fig.suptitle(f'Preselected seed{p["test_seeds"][0]}, t=0.5; red = future high deformation')
        fig.savefig(HERE/'locations.png',dpi=180);plt.close(fig)
    table=['| Model | RMSE | MAE | R² | Spatial AP | Pooled AP | Brier | ECE10 |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for name,model in models.items():
        m=model['macro'];sp=float(np.mean([v['mean_AP'] for v in model['spatial'].values()]))
        table.append(f'|{name}|{m["RMSE"]:.7f}|{m["MAE"]:.7f}|{m["R2"]:.5f}|{sp:.6f}|{model["pooled"]["AUPRC_AP"]:.6f}|{m["Brier"]:.6f}|{m["ECE10"]:.6f}|')
    report=f'''# Fresh independent marker-history test

**Primary conclusion: {r['primary_conclusion']}.** Numerical audit: {r['numerically_audited']}.

The fresh B-minus-A spatial AP gain is {ap['mean']:.7f}, paired whole-run 95% interval {ap['ci95']}. The useful threshold was fixed at .005 absolute before new simulations. Secondary mean relative RMSE reduction is {rmse['mean']:.7f}, interval {rmse['ci95']}; its useful threshold is .01. Secondary conclusion: {r['secondary_continuous_conclusion']}.

There are {len(p['test_seeds'])} independent untouched test initial conditions. All five anchors and markers from a run stay together in the bootstrap. These are distinct from the previously examined eight runs and the four numerical pilot fields. Train32 and validation16 are separate complete runs. Every selected model, threshold and scaler was frozen and published before test trajectories were generated. Validation selects regularization; no train+validation refit occurs. Event cutoff is train-only: {primary['cutoff']:.9f}. The primary measure ranks places within each anchor before averaging across anchors and independent runs. No-event anchors receive AP0 for both models and are counted in the detailed results.

{chr(10).join(table)}

![Paired comparisons](comparisons.png)
![Calibration](calibration.png)
![Preselected spatial example](locations.png)

## Exact scope and provenance

[Protocol](protocol.json), [model lock](test-lock.json), [results](results.json), [test numerical audit](test-numerical-audit.json), [execution inventory](execution.json), [independent verification](verification.json). Primary grid {r['grid']}; protocol SHA256 `{r['protocol_sha256']}`; model manifest SHA256 `{r['manifest_sha256']}`. Raw archives and per-run JSON metadata record source hashes, starts, configurations, numerical diagnostics and archive hashes. Executed source, data, selected models and predictions are published with a file inventory. [Parent provenance](../README.md).

The Eulerian Navier–Stokes equations and joint Heun arithmetic match the parent smooth random low-mode flow family: side6, zero mean/force, primary viscosity.02. The separately documented 4-point tensor Lagrange tracer sampler and its exact piecewise Jacobian replace trilinear sampling only in this new phase. Derivatives can jump at cell faces and the interpolant is not exactly divergence-free; the pilot, grid, time and accumulated-volume screens are essential. The original sharp-field periodic-join and off-central-peak limitations remain a distinct earlier experiment.

A has current velocity, all nine local gradient entries, strain and rotation summaries and current neighbor geometry (31 inputs). B has exactly those plus23 causal historical descriptors. Inputs end at the anchor, with a .025 gap before the forward .2 deformation interval. Full feature names are in the parent features.py. Predictions concern future finite-time deformation at current material-marker locations. They do not identify the exact Eulerian global maximum or first onset time. Finite-time stretch is not an asymptotic Lyapunov exponent.

## Controls and limitations

Matched-current dimensions, equal tree complexity, past-only/no-past/no-exchange ablations, independently refitted reversed/shuffled/irrelevant/mirrored histories and shuffled training labels are separately reported. Consistent marker relabeling is tested without supplying IDs as predictors. Reversal and mirror input controls do not physically reverse the fluid. A refitted reversed representation can preserve information; its performance does not establish chronology's necessity.

The ten observation sensitivities were fitted and locked before test generation; their intervals are exploratory and cannot select the primary conclusion. The viscosity variants reuse four test groups with frozen primary models and are dependent comparisons. The saved numerical failures remain part of the evidence. A failed test accuracy screen prevents usefulness or exclusion claims even if a nominal confidence interval is narrow.

Any out-of-sample benefit applies to this flow distribution, model class, observation method and time window. It does not establish material memory, causation, leadership or a Navier–Stokes breakthrough. History cannot add fundamental information beyond a complete deterministic present state; it can add predictive information to these incomplete local observations. A null or small result cannot rule out all histories, fluids or predictors.
'''
    (HERE/'README.md').write_text(report,encoding='utf8')


if __name__=='__main__':main()
