"""Rebuild completed spin/strain comparisons without rerunning the fluid."""
from pathlib import Path
import os,json
import numpy as np
from scipy import signal
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
LABELS={'baseline-n64-base':'Baseline','spin_factor-n64-v0.9':'Tube spin x0.9','spin_factor-n64-v1.1':'Tube spin x1.1','strain_factor-n64-v0.9':'Background strain x0.9'}
def summarize(r):
    rows=r['series'];first,last=rows[0],rows[-1];peak=max(rows,key=lambda x:x['Wmax'])
    out={'id':r['job']['id'],'job':r['job'],'dt':r['dt'],'sampled_peak_W':peak['Wmax'],'sampled_peak_time':peak['t'],'final_W':last['Wmax'],'final_I':last['I'],'initial_energy':first['energy'],'energy_change_percent':100*(last['energy']/first['energy']-1),'energy_budget_max_percent':100*max(abs(x['energy_budget_relative_residual']) for x in rows),'enstrophy_budget_max_percent':100*max(abs(x['enstrophy_budget_relative_residual']) for x in rows),'final_width_cells':last['width_at_global_peak']['minimum_chord_cells'],'final_tail_enstrophy_percent':100*last['high_band_enstrophy_fraction'],'passed_resolution_rows':sum(bool(x['resolved_screen']) for x in rows),'recurrence':{}}
    t=np.array([v['t'] for v in rows]);assert len(t)==41 and np.allclose(np.diff(t),.01)
    for key in ('Wmax','enstrophy'):
        v=np.array([r[key] for r in rows]);d=signal.detrend(v);c=signal.correlate(d,d,mode='full')[len(v)-1:]
        if c[0]>0:c/=c[0]
        f,p=signal.periodogram(v,fs=1/(t[1]-t[0]),detrend='linear',window='hann',scaling='spectrum')
        out['recurrence'][key]={'return_x':v[:-1].tolist(),'return_y':v[1:].tolist(),'lag_times':(t-t[0]).tolist(),'detrended_autocorrelation':c.tolist(),'frequency':f.tolist(),'detrended_power':p.tolist(),'interpretation':'Short finite-window diagnostic. Drift or a spectral maximum alone does not establish recurrence or chaos.'}
    return out
def main():
    data=json.loads((ROOT/'measurements.json').read_text());runs=data['runs'];out={rid:summarize(r) for rid,r in runs.items()};baseline=runs['baseline-n64-base'];comparisons={}
    for rid in data['new_completed_runs']:
        r=runs[rid];comparisons[rid]={}
        for key in ('Wmax','I','energy','enstrophy','max_to_mean_spin'):
            a=np.array([x[key] for x in r['series']]);b=np.array([x[key] for x in baseline['series']]);comparisons[rid][key]={'relative_curve_l2':float(np.linalg.norm(a-b)/np.linalg.norm(b)),'signed_endpoint_difference_percent':float(100*(a[-1]/b[-1]-1))}
    result={'new_completed_runs':data['new_completed_runs'],'reused_completed_runs':data['reused_completed_runs'],'saved_fields_verified':82,'end_time':.4,'summaries':out,'comparisons_to_baseline':comparisons,'scope':'64-grid completed parameter runs; baseline and x0.9 tube spin reused; no new evolution'}
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    groups=[['spin_factor-n64-v0.9','baseline-n64-base','spin_factor-n64-v1.1'],['strain_factor-n64-v0.9','baseline-n64-base']]
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,ids in zip(axes,groups):
        for rid in ids:
            rows=runs[rid]['series'];t=[r['t'] for r in rows]
            for ax,key in zip(row,['Wmax','I','max_to_mean_spin']):ax.plot(t,[r[key] for r in rows],label=LABELS[rid])
        for ax,title in zip(row,['Maximum spin W','Accumulated maximum I','Maximum / mean spin']):ax.set(title=title,xlabel='Model time');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Original starting-field parameter changes | 64 grid | viscosity 0.001');fig.savefig(ROOT/'parameter-curves.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for rid in ['baseline-n64-base']+data['new_completed_runs']:
        rows=runs[rid]['series'];t=[r['t'] for r in rows];first=rows[0]
        values=[[r['energy']/first['energy'] for r in rows],[100*r['energy_budget_relative_residual'] for r in rows],[100*r['enstrophy_budget_relative_residual'] for r in rows],[r['width_at_global_peak']['minimum_chord_cells'] for r in rows],[100*r['high_band_enstrophy_fraction'] for r in rows],[r['enstrophy_production'] for r in rows]]
        for ax,y in zip(axes.flat,values):ax.plot(t,y,label=LABELS[rid])
    for ax,title,ylabel in zip(axes.flat,['Kinetic energy','Energy budget residual','Enstrophy budget residual','Half-peak width','Enstrophy near cutoff','Enstrophy production'],['E / initial E','Initial-energy fraction (%)','Fraction of max(Z0, Z) (%)','Grid cells','Upper retained band (%)','Model units']):ax.set(title=title,xlabel='Model time',ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=7)
    axes[1,0].axhline(6,color='.3',ls=':',lw=1);axes[1,1].axhline(1,color='.3',ls=':',lw=1)
    fig.suptitle('Small discrete budget residuals coexist with unresolved peaks');fig.savefig(ROOT/'budgets-resolution.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,key in zip(axes,['Wmax','enstrophy']):
        for rid in ['baseline-n64-base']+data['new_completed_runs']:
            r=out[rid]['recurrence'][key];row[0].plot(np.array(r['return_x'])/(1e6 if key=='enstrophy' else 1),np.array(r['return_y'])/(1e6 if key=='enstrophy' else 1),'.-',ms=3,label=LABELS[rid]);row[1].plot(r['lag_times'],r['detrended_autocorrelation'],label=LABELS[rid]);row[2].semilogy(r['frequency'][1:],np.maximum(r['detrended_power'][1:],1e-300),label=LABELS[rid])
        title='W' if key=='Wmax' else 'Enstrophy';units='' if key=='Wmax' else ' / 10^6';row[0].set(title=title+': adjacent observations',xlabel='Value at t'+units,ylabel='Value at t + 0.01'+units);row[1].set(title=title+': autocorrelation',xlabel='Lag',ylabel='Detrended correlation');row[2].set(title=title+': temporal spectrum',xlabel='Frequency (1 / model time)',ylabel='Detrended power')
        for ax in row:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Only 41 observations over 0.40 | these curves do not establish recurrence');fig.savefig(ROOT/'recurrence.png',dpi=150);plt.close(fig)
    text=['# Completed tube-spin and background-strain comparisons','','**Two additional 64-grid runs reached time 0.40. All 82 saved fields passed numerical-record checks; both runs fail the spatial-resolution screen at every saved time.**','','[Main study](../README.md) · [Earlier completed parameters](../completed-parameters/README.md) · [Marker roles](../marker-stress/separation-roles.md)','','## What changed','','This batch adds tube spin multiplied by 1.1 and background strain multiplied by 0.9. The completed baseline and tube-spin 0.9 records are reused. Only the named initial-field parameter changes; the supplied construction, periodic cube side 6, viscosity 0.001, per-run preserved mean, Heun method and zero force remain. The prescribed starting-speed CFL rule determines each run’s step size. No energy normalization or periodic-join repair was added. These runs change the initial field; they are not external forcing. The higher-strain and larger-grid controls remain in the existing queue.','','| Run | Largest saved W | Time of saved peak | Final W | Final I | Change in final I from baseline |','| --- | ---: | ---: | ---: | ---: | ---: |']
    for rid in data['new_completed_runs']:
        r=out[rid];text.append(f"| {LABELS[rid]} | {r['sampled_peak_W']:.3f} | {r['sampled_peak_time']:.2f} | {r['final_W']:.3f} | {r['final_I']:.3f} | {comparisons[rid]['I']['signed_endpoint_difference_percent']:+.3f}% |")
    text+=['','W is the largest grid-point vorticity magnitude. I is its recorded integral accumulated at every integration step. Peak time refers to saved observations spaced by 0.01. These are model units.','','![Parameter curves](parameter-curves.png)','','The reduced-strain run accumulates less maximum spin over this interval. Increasing the tube contribution changes the detailed maximum-spin curve while changing the final accumulated maximum much less. The tube contribution and background strain are different parts of the supplied starting field; a 10% tube change does not scale whole-box initial W by 10%. These observations are confined to the tested discrete construction.','','## Starting energy, budgets and resolution','','| Run | Initial energy | Energy change | Largest energy-budget residual | Largest enstrophy-budget residual | Final half-peak width | Final upper-band enstrophy |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for rid in ['baseline-n64-base']+data['new_completed_runs']:
        r=out[rid];text.append(f"| {LABELS[rid]} | {r['initial_energy']:.3f} | {r['energy_change_percent']:+.4f}% | {r['energy_budget_max_percent']:.6f}% | {r['enstrophy_budget_max_percent']:.6f}% | {r['final_width_cells']:.3f} cells | {r['final_tail_enstrophy_percent']:.2f}% |")
    text+=['','![Budgets and resolution](budgets-resolution.png)','','The energy residual is normalized by initial energy. The enstrophy residual is normalized by the larger of initial and current enstrophy, max(Z0, Z). Small residuals establish consistency of the discrete budgets, not adequate spatial resolution. Both new runs pass 0/41 combined resolution screens; the width requirement is six cells and the enstrophy-tail requirement is below 1%. The raw starting strain is not smooth across periodic joins, and its initial maxima are outside the central tube. Whole-box peak growth cannot be treated as verified central-tube concentration.','','## Recurrence diagnostics','','![Recurrence diagnostics](recurrence.png)','','The adjacent-observation plots, linearly detrended autocorrelations and Hann-window temporal spectra use the same diagnostics as the existing study reporter. Autocorrelation is normalized by its zero-lag value; it is not adjusted for the decreasing number of pairs at larger lags. Forty-one observations over this transient interval do not establish a closed orbit, repeated cycle, chaos or an asymptotic Lyapunov exponent. A spectral maximum or a returning autocorrelation by itself is insufficient.','','## Verification and reproduction','','Each new run is complete through 0.40 with 41 immutable fields. All were finite-checked and hashed; W, kinetic energy and enstrophy were independently remeasured using NumPy inverse transforms. Starting fields match their prescribed parameters exactly; mean drift, divergence, saved budgets and the direct-versus-RHS enstrophy-production identity were checked. Recorded step-integrated quantities were preserved. The fluid evolution was not repeated. Source hashes still match the running approved solver.','','`python report.py` rebuilds this report, all three figures and summary from `measurements.json`. The independent `verify.py` needs the original saved fields in the study workspace. The baseline and lower-spin record are reused from already verified completed batches.','','[Measurements](measurements.json) · [Summary and recurrence arrays](summary.json) · [Field verification](verification.json) · [Original protocol](../protocol.json)','']
    (ROOT/'README.md').write_text('\n'.join(text),encoding='utf-8');print(json.dumps({rid:{k:v for k,v in r.items() if k!='recurrence'} for rid,r in out.items()},indent=2))
if __name__=='__main__':main()
