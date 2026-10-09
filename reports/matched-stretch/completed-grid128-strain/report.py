"""Report completed 128-grid and higher-strain runs with preserved initial data."""
from pathlib import Path
import json,os
import numpy as np
from scipy import signal
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
LABELS={'baseline-n64-base':'64 grid, base step','baseline-n64-half':'64 grid, half step','baseline-n128-base':'128 grid, base step','strain_factor-n64-v0.9':'64 grid, strain x0.9','strain_factor-n64-v1.1':'64 grid, strain x1.1'}
def summary(r):
    rows=r['series'];a,b=rows[0],rows[-1];peak=max(rows,key=lambda x:x['Wmax']);t=np.array([x['t'] for x in rows])
    out={'job':r['job'],'dt':r['dt'],'initial_W':a['Wmax'],'initial_energy':a['energy'],'initial_mean_velocity':a['mean_velocity'],'sampled_peak_W':peak['Wmax'],'sampled_peak_time':peak['t'],'final_W':b['Wmax'],'final_I':b['I'],'energy_change_percent':100*(b['energy']/a['energy']-1),'energy_budget_max_percent':100*max(abs(x['energy_budget_relative_residual']) for x in rows),'enstrophy_budget_max_percent':100*max(abs(x['enstrophy_budget_relative_residual']) for x in rows),'final_width_cells':b['width_at_global_peak']['minimum_chord_cells'],'final_tail_enstrophy_percent':100*b['high_band_enstrophy_fraction'],'passed_resolution_rows':sum(bool(x['resolved_screen']) for x in rows),'final_global_peak':b['peak_location'],'final_central_ROI_W':b['central_roi_Wmax'],'recurrence':{}}
    for key in ('Wmax','enstrophy'):
        v=np.array([x[key] for x in rows]);d=signal.detrend(v);c=signal.correlate(d,d,mode='full')[len(v)-1:]
        if c[0]>0:c/=c[0]
        f,p=signal.periodogram(v,fs=1/(t[1]-t[0]),detrend='linear',window='hann',scaling='spectrum')
        out['recurrence'][key]={'return_x':v[:-1].tolist(),'return_y':v[1:].tolist(),'lag_times':(t-t[0]).tolist(),'detrended_autocorrelation':c.tolist(),'frequency':f.tolist(),'detrended_power':p.tolist(),'interpretation':'Short finite-window diagnostic. Drift or a spectral maximum alone does not establish recurrence or chaos.'}
    return out
def compare(a,b):
    assert [r['t'] for r in a['series']]==[r['t'] for r in b['series']]
    out={'through':.4,'observations':41,'reference':b['job']['id'],'fields':{}}
    for key in ('Wmax','I','energy','enstrophy','max_to_mean_spin'):
        x=np.array([r[key] for r in a['series']]);y=np.array([r[key] for r in b['series']])
        out['fields'][key]={'relative_curve_l2':float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-300)),'relative_max_error':float(np.max(abs(x-y))/max(np.max(abs(y)),1e-300)),'relative_endpoint_difference':float(abs(x[-1]-y[-1])/max(abs(y[-1]),1e-300)),'signed_endpoint_difference_percent':float(100*(x[-1]/y[-1]-1))}
    return out
def main():
    data=json.loads((ROOT/'measurements.json').read_text());runs=data['runs'];summaries={n:summary(r) for n,r in runs.items()}
    comparisons={'grid-64-128':compare(runs['baseline-n64-base'],runs['baseline-n128-base']),'timestep-64':compare(runs['baseline-n64-base'],runs['baseline-n64-half']),'strain1.1-vs-baseline':compare(runs['strain_factor-n64-v1.1'],runs['baseline-n64-base'])}
    (ROOT/'summary.json').write_text(json.dumps({'new_completed_runs':data['new_completed_runs'],'reused_completed_runs':data['reused_completed_runs'],'saved_fields_verified':82,'end_time':.4,'summaries':summaries,'comparisons':comparisons,'limitations':'Grid starts and means differ; both grids fail the peak-resolution screen. 128 half-step and 256-grid results are not complete.'},indent=2,allow_nan=False)+'\n')
    grid=['baseline-n64-base','baseline-n64-half','baseline-n128-base'];strain=['strain_factor-n64-v0.9','baseline-n64-base','strain_factor-n64-v1.1']
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for rid in grid:
        rows=runs[rid]['series'];t=[r['t'] for r in rows]
        values=[[r['Wmax'] for r in rows],[r['I'] for r in rows],[r['energy'] for r in rows],[r['width_at_global_peak']['minimum_chord_cells'] for r in rows],[100*r['high_band_enstrophy_fraction'] for r in rows],[r['central_roi_Wmax'] for r in rows]]
        for ax,y in zip(axes.flat,values):ax.plot(t,y,'--' if rid.endswith('half') else '-',label=LABELS[rid])
    for ax,title,ylabel in zip(axes.flat,['Whole-box maximum spin','Accumulated maximum spin','Kinetic energy','Half-peak width','Enstrophy near cutoff','Central cylinder maximum'],['W','I','E','Grid cells','Upper retained band (%)','W in radius-0.4 cylinder']):ax.set(title=title,xlabel='Model time',ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=7)
    axes[1,0].axhline(6,ls=':',color='.3');axes[1,1].axhline(1,ls=':',color='.3');fig.suptitle('Completed 64/128 comparison through 0.40 | grid convergence not established');fig.savefig(ROOT/'grid-comparison.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for rid in strain:
        rows=runs[rid]['series'];t=[r['t'] for r in rows];first=rows[0]
        values=[[r['Wmax'] for r in rows],[r['I'] for r in rows],[r['energy']/first['energy'] for r in rows],[r['width_at_global_peak']['minimum_chord_cells'] for r in rows],[100*r['energy_budget_relative_residual'] for r in rows],[100*r['enstrophy_budget_relative_residual'] for r in rows]]
        for ax,y in zip(axes.flat,values):ax.plot(t,y,label=LABELS[rid])
    for ax,title,ylabel in zip(axes.flat,['Maximum spin','Accumulated maximum','Kinetic energy','Half-peak width','Energy budget residual','Enstrophy budget residual'],['W','I','E / initial E','Grid cells','Fraction of initial E (%)','Fraction of max(Z0, Z) (%)']):ax.set(title=title,xlabel='Model time',ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=7)
    axes[1,0].axhline(6,ls=':',color='.3');fig.suptitle('Completed background-strain comparisons | original 64-grid construction');fig.savefig(ROOT/'strain-comparison.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,key in zip(axes,['Wmax','enstrophy']):
        for rid in ['baseline-n64-base','baseline-n128-base','strain_factor-n64-v1.1']:
            a=summaries[rid]['recurrence'][key];scale=1 if key=='Wmax' else 1e6
            row[0].plot(np.array(a['return_x'])/scale,np.array(a['return_y'])/scale,'.-',ms=3,label=LABELS[rid]);row[1].plot(a['lag_times'],a['detrended_autocorrelation'],label=LABELS[rid]);row[2].semilogy(a['frequency'][1:],np.maximum(a['detrended_power'][1:],1e-300),label=LABELS[rid])
        name='W' if key=='Wmax' else 'Enstrophy';unit='' if key=='Wmax' else ' / 10^6'
        row[0].set(title=name+': adjacent observations',xlabel='Value at t'+unit,ylabel='Value at t + 0.01'+unit);row[1].set(title=name+': autocorrelation',xlabel='Lag',ylabel='Detrended correlation');row[2].set(title=name+': temporal spectrum',xlabel='Frequency (1 / model time)',ylabel='Detrended power')
        for ax in row:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('41 saved observations per run | transient diagnostics, no recurrence proof');fig.savefig(ROOT/'recurrence.png',dpi=150);plt.close(fig)
    text=['# Completed 128-grid baseline and higher-strain run','','**The 128-grid baseline and the 64-grid strain x1.1 case both completed through time 0.40. Their 82 saved fields passed numerical-record verification. Neither passes the spatial-resolution screen.**','','[Main study](../README.md) · [Earlier spin/strain tests](../completed-spin-strain/README.md)','','## Completed runs','','| New run | Largest saved W | Time of saved peak | Final W | Final I | Final central-cylinder W |','| --- | ---: | ---: | ---: | ---: | ---: |']
    for rid in data['new_completed_runs']:
        s=summaries[rid];text.append(f"| {LABELS[rid]} | {s['sampled_peak_W']:.3f} | {s['sampled_peak_time']:.2f} | {s['final_W']:.3f} | {s['final_I']:.3f} | {s['final_central_ROI_W']:.3f} |")
    text+=['','W is maximum grid-point vorticity magnitude. I is its recorded per-step time integral. Saved observations are spaced by 0.01; peak times are sampled. All quantities use model units. The central diagnostic is the cylinder x²+y² ≤ 0.4² across all z, including its periodic ends. It is not a guarantee that an interior tube is resolved.','','## What the completed grid comparison says','','![Grid comparison](grid-comparison.png)','','| Quantity | 64/128 relative curve difference | 64-grid base/half-step relative curve difference |','| --- | ---: | ---: |']
    for key in ('Wmax','I','energy','enstrophy','max_to_mean_spin'):
        text.append(f"| {key} | {100*comparisons['grid-64-128']['fields'][key]['relative_curve_l2']:.4f}% | {100*comparisons['timestep-64']['fields'][key]['relative_curve_l2']:.6f}% |")
    text+=['','A relative curve difference is `norm(a-b)/norm(b)` over all 41 common saved times. The reference b is the 128-grid run for the grid comparison and the half-step run for the 64-grid timestep comparison. These are the same definitions used by the existing study reporter. Small 64-grid timestep differences coexist with large grid differences. The 128-grid half-step run and 256-grid controls are still unfinished, so no finest-grid convergence trend is available.','','| Grid | Initial W | Initial energy | Preserved mean velocity |','| --- | ---: | ---: | --- |']
    for rid in ['baseline-n64-base','baseline-n128-base']:
        s=summaries[rid];text.append(f"| {s['job']['n']} | {s['initial_W']:.6f} | {s['initial_energy']:.6f} | "+', '.join(f'{v:.9f}' for v in s['initial_mean_velocity'])+' |')
    text+=['','Both grids use the supplied construction, cube side 6, viscosity 0.001, Heun method, preserved per-grid mean and zero external force. Their sampled initial fields, energies and means differ. Raw strain does not match smoothly across periodic joins, and the initial global maxima lie outside the central tube. The grid difference combines initial sampling and spatial-discretization effects; it is not a clean convergence test from a demonstrated common smooth periodic field. No field replacement or mean adjustment was applied.','','## Background strain x1.1','','![Strain comparison](strain-comparison.png)','',f"The final accumulated maximum I is {comparisons['strain1.1-vs-baseline']['fields']['I']['signed_endpoint_difference_percent']:+.3f}% relative to the 64-grid baseline. The lower-strain record is reused to complete the 0.9/1.0/1.1 comparison. Changing strain changes the initial energy; it is not an externally driven run. Each uses the prescribed initial-speed CFL rule, so actual timesteps differ.",'','## Resolution and budgets','','| Run | Final half-peak width | Final upper-band enstrophy | Passed screens | Energy change | Max energy residual | Max enstrophy residual |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for rid in ['baseline-n64-base']+data['new_completed_runs']:
        s=summaries[rid];text.append(f"| {LABELS[rid]} | {s['final_width_cells']:.3f} cells | {s['final_tail_enstrophy_percent']:.2f}% | {s['passed_resolution_rows']}/41 | {s['energy_change_percent']:+.4f}% | {s['energy_budget_max_percent']:.6f}% | {s['enstrophy_budget_max_percent']:.6f}% |")
    text+=['','The declared screen requires at least six cells across the minimum transverse half-peak chord, less than 0.1% high-band energy and less than 1% high-band enstrophy. The unresolved peaks preclude a physical peak-convergence claim. Energy residuals are normalized by initial energy; enstrophy residuals use max(initial enstrophy, current enstrophy). Small budget residuals check consistency of the discrete evolution. They do not repair the starting-field or spatial-resolution limitations.','','## Recurrence diagnostics','','![Recurrence diagnostics](recurrence.png)','','Adjacent observations, linearly detrended autocorrelation and a Hann-window temporal spectrum are provided for W and enstrophy. Autocorrelation is normalized at lag zero without correcting for the number of remaining pairs. The 0.40 transient window does not establish a closed orbit, sustained cycle or chaos. These diagnostics are not asymptotic Lyapunov exponents.','','## Verification and reproduction','','All 82 new fields were finite-checked and hashed. W, energy and enstrophy were independently remeasured using NumPy inverse transforms at every saved time. Prescribed starting fields, conserved means, endpoint spectral and divergence diagnostics, discrete budgets, and direct-versus-RHS enstrophy production were checked. Integrated quantities retain the solver’s recorded accumulations. Numerical source hashes match the approved versions; any historical difference is restricted to the documented metadata-write retry. No fluid simulation was repeated.','','Rebuild the summary, three charts and this report with `python report.py` beside `measurements.json`. `verify.py` requires the original saved fields. Reused completed 64-grid and lower-strain records were already verified in prior batches.','','[Measurements](measurements.json) · [Summary and curve diagnostics](summary.json) · [Field verification](verification.json) · [Original protocol](../protocol.json)','']
    (ROOT/'README.md').write_text('\n'.join(text),encoding='utf-8');print(json.dumps({n:{k:v for k,v in s.items() if k!='recurrence'} for n,s in summaries.items()},indent=2))
if __name__=='__main__':main()
