"""Rebuild the completed -5 degree initial-strain-axis comparison from saved records."""
from pathlib import Path
import json, os
import numpy as np
from scipy import signal
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE='baseline-n64-base'; TILT='angle_degrees-n64-v-5.0'
LABELS={BASE:'Baseline, 0 degrees',TILT:'Strain axis, -5 degrees'}

def summarize(run):
    rows=run['series'];first,last=rows[0],rows[-1];peak=max(rows,key=lambda r:r['Wmax'])
    result={'initial_W':first['Wmax'],'initial_energy':first['energy'],
        'initial_enstrophy':first['enstrophy'],'initial_tail_enstrophy_percent':100*first['high_band_enstrophy_fraction'],
        'initial_mean_velocity':first['mean_velocity'],'initial_peak_location':first['peak_location'],
        'dt':run['dt'],'sampled_peak_W':peak['Wmax'],'sampled_peak_time':peak['t'],
        'sampled_peak_over_initial':peak['Wmax']/first['Wmax'],
        'final_W':last['Wmax'],'final_I':last['I'],'final_central_ROI_W':last['central_roi_Wmax'],
        'energy_change_percent':100*(last['energy']/first['energy']-1),
        'max_energy_budget_percent':100*max(abs(r['energy_budget_relative_residual']) for r in rows),
        'max_enstrophy_budget_percent':100*max(abs(r['enstrophy_budget_relative_residual']) for r in rows),
        'final_width_cells':last['width_at_global_peak']['minimum_chord_cells'],
        'final_tail_enstrophy_percent':100*last['high_band_enstrophy_fraction'],
        'passed_resolution_rows':sum(bool(r['resolved_screen']) for r in rows),'recurrence':{}}
    for key in ['Wmax','enstrophy']:
        y=np.array([r[key] for r in rows]);d=signal.detrend(y)
        ac=signal.correlate(d,d,mode='full')[len(y)-1:];ac/=ac[0]
        f,p=signal.periodogram(y,fs=100,detrend='linear',window='hann',scaling='spectrum')
        result['recurrence'][key]={'return_x':y[:-1].tolist(),'return_y':y[1:].tolist(),
            'lag_times':[r['t'] for r in rows],'detrended_autocorrelation':ac.tolist(),
            'frequency':f.tolist(),'detrended_power':p.tolist()}
    return result

def main():
    data=json.loads((ROOT/'measurements.json').read_text());runs=data['runs'];out={k:summarize(r) for k,r in runs.items()}
    comparison={}
    for key in ['Wmax','I','energy','enstrophy','max_to_mean_spin']:
        x=np.array([r[key] for r in runs[TILT]['series']]);y=np.array([r[key] for r in runs[BASE]['series']])
        comparison[key]={'relative_curve_l2':float(np.linalg.norm(x-y)/np.linalg.norm(y)),
            'signed_endpoint_difference_percent':float(100*(x[-1]/y[-1]-1))}
    # Raw strain is S*x*exp(-0.04*r^2); rotation adds off-diagonal S_xz=S_zx.
    cross=120*np.sin(np.deg2rad(-5))*np.cos(np.deg2rad(-5))
    join={'offdiagonal_S_xz':float(cross),'max_raw_tangential_velocity_jump_at_x_or_z_face':float(6*abs(cross)*np.exp(-.04*9)),
        'formula':'Jump magnitude = 6*abs(S_xz)*exp(-0.04*(9+other_coordinate_squares)); maximum at face center. This describes the raw construction before periodic Fourier filtering/projection.'}
    result={'new_completed_runs':[TILT],'reused_completed_runs':[BASE],'saved_fields_verified':41,
        'end_time':.4,'summaries':out,'comparisons_to_baseline':comparison,'raw_join_calculation':join,
        'scope':'Completed 64-grid -5 degree tilt; +5 degree and 128-grid tilt controls are unfinished as of this batch.',
        'limitations':'Raw periodic joins are not smooth; tilt adds tangential velocity jumps; initial maxima outside central tube; no passed spatial resolution or physical peak-convergence claim.'}
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for rid in [BASE,TILT]:
        rows=runs[rid]['series'];t=[r['t'] for r in rows];w0=rows[0]['Wmax']
        ys=[[r['Wmax'] for r in rows],[r['Wmax']/w0 for r in rows],[r['I'] for r in rows],
            [r['energy']/rows[0]['energy'] for r in rows],[r['enstrophy'] for r in rows],[r['central_roi_Wmax'] for r in rows]]
        for ax,y in zip(axes.flat,ys):ax.plot(t,y,label=LABELS[rid])
    for ax,title,label in zip(axes.flat,['Maximum spin','Growth relative to own start','Accumulated maximum','Kinetic energy','Enstrophy','Central cylinder maximum'],['W','W / initial W','I','E / initial E','Z','W in radius-0.4 cylinder']):
        ax.set(title=title,xlabel='Model time',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Completed -5 degree tilt | the initial peak already differs strongly');fig.savefig(ROOT/'tilt-comparison.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for rid in [BASE,TILT]:
        rows=runs[rid]['series'];t=[r['t'] for r in rows]
        ys=[[100*r['energy_budget_relative_residual'] for r in rows],[100*r['enstrophy_budget_relative_residual'] for r in rows],
            [r['enstrophy_production'] for r in rows],[r['width_at_global_peak']['minimum_chord_cells'] for r in rows],
            [100*r['high_band_enstrophy_fraction'] for r in rows],[100*r['high_band_energy_fraction'] for r in rows]]
        for ax,y in zip(axes.flat,ys):ax.plot(t,y,label=LABELS[rid])
    for ax,title,label in zip(axes.flat,['Energy budget residual','Enstrophy budget residual','Enstrophy production','Half-peak width','Enstrophy near cutoff','Energy near cutoff'],['Initial E fraction (%)','Fraction of max(Z0, Z) (%)','Model units','Grid cells','Upper retained band (%)','Upper retained band (%)']):
        ax.set(title=title,xlabel='Model time',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=7)
    axes[1,0].axhline(6,ls=':',color='.3');axes[1,1].axhline(1,ls=':',color='.3');axes[1,2].axhline(.1,ls=':',color='.3')
    fig.suptitle('Discrete budget checks pass; the spatial-resolution screen fails');fig.savefig(ROOT/'budgets-resolution.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,key in zip(axes,['Wmax','enstrophy']):
        scale=1 if key=='Wmax' else 1e6
        for rid in [BASE,TILT]:
            r=out[rid]['recurrence'][key]
            row[0].plot(np.array(r['return_x'])/scale,np.array(r['return_y'])/scale,'.-',ms=3,label=LABELS[rid])
            row[1].plot(r['lag_times'],r['detrended_autocorrelation'],label=LABELS[rid])
            row[2].semilogy(r['frequency'][1:],np.maximum(r['detrended_power'][1:],1e-300),label=LABELS[rid])
        unit='' if scale==1 else ' / 10^6'
        row[0].set(title=key+': adjacent observations',xlabel='Value at t'+unit,ylabel='Value at t + 0.01'+unit)
        row[1].set(title=key+': autocorrelation',xlabel='Lag',ylabel='Detrended correlation')
        row[2].set(title=key+': temporal spectrum',xlabel='Frequency (1 / model time)',ylabel='Detrended power')
        for ax in row:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('41 observations through 0.40 | no demonstrated recurrence or chaos');fig.savefig(ROOT/'recurrence.png',dpi=150);plt.close(fig)
    a,b=out[BASE],out[TILT]
    text=['# Completed minus-five-degree strain-axis tilt','','**The 64-grid -5 degree run completed through time 0.40. All 41 saved fields passed numerical-record checks. It passes 0/41 spatial-resolution screens.**','',
        '[Main study](../README.md) · [Completed grid and strain comparison](../completed-grid128-strain/README.md)','','## What this tests','',
        'Only the initial background-strain axis is tilted by -5 degrees in the x-z plane. The central tube stays on its original axis. The supplied field construction, radius, spin, strain magnitude, periodic side-6 domain, viscosity 0.001, Heun integration and zero external force are preserved. Each run retains its own sampled mean and uses the prescribed starting-speed CFL rule. The completed baseline is reused; no fluid evolution was repeated.','','## A substantial starting-field change','',
        '| Quantity | Baseline | -5 degree tilt |','| --- | ---: | ---: |',
        f"| Initial W | {a['initial_W']:.6f} | {b['initial_W']:.6f} |",
        f"| Initial energy | {a['initial_energy']:.6f} | {b['initial_energy']:.6f} |",
        f"| Initial enstrophy | {a['initial_enstrophy']:.6f} | {b['initial_enstrophy']:.6f} |",
        f"| Initial upper-band enstrophy | {a['initial_tail_enstrophy_percent']:.3f}% | {b['initial_tail_enstrophy_percent']:.3f}% |",'',
        f"The initial W is already {b['initial_W']/a['initial_W']:.3f} times the baseline. A small angular change is not a small initial-vorticity perturbation for this periodically sampled construction.",'',
        'The raw strain is `u_raw = S x exp(-0.04 |x|^2)`, with `S = 40(-Id + 3 a a^T)` and `a = (sin(theta), 0, cos(theta))`. Tilting introduces off-diagonal S_xz and S_zx. Across the x and z faces these create tangential velocity jumps in the raw field before filtering/projection. Their maximum face-center magnitude is `6 |S_xz| exp(-0.36)`.',
        f"For -5 degrees, S_xz = {cross:.6f} and that raw jump is {join['max_raw_tangential_velocity_jump_at_x_or_z_face']:.6f}. At zero tilt this particular tangential jump is zero; the original raw strain still has the previously documented lack of smooth periodic joins. The filtered fields are finite, but the construction has not established a common smooth periodic limiting start.",'',
        '## Completed trajectory comparison','','| Run | Largest saved W | Peak time | Peak / initial W | Final W | Final I |','| --- | ---: | ---: | ---: | ---: | ---: |']
    for rid in [BASE,TILT]:
        r=out[rid];text.append(f"| {LABELS[rid]} | {r['sampled_peak_W']:.3f} | {r['sampled_peak_time']:.2f} | {r['sampled_peak_over_initial']:.3f} | {r['final_W']:.3f} | {r['final_I']:.3f} |")
    text+=['',f"Final I changes by {comparison['I']['signed_endpoint_difference_percent']:+.3f}% relative to the baseline. This combines the different starting field and its subsequent evolution; it cannot be attributed to a clean central-tube alignment effect.",'',
        '![Tilt curves](tilt-comparison.png)','','W is maximum grid-point vorticity magnitude; I is its recorded per-step integral. Peak times refer to observations saved every 0.01, not continuous-time maxima. The central diagnostic is a radius-0.4 cylinder spanning all z, including the periodic ends. Initial global maxima lie outside the central tube. All quantities use model units.','','## Budgets, spectrum and widths','','| Run | Energy change | Max energy residual | Max enstrophy residual | Final half-peak width | Final upper-band enstrophy |','| --- | ---: | ---: | ---: | ---: | ---: |']
    for rid in [BASE,TILT]:
        r=out[rid];text.append(f"| {LABELS[rid]} | {r['energy_change_percent']:+.4f}% | {r['max_energy_budget_percent']:.6f}% | {r['max_enstrophy_budget_percent']:.6f}% | {r['final_width_cells']:.3f} cells | {r['final_tail_enstrophy_percent']:.2f}% |")
    text+=['','![Budgets and resolution](budgets-resolution.png)','','Energy residuals use initial energy; enstrophy residuals use max(initial Z, current Z). The resolution screen requires a minimum transverse half-peak chord of six cells, upper-band energy below 0.1%, and upper-band enstrophy below 1%. Discrete budget consistency does not establish spatial resolution. No physical peak-convergence or regularity claim follows. The +5 degree counterpart and higher-grid tilt controls remain unfinished as of this batch.','','## Recurrence diagnostics','','![Recurrence](recurrence.png)','','These are adjacent-observation plots, linearly detrended autocorrelations normalized at lag zero, and Hann-window temporal spectra for W and enstrophy. Autocorrelation is not corrected for the decreasing number of pairs at long lags. Forty-one transient observations do not establish recurrence, a closed orbit or chaos, and are not asymptotic Lyapunov exponents.','','## Verification and reproduction','','All 41 new saved spectral fields were finite-checked and hashed. NumPy inverse transforms independently reproduce W, energy and enstrophy. The initial field matches its prescribed formula exactly. The preserved mean, endpoint divergence and spectral measurements, energy budget gate, and direct-versus-RHS enstrophy production were verified. Step-integrated I and budgets retain the recorded accumulations. Numerical source hashes match the approved solver and documented metadata-only revision.','','Run `python report.py` beside `measurements.json` to regenerate this report, three charts and summary. `verify.py` requires the original saved fields in the study workspace. The baseline record was previously verified.','','[Measurements](measurements.json) · [Summary](summary.json) · [Field verification](verification.json) · [Original protocol](../protocol.json)','']
    (ROOT/'README.md').write_text('\n'.join(text),encoding='utf-8')
    print(json.dumps({k:{f:v for f,v in r.items() if f!='recurrence'} for k,r in out.items()},indent=2))

if __name__=='__main__':main()
