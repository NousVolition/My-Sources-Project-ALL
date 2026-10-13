from pathlib import Path
import os,json,hashlib,shutil
HERE=Path(__file__).resolve().parent
STUDY=Path('C:/Users/ziman/Documents/Codex/2026-10-05/ca/work/matched-stretch-study')
REPORT=HERE/'publication/reports/matched-stretch/completed-n128-viscosity0p01'
REPORT.mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(HERE/'mpl-cache')
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ids=['viscosity-n128-nu0.01','viscosity-n64-nu0.01','baseline-n128-base','baseline-n128-half','baseline-n64-base','baseline-n64-half']
results={k:read(STUDY/'runs'/k/'result.json') for k in ids}
assert all(x['status']=='complete' and len(x['series'])==41 and x['series'][-1]['t']==.4 for x in results.values())
r=results[ids[0]];rows=r['series'];t=np.array([x['t'] for x in rows]);last=rows[-1]
labels=['128, nu=0.01','64, nu=0.01','128, nu=0.001','128, nu=0.001, half dt','64, nu=0.001','64, nu=0.001, half dt']
summary={}
for rid,x in results.items():
    rr=x['series'];peak=max(rr,key=lambda a:a['Wmax'])
    summary[rid]={'grid':x['job']['n'],'viscosity':x['job']['nu'],'dt':x['dt'],'status':x['status'],'through':.4,'peak_W':peak['Wmax'],'peak_saved_t':peak['t'],'final_W':rr[-1]['Wmax'],'final_I':rr[-1]['I'],'final_energy':rr[-1]['energy'],'energy_fraction_remaining':rr[-1]['energy']/rr[0]['energy'],'initial_mean':rr[0]['mean_velocity'],'failed_spatial_screens':sum(not a['resolved_screen'] for a in rr),'max_energy_budget_residual':max(abs(a['energy_budget_relative_residual']) for a in rr),'max_enstrophy_budget_residual':max(abs(a['enstrophy_budget_relative_residual']) for a in rr),'result_sha256':sha(STUDY/'runs'/rid/'result.json')}
comparisons={}
for a,b in [(ids[0],ids[1]),(ids[0],ids[2]),(ids[2],ids[3]),(ids[4],ids[5])]:
    comparisons[a+' vs '+b]={k:float(np.linalg.norm(np.array([q[k] for q in results[a]['series']])-np.array([q[k] for q in results[b]['series']]))/np.linalg.norm([q[k] for q in results[b]['series']])) for k in ['Wmax','energy','enstrophy','I']}
rec={}
for key in ['Wmax','enstrophy']:
    v=np.array([x[key] for x in rows]);d=signal.detrend(v);ac=signal.correlate(d,d,mode='full')[len(v)-1:];ac/=ac[0]
    freq,power=signal.periodogram(v,fs=100,detrend='linear',window='hann',scaling='spectrum')
    rec[key]={'return_x':v[:-1].tolist(),'return_y':v[1:].tolist(),'lag_times':t.tolist(),'detrended_autocorrelation':ac.tolist(),'frequency':freq.tolist(),'detrended_power':power.tolist()}
save(REPORT/'measurements.json',results)
save(REPORT/'analysis.json',{'summary':summary,'relative_curve_L2_denominator_second_run':comparisons,'recurrence':rec,'recurrence_caution':'Short scalar window; no recurrence or chaos established.'})
shutil.copy2(STUDY/'protocol.json',REPORT/'protocol.json')
(REPORT/'source').mkdir(exist_ok=True)
for name in ['numerics.py','run_suite.py']:shutil.copy2(STUDY/name,REPORT/'source'/name)
shutil.copy2(HERE/'verify.py',REPORT/'verify.py')
shutil.copy2(HERE/'package.py',REPORT/'build_publication.py')
fig,axes=plt.subplots(2,3,figsize=(13,7.5),constrained_layout=True)
for ax,key,title in zip(axes.flat,['Wmax','I','energy','enstrophy','central_roi_Wmax','max_to_mean_spin'],['Maximum vorticity','Integrated maximum vorticity','Kinetic energy','Enstrophy','Central-region maximum vorticity','Maximum / mean vorticity']):
    for rid,label in zip(ids,labels):ax.plot(t,[x[key] for x in results[rid]['series']],label=label,lw=1.5)
    ax.set(title=title,xlabel='Model time');ax.grid(alpha=.2)
axes[0,0].legend(fontsize=7)
fig.suptitle('Completed controls through 0.40 — supplied matched-strain start')
fig.savefig(REPORT/'curves.png',dpi=140);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(13,7.5),constrained_layout=True)
axes[0,0].semilogy(t,np.maximum([abs(x['energy_budget_relative_residual']) for x in rows],1e-16),label='Energy')
axes[0,0].semilogy(t,np.maximum([abs(x['enstrophy_budget_relative_residual']) for x in rows],1e-16),label='Enstrophy');axes[0,0].legend();axes[0,0].set(title='Absolute budget residuals')
for key,label in [('high_band_energy_fraction','Energy'),('high_band_enstrophy_fraction','Enstrophy')]:axes[0,1].semilogy(t,[100*x[key] for x in rows],label=label)
axes[0,1].axhline(.1,ls=':',c='gray');axes[0,1].axhline(1,ls='--',c='gray');axes[0,1].legend();axes[0,1].set(title='Outer retained modes (%)')
for key,label in [('width_at_global_peak','Global peak'),('width_at_central_roi_peak','Central peak')]:axes[0,2].plot(t,[x[key]['minimum_chord_cells'] for x in rows],label=label)
axes[0,2].axhline(6,ls=':',c='gray');axes[0,2].legend();axes[0,2].set(title='Minimum transverse chord (cells)')
for key,label in [('stretch_parallel_at_Wmax','Stretch'),('viscosity_parallel_at_Wmax','Viscosity')]:axes[1,0].plot(t,[x[key] for x in rows],label=label)
axes[1,0].legend();axes[1,0].set(title='Signed local vorticity terms')
for key,label in [('enstrophy_production','Production'),('viscous_enstrophy_loss_rate','Viscous loss')]:axes[1,1].plot(t,[x[key] for x in rows],label=label)
axes[1,1].legend();axes[1,1].set(title='Enstrophy rates')
axes[1,2].plot(t,[x['half_peak_volume'] for x in rows]);axes[1,2].set(title='Half-peak volume')
for ax in axes.flat:ax.set(xlabel='Model time');ax.grid(alpha=.2)
fig.suptitle('128 grid, viscosity 0.01 — numerical budget and resolution checks')
fig.savefig(REPORT/'budgets-resolution.png',dpi=140);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
for row,key in enumerate(rec):
    x=rec[key];axes[row,0].plot(x['return_x'],x['return_y'],'.-');axes[row,0].set(xlabel=f'{key}(t)',ylabel=f'{key}(t+0.01)',title='Scalar return plot')
    axes[row,1].plot(x['lag_times'],x['detrended_autocorrelation']);axes[row,1].set(xlabel='Lag',title='Detrended autocorrelation')
    axes[row,2].semilogy(x['frequency'][1:],x['detrended_power'][1:]);axes[row,2].set(xlabel='Frequency (cycles / model time)',title='Hann-window periodogram')
fig.suptitle('Short-window diagnostics — no recurrence or chaos established')
fig.savefig(REPORT/'recurrence.png',dpi=140);plt.close(fig)
s=summary[ids[0]]
table='\n'.join(f"| {label} | {summary[rid]['peak_W']:.6g} | {summary[rid]['peak_saved_t']:.2f} | {summary[rid]['final_W']:.6g} | {summary[rid]['final_I']:.6g} | {summary[rid]['failed_spatial_screens']}/41 |" for rid,label in zip(ids,labels))
text=f'''# Completed 128-grid viscosity 0.01 control

[Study index](../README.md) · [All reports](../../README.md)

The existing unperturbed viscosity-0.01 run completed all **41 saved outputs through model time 0.40**. Its sampled peak vorticity was **{s['peak_W']:.6g} at t={s['peak_saved_t']:.2f}**; final vorticity was **{s['final_W']:.6g}**, and the accumulated maximum-vorticity integral was **{s['final_I']:.6g}**. Final kinetic energy was **{s['energy_fraction_remaining']:.3%}** of its initial value.

This is an unperturbed viscosity sweep. It has no disturbance-separation diagnostic. Peak-vorticity growth cannot answer whether a disturbance shrinks. The separately authorized disturbed/reference comparisons remain in progress.

## Scope and numerical controls

The supplied initial field, periodic cube of side 6, fixed discrete mean, zero force and frozen Heun integrator are retained. Grid: 128 cubed. Viscosity: 0.01. Step: {r['dt']:.15g}, with 205 substeps per 0.01 saved interval. The initial mean is {s['initial_mean']}. Per-grid discrete means differ in the sampled nonsmooth initial field and are preserved without renormalizing controls.

All comparisons below are existing completed runs at the identical 41 saved times. Each retains its own viscosity, grid, initial norm and step. The 256-grid run is unfinished and is excluded from these completed comparisons. There is no viscosity-0.01 half-step pair, and no time-convergence claim for this new result.

| Completed run | Sampled peak W | Peak time | Final W | Final I | Failed spatial screens |
| --- | ---: | ---: | ---: | ---: | ---: |
{table}

![Completed curve comparisons](curves.png)

The 64/128 viscosity-0.01 W-curve relative L2 difference is {comparisons[ids[0]+' vs '+ids[1]]['Wmax']:.3%}, with the 64-grid curve as denominator. Changing viscosity changes the model being integrated; a smaller peak or smaller budget residual does not demonstrate a physical stability transition.

## Budgets, spectra and widths

The maximum absolute relative energy-budget residual is {s['max_energy_budget_residual']:.6g}; the enstrophy-budget residual is {s['max_enstrophy_budget_residual']:.6g}. The declared energy hard stop is 0.005; the maximum recorded speed CFL is {r['max_speed_cfl']:.6g}, below 0.5. No declared hard stop occurred. Saved stage integrals are preserved; this verification does not repeat integration.

The spatial screen fails at **{s['failed_spatial_screens']}/41** saved outputs. At 0.40, the outer retained modes contain **{100*last['high_band_energy_fraction']:.4g}% of energy** and **{100*last['high_band_enstrophy_fraction']:.4g}% of enstrophy**. The minimum global-peak chord spans **{last['width_at_global_peak']['minimum_chord_cells']:.4g} cells**; the central-peak chord spans **{last['width_at_central_roi_peak']['minimum_chord_cells']:.4g} cells**. The screen requires less than 0.1% outer-band energy, less than 1% outer-band enstrophy and at least 6 cells across the global-peak chord. It is an exploratory screen, not proof of resolution.

![Budget and resolution checks](budgets-resolution.png)

![Saved velocity and vorticity spectra](spectra.png)

**The raw strain is nonsmooth at periodic joins, and the initial maximum lies outside the central tube. Spatial and time convergence are unestablished.** Small budget errors check consistency of the discrete calculation; they do not repair unresolved structure or establish continuum behavior.

## Recurrence diagnostics

Return plots, linearly detrended autocorrelations and Hann-window periodograms use the saved W and enstrophy series. Forty-one samples over 0.40 cannot establish recurrence, chaos or asymptotic Lyapunov behavior. No perturbation rate is calculated for this unperturbed run.

![Finite-window recurrence diagnostics](recurrence.png)

## Verification and records

All 41 fields were independently remeasured using NumPy physical-space transforms for vorticity, energy, enstrophy, mean, divergence, stretching production and viscous losses. Spectral fractions, global and central widths and budget residuals were rechecked from saved fields with the frozen observation code. The prescribed initial field and original source fingerprints were verified. Controls in the comparison table reuse their completed saved diagnostic records; no new control trajectory was computed.

- [Raw diagnostics for the six completed comparisons](measurements.json)
- [Curve differences, summaries and recurrence arrays](analysis.json)
- [Independent verification, field hashes and spectra](verification.json)
- [Verifier](verify.py), [publication builder](build_publication.py), [spectrum plotter](plot_spectra.py)
- [Frozen protocol](protocol.json), [numerical source](source/numerics.py), [controller source](source/run_suite.py)

Run the verifier with the original study directory as its argument; the large saved fields remain locally preserved. The publication builder expects this workspace layout and records provenance, rather than serving as a standalone simulator. The frozen protocol retains its historical publication language; later explicit authorization covers this completed numerical batch without changing numerical settings.

At publication preparation, the original matched-stretch matrix has **35/48 completed jobs, 2 active jobs and 11 queued jobs**. The remaining matrix and both disturbance comparisons continue. The separate vortex matrix is already complete and published. Historical intake reviews are unchanged.
'''
(REPORT/'README.md').write_text(text,encoding='utf-8')
print(json.dumps({'summary':s,'grid_W_curve_difference':comparisons[ids[0]+' vs '+ids[1]]['Wmax'],'report':str(REPORT)},indent=2))
