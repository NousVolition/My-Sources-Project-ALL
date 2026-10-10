"""Rebuild the completed-batch figures from bundled saved measurements."""
import json,os,sys
from pathlib import Path
import numpy as np
P=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(P/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
d=json.loads((P/'measurements.json').read_text(encoding='utf-8'));a=json.loads((P/'analysis.json').read_text(encoding='utf-8'))
seed=d['new']['series'];base=d['baseline']['series'];t=np.array([r['t'] for r in seed])
fig,axs=plt.subplots(2,3,figsize=(13,7.5),layout='constrained')
for rows,label in [(base,'128 baseline'),(seed,'128 seed 101, 0.5%')]:
    axs[0,0].plot(t,[r['Wmax'] for r in rows],label=label)
    axs[0,2].plot(t,[100*r['high_band_enstrophy_fraction'] for r in rows])
    axs[1,0].plot(t,[r['width_at_global_peak']['minimum_chord_cells'] for r in rows])
    axs[1,1].semilogy(t,np.maximum([abs(r['energy_budget_relative_residual']) for r in rows],1e-18))
    axs[1,2].semilogy(t,np.maximum([abs(r['enstrophy_budget_relative_residual']) for r in rows],1e-18))
axs[0,1].semilogy(t,[r['D_relative_initial_baseline_l2'] for r in seed],color='tab:orange')
axs[0,2].axhline(1,color='black',ls=':',label='1% screen');axs[1,0].axhline(6,color='black',ls=':',label='Six-cell screen')
titles=['Peak vorticity','Separation from baseline','Highest retained band','Width at global maximum','Energy budget','Enstrophy budget']
ylabels=['W','Difference L2 / initial baseline L2','Enstrophy (%)','Grid cells','Absolute relative residual','Absolute relative residual']
for ax,title,y in zip(axs.flat,titles,ylabels):ax.set(title=title,xlabel='Model time',ylabel=y,xlim=(0,.4));ax.grid(alpha=.18)
for ax in [axs[0,0],axs[0,2],axs[1,0]]:ax.legend(fontsize=8)
fig.suptitle('First completed 128-grid 0.5% perturbation: resolution remains inadequate')
fig.savefig(P/'curves.png',dpi=140);plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for key,label in [('new','128, 0.5%'),('small_amplitude','128, 0.1%'),('coarse','64, 0.5%')]:
    rows=d[key]['series'];sep=np.array([r['D_relative_initial_baseline_l2'] for r in rows])
    axs[0].semilogy(t,sep,label=label);axs[1].semilogy(t,sep/sep[0],label=label);axs[2].plot(t,[r['Wmax'] for r in rows],label=label)
for ax,title,y in zip(axs,['Separation from own baseline','Amplification from initial difference','Peak vorticity'],['D','D / D(0)','W']):ax.set(title=title,xlabel='Model time',ylabel=y);ax.grid(alpha=.18);ax.legend(fontsize=8)
fig.suptitle('Seed 101: saved amplitude and grid comparisons; no common asymptotic rate inferred')
fig.savefig(P/'amplitude-comparison.png',dpi=140);plt.close(fig)
rec=a['new']['recurrence']['Wmax'];fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
axs[0].scatter(rec['return_x'],rec['return_y'],c=t[:-1],s=14,cmap='viridis');axs[0].set(title='Adjacent saved values',xlabel='W(t)',ylabel='W(t+0.01)')
axs[1].plot(rec['lag_times'],rec['detrended_autocorrelation']);axs[1].set(title='Detrended autocorrelation',xlabel='Lag in model time',ylabel='Correlation')
axs[2].semilogy(rec['frequency'],np.maximum(rec['detrended_power'],1e-20));axs[2].set(title='Detrended temporal spectrum',xlabel='Cycles per model time',ylabel='Power')
for ax in axs:ax.grid(alpha=.18)
fig.suptitle('Short-window diagnostics: no established recurrence, periodicity or chaos')
fig.savefig(P/'recurrence.png',dpi=140);plt.close(fig)
