"""Rebuild charts using only the bundled measurements and analysis."""
from pathlib import Path
import sys,json,os
import numpy as np
HERE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data=json.loads((HERE/'measurements.json').read_text())
analysis=json.loads((HERE/'analysis.json').read_text())
seed=data['seed']['series'];base=data['baseline']['series']
t=np.array([r['t'] for r in seed])
fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
for rows,label,color in [(base,'128 baseline','#247c91'),(seed,'128 seed 101, 0.1%','#c76825')]:
    axes[0,0].plot(t,[r['Wmax'] for r in rows],color=color,label=label)
    axes[0,2].plot(t,[100*r['high_band_enstrophy_fraction'] for r in rows],color=color)
    axes[1,0].plot(t,[r['width_at_global_peak']['minimum_chord_cells'] for r in rows],color=color)
    axes[1,1].semilogy(t,np.maximum([abs(r['energy_budget_relative_residual']) for r in rows],1e-18),color=color)
    axes[1,2].semilogy(t,np.maximum([abs(r['enstrophy_budget_relative_residual']) for r in rows],1e-18),color=color)
axes[0,1].semilogy(t,[r['D_relative_initial_baseline_l2'] for r in seed],color='#c76825')
axes[0,2].axhline(1,color='#a23b36',ls=':',label='1% screen')
axes[1,0].axhline(6,color='#a23b36',ls=':',label='Six-cell screen')
titles=['Baseline and perturbed peak','Separation from baseline','Highest retained band','Width at global maximum','Absolute energy-budget residual','Absolute enstrophy-budget residual']
labels=['Maximum vorticity W','Difference L2 / initial baseline L2','Percent of enstrophy','Grid cells','Fraction of initial energy','Fraction of initial enstrophy']
for ax,title,label in zip(axes.flat,titles,labels):ax.set(title=title,xlabel='Model time',ylabel=label,xlim=(0,.4));ax.grid(alpha=.18)
for ax in (axes[0,0],axes[0,2],axes[1,0]):ax.legend(fontsize=8)
fig.suptitle('128-grid perturbation, seed 101: completed through 0.40',fontsize=15)
fig.savefig(HERE/'curves.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
rec=analysis['recurrence']['Wmax']
axes[0].scatter(rec['return_x'],rec['return_y'],s=13,c=t[:-1],cmap='viridis')
axes[0].set(title='Adjacent saved values',xlabel='W(t)',ylabel='W(t + 0.01)')
axes[1].plot(rec['lag_times'],rec['detrended_autocorrelation']);axes[1].set(title='Detrended autocorrelation',xlabel='Lag in model time',ylabel='Correlation')
axes[2].semilogy(rec['frequency'],np.maximum(rec['detrended_power'],1e-20));axes[2].set(title='Detrended temporal spectrum',xlabel='Cycles per model time',ylabel='Power')
for ax in axes:ax.grid(alpha=.18)
fig.suptitle('Short-window diagnostics: these curves do not establish recurrence or chaos')
fig.savefig(HERE/'recurrence.png',dpi=150);plt.close(fig)

