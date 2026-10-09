"""Render the published measurements without running a simulation."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
d=json.loads((P/'measurements.json').read_text(encoding='utf-8'))
v=json.loads((P/'verification.json').read_text(encoding='utf-8'))
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22,'figure.facecolor':'white'})
colors=['#155a9c','#ba4c22'];labels=['Ordinary step','Half step']
fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
for (rid,r),color,label in zip(d.items(),colors,labels):
 s=r['series'];t=[z['t'] for z in s]
 axs[0,0].plot(t,[z['Wmax'] for z in s],color=color,label=label)
 axs[0,1].plot(t,[z['I'] for z in s],color=color,label=label)
a=np.array([z['Wmax'] for z in d['baseline-n128-base']['series']]);b=np.array([z['Wmax'] for z in d['baseline-n128-half']['series']]);t=np.arange(41)*.01
axs[1,0].plot(t,100*np.maximum.accumulate(abs(a-b))/np.maximum.accumulate(b),color='#72549c')
axs[1,0].axhline(1,color='black',ls='--',label='1% comparison screen')
axs[1,0].axvline(.31,color='gray',ls=':',label='First failure: t=0.31')
axs[1,1].plot(t,[100*z['velocity_relative_L2_to_half'] for z in v['velocity_comparisons']],color='#267865')
for ax,title,y in zip(axs.flat,['Peak spin','Accumulated peak spin','Largest curve discrepancy so far','Whole-field timestep difference'],['W','I','Difference / half-step peak so far (%)','Velocity L2 difference / half-step L2 (%)']):
 ax.set(title=title,xlabel='Model time',ylabel=y)
axs[0,0].legend();axs[1,0].legend(fontsize=9)
fig.suptitle('128-grid timestep control: completed through 0.40',fontsize=16)
fig.savefig(P/'timestep-comparison.png',dpi=150);plt.close(fig)
fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
for (rid,r),color,label in zip(d.items(),colors,labels):
 s=r['series'];t=[z['t'] for z in s]
 axs[0,0].plot(t,[100*z['high_band_enstrophy_fraction'] for z in s],color=color,label=label)
 axs[0,1].plot(t,[z['width_at_global_peak']['minimum_chord_cells'] for z in s],color=color,label=label)
 axs[0,2].plot(t,[z['enstrophy_production'] for z in s],color=color,label=label+' production')
 axs[0,2].plot(t,[z['viscous_enstrophy_loss_rate'] for z in s],color=color,ls='--',label=label+' viscous loss')
 axs[1,0].plot(t,[z['energy_budget_relative_residual'] for z in s],color=color,label=label)
 axs[1,1].plot(t,[z['enstrophy_budget_relative_residual'] for z in s],color=color,label=label)
 w=np.array([z['Wmax'] for z in s]);axs[1,2].plot(w[:-1],w[1:],'.-',color=color,label=label,ms=3)
axs[0,0].axhline(1,color='black',ls='--');axs[0,1].axhline(6,color='black',ls='--')
titles=['Enstrophy in highest retained band','Width at the global peak','Enstrophy production and loss','Energy budget residual','Enstrophy budget residual','One-output-lag return plot']
ylabels=['Enstrophy fraction (%)','Minimum half-peak chord (cells)','Rate in model units','Relative residual','Relative residual','W(t + 0.01)']
for ax,title,y in zip(axs.flat,titles,ylabels):ax.set(title=title,xlabel='Model time',ylabel=y)
axs[1,2].set_xlabel('W(t)');axs[0,0].legend();axs[0,2].legend(fontsize=8)
fig.suptitle('Spatial screens fail; budgets and finite-window diagnostics retained',fontsize=15)
fig.savefig(P/'resolution-budgets.png',dpi=150);plt.close(fig)
