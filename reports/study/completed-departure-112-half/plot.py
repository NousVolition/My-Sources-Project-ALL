"""Rebuild departure-case charts from the bundled completed records."""
from pathlib import Path
import sys,json,os
import numpy as np
HERE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data=json.loads((HERE/'measurements.json').read_text())
fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
colors={48:'#8e61ad',64:'#247c91',80:'#c76825'}
for rid,r in sorted(data.items()):
    rows=r['rows'];n=r['job']['n'];method=r['job']['method'];t=[x['t'] for x in rows];style='--' if r['job']['half'] else '-';label=method+(' half step' if r['job']['half'] else ' ordinary step');color='#247c91' if method=='fourier' else '#c76825'
    for ax,key in [(axes[0,0],'Wmax'),(axes[0,1],'I')]:ax.plot(t,[x[key] for x in rows],style,color=color,label=label)
    axes[0,2].plot(t,[100*x['high_band_enstrophy_fraction'] for x in rows],style,color=color)
    axes[1,0].plot(t,[x['global_width']['minimum_chord_cells'] for x in rows],style,color=color)
    for ax,key in [(axes[1,1],'energy_balance_relative'),(axes[1,2],'enstrophy_balance_relative')]:ax.semilogy(t,np.maximum([abs(x[key]) for x in rows],1e-18),style,color=color)
axes[0,2].axhline(1,color='#a23b36',ls=':',label='1% screen');axes[1,0].axhline(6,color='#a23b36',ls=':',label='Six-cell screen')
titles=['Peak-vorticity curves','Accumulated maximum vorticity','Highest retained band','Width at global maximum','Absolute energy-budget residual','Absolute enstrophy-budget residual']
labels=['Maximum vorticity W','I','Percent of physical enstrophy','Grid cells','Fraction of initial energy','Fraction of initial native enstrophy']
for ax,title,label in zip(axes.flat,titles,labels):ax.set(title=title,xlabel='Model time',ylabel=label,xlim=(0,.4));ax.grid(alpha=.18)
for ax in (axes[0,0],axes[0,2],axes[1,0]):ax.legend(fontsize=8)
fig.suptitle('Departure case: 112-grid timestep controls through 0.40',fontsize=15)
fig.savefig(HERE/'curves.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
for method,color in [('fourier','#247c91'),('fd4','#c76825')]:
 a=data[f'exodus-{method}-n112-base']['rows'];b=data[f'exodus-{method}-n112-half']['rows'];peak=max(r['Wmax'] for r in b)
 ax.plot([r['t'] for r in a],[100*abs(x['Wmax']-y['Wmax'])/peak for x,y in zip(a,b)],label=method,color=color)
ax.set(title='Ordinary-step versus half-step W curves',xlabel='Model time',ylabel='Difference (% of half-step peak W)',xlim=(0,.4));ax.grid(alpha=.18);ax.legend()
fig.savefig(HERE/'timestep-difference.png',dpi=150);plt.close(fig)

