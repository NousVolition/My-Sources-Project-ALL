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
colors={112:'#247c91',160:'#c76825'}
for rid,r in sorted(data.items()):
    rows=r['rows'];n=r['job']['n'];method=r['job']['method'];t=[x['t'] for x in rows];style='-' if method=='fourier' else '--';label=f'{method} {n}'
    for ax,key in [(axes[0,0],'Wmax'),(axes[0,1],'I')]:ax.plot(t,[x[key] for x in rows],style,color=colors[n],label=label)
    axes[0,2].plot(t,[100*x['high_band_enstrophy_fraction'] for x in rows],style,color=colors[n])
    axes[1,0].plot(t,[x['global_width']['minimum_chord_cells'] for x in rows],style,color=colors[n])
    for ax,key in [(axes[1,1],'energy_balance_relative'),(axes[1,2],'enstrophy_balance_relative')]:ax.semilogy(t,np.maximum([abs(x[key]) for x in rows],1e-18),style,color=colors[n])
axes[0,2].axhline(1,color='#a23b36',ls=':',label='1% screen');axes[1,0].axhline(6,color='#a23b36',ls=':',label='Six-cell screen')
titles=['Peak-vorticity curves','Accumulated maximum vorticity','Highest retained band','Width at global maximum','Absolute energy-budget residual','Absolute enstrophy-budget residual']
labels=['Maximum vorticity W','I','Percent of physical enstrophy','Grid cells','Fraction of initial energy','Fraction of initial native enstrophy']
for ax,title,label in zip(axes.flat,titles,labels):ax.set(title=title,xlabel='Model time',ylabel=label,xlim=(0,.4));ax.grid(alpha=.18)
for ax in (axes[0,0],axes[0,2],axes[1,0]):ax.legend(fontsize=8)
fig.suptitle('Departure case: 112 and 160 grids through 0.40',fontsize=15)
fig.savefig(HERE/'curves.png',dpi=150);plt.close(fig)


fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j,method in enumerate(['fourier','fd4']):
    for rid,r in sorted(data.items()):
        if r['job']['method']!=method:continue
        rows=r['rows'];n=r['job']['n'];t=[z['t'] for z in rows]
        axes[j,0].plot(t,[z['Wmax'] for z in rows],color=colors[n],label=f'Global, grid {n}')
        axes[j,0].plot(t,[z['W_central_roi'] for z in rows],color=colors[n],ls='--',label=f'Central region, grid {n}')
        axes[j,1].plot(t,[z['global_width']['minimum_chord_cells'] for z in rows],color=colors[n],label=f'Global peak, grid {n}')
        axes[j,1].plot(t,[z['core_width']['minimum_chord_cells'] for z in rows],color=colors[n],ls='--',label=f'Central core, grid {n}')
    axes[j,1].axhline(6,color='gray',ls=':',label='Six-cell screen')
    for ax,title,y in zip(axes[j],[f'{method}: global and central peaks',f'{method}: separately measured widths'],['Maximum vorticity','Minimum half-peak chord (cells)']):
        ax.set(title=title,xlabel='Model time',ylabel=y,xlim=(0,.4));ax.legend(fontsize=8);ax.grid(alpha=.18)
fig.suptitle('Departure case: central cores are wider than the unresolved global peaks')
fig.savefig(HERE/'central-global.png',dpi=150);plt.close(fig)
