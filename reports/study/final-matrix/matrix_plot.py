"""Rebuild final charts from the published saved measurements; no simulation."""
from pathlib import Path
import json, os, sys
import numpy as np
HERE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
STUDY=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else HERE.parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
records={p.parent.name:read(p) for p in (STUDY/'runs').glob('*/result.json')}
comparisons=read(HERE/'comparisons.json')
cases=['aligned','compressive','exodus'];colors=['#246e8a','#b56036','#548446']
fig,axes=plt.subplots(3,3,figsize=(15,11),layout='constrained')
for i,case in enumerate(cases):
    for method,style in [('fourier','-'),('fd4','--')]:
        rows=records[f'{case}-{method}-n160-half']['rows'];t=[r['t'] for r in rows]
        axes[i,0].plot(t,[100*r['high_band_enstrophy_fraction'] for r in rows],style,label=method)
        axes[i,1].plot(t,[r['global_width']['minimum_chord_cells'] for r in rows],style,label=method+' global')
        axes[i,1].plot(t,[r['core_width']['minimum_chord_cells'] for r in rows],style,alpha=.55,label=method+' central')
        axes[i,2].semilogy(t,np.maximum([abs(r['energy_balance_relative']) for r in rows],1e-18),style,label=method+' energy')
        axes[i,2].semilogy(t,np.maximum([abs(r['enstrophy_balance_relative']) for r in rows],1e-18),style,alpha=.55,label=method+' native enstrophy')
    axes[i,0].axhline(1,ls=':',color='#9a3441',label='1% screen')
    axes[i,1].axhline(6,ls=':',color='#9a3441',label='Six-cell screen')
    for j,(title,ylabel) in enumerate([('High-band enstrophy','Fraction (%)'),('Global and central widths','Minimum half-peak chord (cells)'),('Absolute budget residuals','Fraction of respective initial total')]):
        axes[i,j].set(title=case+': '+title,xlabel='Model time',ylabel=ylabel,xlim=(0,.4));axes[i,j].grid(alpha=.2);axes[i,j].legend(fontsize=7)
fig.suptitle('All three cases: grid 160 half-step controls\nSmall budget residuals coexist with spatial-resolution warnings',fontsize=15)
fig.savefig(HERE/'finest-diagnostics.png',dpi=140);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
for case,color in zip(cases,colors):
    for method,style in [('fourier','-'),('fd4','--')]:
        for j,kind in enumerate(['grid','timestep']):
            pairs=[(48,64),(64,80),(80,112),(112,160)] if kind=='grid' else [(80,80),(112,112),(160,160)]
            for i,key in enumerate(['Wmax','I']):
                vals=[]
                for lo,hi in pairs:
                    a=f'{case}-{method}-n{lo}-base';b=f'{case}-{method}-n{hi}-'+('base' if kind=='grid' else 'half')
                    vals.append(100*comparisons[a+'__'+b]['diagnostics'][key]['max_curve_difference_relative_to_reference_peak'])
                axes[i,j].semilogy([b for a,b in pairs],vals,style,marker='o',color=color,label=case+' '+method)
for i,key in enumerate(['Global W','Integral I']):
    for j,kind in enumerate(['Adjacent grid comparison','Step halving on the same grid']):
        ax=axes[i,j];ax.set(title=key+': '+kind,xlabel='Reference grid N',ylabel='Maximum curve difference / reference peak (%)');ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Grid and timestep sensitivity through 0.40\nEach diagnostic uses its own reference-curve maximum',fontsize=15)
fig.savefig(HERE/'control-comparisons.png',dpi=140);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
for ax,case in zip(axes,cases):
    for method,style in [('fourier','-'),('fd4','--')]:
        w=np.array([r['Wmax'] for r in records[f'{case}-{method}-n160-half']['rows']])
        ax.plot(w[:-1],w[1:],style,marker='o',markersize=3,label=method)
    lim=ax.get_xlim();ax.plot(lim,lim,':',color='gray',label='Equal consecutive values')
    ax.set(title=case,xlabel='W at saved time k',ylabel='W at saved time k+1');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Consecutive scalar samples: grid 160 half steps\nTwenty pairs spaced by 0.02; these are not Poincare sections or evidence of field recurrence',fontsize=13)
fig.savefig(HERE/'successive-samples.png',dpi=140);plt.close(fig)
