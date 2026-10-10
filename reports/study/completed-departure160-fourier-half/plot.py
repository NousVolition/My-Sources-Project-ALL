"""Rebuild the completed Fourier 160 timestep comparison from bundled records."""
from pathlib import Path
import json
import os
import sys
import numpy as np
HERE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data=json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
proof=json.loads((HERE/'comparison-verification.json').read_text(encoding='utf-8'))
base=data['exodus-fourier-n160-base']['rows'];half=data['exodus-fourier-n160-half']['rows']
t=np.array([r['t'] for r in half])
fig,axes=plt.subplots(3,2,figsize=(13,12),layout='constrained')
for rows,style,label,color in [(base,'-','dt 0.0004','#217f97'),(half,'--','dt 0.0002','#bc552d')]:
    axes[0,0].plot(t,[r['Wmax'] for r in rows],style,color=color,label='Global, '+label)
    axes[0,0].plot(t,[r['W_central_roi'] for r in rows],style,color=color,alpha=.65,label='Central, '+label)
    axes[1,0].plot(t,[100*r['high_band_enstrophy_fraction'] for r in rows],style,color=color,label=label)
    axes[1,1].plot(t,[r['global_width']['minimum_chord_cells'] for r in rows],style,color=color,label='Global, '+label)
    axes[1,1].plot(t,[r['core_width']['minimum_chord_cells'] for r in rows],style,color=color,alpha=.65,label='Central, '+label)
    axes[2,0].semilogy(t,np.maximum([abs(r['energy_balance_relative']) for r in rows],1e-18),style,color=color,label='Energy, '+label)
    axes[2,0].semilogy(t,np.maximum([abs(r['enstrophy_balance_relative']) for r in rows],1e-18),style,color=color,alpha=.65,label='Enstrophy, '+label)
for key,label in [('Wmax','Global W'),('W_central_roi','Central W'),('I','Integral I')]:
    a=np.array([r[key] for r in base]);b=np.array([r[key] for r in half])
    difference=100*abs(a-b)/max(abs(b));difference[difference==0]=np.nan
    axes[0,1].semilogy(t,difference,label=label)
axes[1,0].axhline(1,color='#923d43',ls=':',label='1% screen')
axes[1,1].axhline(6,color='#923d43',ls=':',label='Six-cell screen')
checks=proof['field_checks'][1:]
for key,label in [('relative_velocity_L2','Velocity L2'),('relative_velocity_Linf','Velocity maximum norm')]:
    axes[2,1].semilogy([r['t'] for r in checks],[100*r[key] for r in checks],marker='o',label=label)
titles=['Global and central peaks (curves nearly overlap)','Timestep difference, normalized by each half-step peak','High-band spectral warning remains','Global and central widths remain different','Absolute budget residuals','Independent same-grid field comparison']
ylabels=['Maximum vorticity W','Difference (%)','Physical enstrophy in high band (%)','Minimum half-peak chord (cells)','Fraction of respective initial energy/enstrophy','Relative velocity difference (%)']
for ax,title,y in zip(axes.flat,titles,ylabels):
    ax.set(title=title,xlabel='Model time',ylabel=y,xlim=(0,.4));ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Departure case: Fourier grid 160, ordinary and half steps through 0.40\nClose time-step agreement; spatial-resolution screens still fail',fontsize=15)
fig.savefig(HERE/'curves.png',dpi=150);plt.close(fig)
