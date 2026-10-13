from pathlib import Path
import os,json
HERE=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(HERE/'mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
v=json.loads((HERE/'verification.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
for s in v['spectra']:
    k=.5*(np.array(s['k_edges'][1:])+s['k_edges'][:-1])
    for ax,key in zip(axes,['energy','enstrophy']):
        y=np.array(s[key]);keep=y>1e-25;ax.semilogy(k[keep],y[keep],label=f't={s["t"]:.2f}')
for ax,title in zip(axes,['Kinetic energy by wave-number shell','Enstrophy by wave-number shell']):
    ax.set(title=title,xlabel='Wave number (radians / domain unit)',ylabel='Shell integral');ax.grid(alpha=.2);ax.legend()
fig.suptitle('128 grid, viscosity 0.01 — saved-field spectra')
fig.savefig(HERE/'spectra.png',dpi=140);plt.close(fig)
