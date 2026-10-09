"""Render retained vorticity slices with a separate, non-overlapping colorbar."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
def main():
    data=ROOT/'data-gpu' if (ROOT/'data-gpu').exists() else ROOT/'data'
    selected=[]
    for n in [32,64,128,256]:
        folder=data/f'kida-n{n}-nu0.001-e0-s1-dt0.01'
        if (folder/'result.json').exists():selected.append((n,folder,json.loads((folder/'result.json').read_text())))
    if not selected:return
    fields=[np.load(p/'diagnostics.npz')['vorticity_slice'][-1] for n,p,r in selected]
    vmax=max(float(x.max()) for x in fields)
    fig,axes=plt.subplots(1,len(selected),figsize=(14.2,4),squeeze=False)
    for ax,(n,p,r),sl in zip(axes[0],selected,fields):
        im=ax.imshow(sl.T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],vmin=0,vmax=vmax,cmap='magma')
        status='passed' if r['rows'][-1]['local_screen_pass'] else 'failed'
        ax.set_title(f'{n}³ · local screen {status}',fontsize=9);ax.set_xlabel('x');ax.set_ylabel('y')
    fig.subplots_adjust(top=.75,bottom=.18,left=.045,right=.88,wspace=.3)
    cax=fig.add_axes([.91,.23,.012,.48])
    fig.colorbar(im,cax=cax,label='Vorticity magnitude')
    fig.suptitle('Same low-viscosity start, same endpoint, different spatial resolution\nz = π slices; common color scale; local screening is not cross-grid convergence',fontsize=11)
    fig.savefig(ROOT/'vorticity-slices.png',dpi=145);plt.close(fig)
if __name__=='__main__':main()
