"""Visualize one completed pilot without implying an ensemble prediction result."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,read
from solver import Flow,L

def main():
 tag='n64_nu0.025_dt0.025_seed300000';q=read('data/screen/'+tag+'.json');rows=q['rows'];c=q['config'];f=Flow(c['n'],c['nu'],c['power'])
 with np.load(ROOT/'data/screen'/(tag+'.npz')) as z:h=z['final']
 w=f.real(f.curl(h));mag=np.sqrt(np.sum(w*w,axis=0));t=[r['t'] for r in rows]
 fig,ax=plt.subplots(2,2,figsize=(10.5,7.4),layout='constrained')
 for key,a,label in [('energy',ax[0,0],'Kinetic energy'),('dissipation',ax[0,1],'Dissipation'),('re_lambda',ax[1,0],'Approximate Taylor Reynolds number')]:
  a.plot(t,[r[key] for r in rows],color='#176b87',lw=1.6);a.axvspan(20,44,color='#e4edf1',zorder=-1);a.set(xlabel='Model time',ylabel=label,title=label);a.spines[['top','right']].set_visible(False)
 ax[0,1].axhline(c['power'],color='#c87530',ls='--',label='Power input');ax[0,1].legend(fontsize=9)
 im=ax[1,1].imshow(mag[:,:,c['n']//2].T,origin='lower',extent=[-L/2,L/2,-L/2,L/2],cmap='magma',interpolation='nearest');ax[1,1].set(xlabel='x',ylabel='y',title=f'Vorticity magnitude at z = 0, time {t[-1]:g}');fig.colorbar(im,ax=ax[1,1],label='Inverse model time',shrink=.85)
 fig.suptitle('Completed 64³ pilot · independent prediction testing is a separate stage',fontsize=12)
 dest=ROOT/'results/figures/pilot.png';dest.parent.mkdir(parents=True,exist_ok=True);fig.savefig(dest,dpi=160);plt.close(fig);print(dest)
if __name__=='__main__':main()
