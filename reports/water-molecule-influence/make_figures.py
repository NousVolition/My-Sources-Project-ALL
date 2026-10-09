"""Regenerate report figures from the saved experimental data."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
D=ROOT/'data'
toy=json.loads((D/'toy_results.json').read_text())
md=json.loads((D/'md_results.json').read_text())
ta=np.load(D/'toy_arrays.npz')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
    'axes.spines.right':False,'axes.titleweight':'bold','axes.titlesize':12,'figure.dpi':150})
teal='#087e83';orange='#cf6829';blue='#4058a5';colors=[teal,orange,blue]
fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
fig.suptitle('Unequal response, changing identities',fontsize=19,fontweight='bold')
a=ax[0,0];s=ta['scores']
a.plot(np.arange(1,65),np.sort(s),color=teal,lw=2.5)
a.axhline(toy['controls']['equal_all_to_all']['mean'],color='#828894',ls='--',label='Equal-coupling control (all tied)')
a.set(xlabel='Particle rank (lowest to highest)',ylabel='Pulse share reaching other particles',
      title='A  |  Toy network: positions create unequal response',ylim=(0,.55))
a.text(.04,.88,f"Largest / smallest = {s.max()/s.min():.2f}",transform=a.transAxes)
a.legend(frameon=False,loc='lower right',fontsize=8)
a=ax[0,1];hi=toy['checks']['swap_old_leader'];lo=toy['checks']['swap_old_lowest']
x=np.arange(2)
a.bar(x-.17,[s[hi],s[lo]],.34,color=teal,label='Before position swap')
a.bar(x+.17,[s[lo],s[hi]],.34,color=orange,label='After position swap')
a.set(xticks=x,xticklabels=[f'Label {hi}',f'Label {lo}'],ylabel='Pulse share reaching other particles',
      title='B  |  Toy network: the rank follows the occupied site',ylim=(0,.65))
a.legend(frameon=False,loc='upper center',ncols=2,fontsize=8)
a=ax[1,0]
for rep in range(3):
    s=np.load(D/f'md_replica_{rep}.npz')['force_scores'][-1]
    a.plot(np.arange(1,217),np.sort(s/s.mean()),color=colors[rep],lw=2,label=f'Run {rep+1}')
a.set(xlabel='Molecule rank (lowest to highest)',ylabel='Force sensitivity / snapshot mean',
      title='C  |  TIP3P water: unequal force sensitivity at 50 ps')
a.axhline(1,color='#828894',ls=':',lw=1)
a.legend(frameon=False,fontsize=9)
a=ax[1,1]
for rep in range(3):
    leaders=md['replicas'][rep]['force_leaders']
    a.scatter(np.arange(5,51,5),[rep+1]*10,c=colors[rep],s=24)
    for t,label in zip(np.arange(5,51,5),leaders):
        a.annotate(str(label),(t,rep+1),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
a.set(xlabel='Production time (ps; snapshots every 5 ps)',ylabel='Independent run',
      title='D  |  TIP3P water: top-ranked molecule changes',
      yticks=[1,2,3],ylim=(.6,3.65),xlim=(1,54))
a.text(.02,.04,'Numbers are molecule labels, not score magnitudes.',transform=a.transAxes,fontsize=8,color='#555')
fig.savefig(ROOT/'results.png',dpi=180)
fig.savefig(ROOT/'results.svg')
plt.close(fig)

fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
fig.suptitle('Water-model diagnostics and temporal persistence',fontsize=17,fontweight='bold')
for rep in range(3):
    data=np.load(D/f'md_replica_{rep}.npz')
    times=np.arange(1,501)/10
    def smooth(v):return np.convolve(v,np.ones(20)/20,mode='valid')
    ax[0,0].plot(times[19:],smooth(data['temperature_K']),color=colors[rep],label=f'Run {rep+1}')
    ax[0,1].plot(times[19:],smooth(data['potential_kJ_mol']/216),color=colors[rep])
    lag=md['replicas'][rep]['coordination_rank_correlation_by_lag_ps']
    ax[1,1].plot([0]+[float(x) for x in lag],[1]+list(lag.values()),'-o',color=colors[rep])
ax[0,0].axhline(300,color='black',ls=':',lw=1)
ax[0,0].set(title='Temperature (2 ps moving average)',xlabel='Production time (ps)',ylabel='Temperature (K)')
ax[0,0].legend(frameon=False)
ax[0,1].set(title='Potential energy (2 ps moving average)',xlabel='Production time (ps)',ylabel='Energy per water (kJ/mol)')
rdf=np.genfromtxt(D/'oxygen_rdf.csv',delimiter=',',names=True)
ax[1,0].plot(rdf['r_nm'],rdf['g_OO'],color=teal)
ax[1,0].axhline(1,color='black',ls=':',lw=1)
ax[1,0].set(title='Oxygen–oxygen radial distribution',xlabel='Distance (nm)',ylabel='g(r)',xlim=(.2,.9))
ax[1,1].set(title='Persistence of oxygen-neighborhood rank',xlabel='Lag (ps)',ylabel='Mean same-label Spearman correlation',ylim=(-.1,1.05))
ax[1,1].axhline(0,color='black',ls=':',lw=1)
fig.savefig(ROOT/'diagnostics.png',dpi=180)
plt.close(fig)
print('Figures saved.')
