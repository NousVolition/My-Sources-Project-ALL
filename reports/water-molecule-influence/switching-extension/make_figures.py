"""Create scientific figures and independently audit the saved arrays."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from switching_test import equilibria, slope

HERE=Path(__file__).resolve().parent
D=HERE/'data'


def main():
    result=json.loads((D/'results.json').read_text())
    phase=json.loads((D/'phase_results.json').read_text())
    data=np.load(D/'trajectories.npz')
    pd=np.load(D/'phase_trajectories.npz')
    for arrays in [data,pd]:
        for name in arrays.files:
            assert np.isfinite(arrays[name]).all(),name
    for row in result['sweeps']:
        key=f'a{row["a"]:g}_n{row["n"]}_d{row["dwell_per_stimulus"]:g}'
        area=np.trapezoid(abs(data[key+'_down']-data[key+'_up']),data[key+'_s'])
        assert np.isclose(area,row['loop_area'],rtol=1e-12)
    for row in result['persistent']['pulse_runs']:
        assert abs(data[row['name']+'_x'][-1]-row['x_after_100_without_stimulus'])<1e-12
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    colors=['#087e83','#cc642c','#4763ab','#855390']
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    fig.suptitle('Can a temporary signal leave a persistent switch?',fontsize=18,fontweight='bold')
    x=np.linspace(0,1.9,1500);s=x-1.8*x*x/(1+x*x)
    stable=slope(x,1.8)<0
    ax[0,0].plot(s,np.where(stable,x,np.nan),color='#222',lw=2,label='Stable equilibrium')
    ax[0,0].plot(s,np.where(~stable,x,np.nan),'--',color='#777',label='Unstable equilibrium')
    key='a1.8_n2_d125'
    ax[0,0].plot(data[key+'_s'],data[key+'_up'],color=colors[0],lw=1.3,label='Increasing stimulus')
    ax[0,0].plot(data[key+'_s'],data[key+'_down'],color=colors[1],lw=1.3,label='Decreasing stimulus')
    for fold in result['reversible']['folds']:
        ax[0,0].axvline(fold['s'],color='#aaa',lw=.7,ls=':')
    ax[0,0].set(title='A | Different switching thresholds (a=1.8, n=2)',xlabel='Stimulus s',ylabel='Activity x',xlim=(0,.3),ylim=(0,1.8))
    ax[0,0].legend(frameon=False,fontsize=8)
    for name,color,label in [('short',colors[0],'Short pulse: returns low'),('long',colors[1],'Long pulse: stays high'),('no_feedback',colors[2],'No feedback: returns low')]:
        ax[0,1].plot(data[name+'_time'],data[name+'_x'],color=color,label=label)
    ax[0,1].axhline(result['persistent']['separator'],color='#777',ls=':',lw=1,label='Zero-input basin boundary')
    long=result['persistent']['pulse_runs'][1]
    ax[0,1].axvspan(5,5+long['pulse_duration'],color='#ddd',alpha=.5)
    ax[0,1].set(title='B | Pulse duration changes the final state (a=3)',xlabel='Scaled time τ',ylabel='Activity x',xlim=(0,30))
    ax[0,1].legend(frameon=False,fontsize=8)
    for a,color,label in [(1.8,colors[1],'Bistable feedback'),(1.4,colors[0],'Weak feedback'),(0.,colors[2],'No feedback')]:
        rows=[r for r in result['sweeps'] if r['a']==a]
        ax[1,0].loglog([r['dwell_per_stimulus'] for r in rows],[r['loop_area'] for r in rows],'-o',color=color,label=label)
    ax[1,0].axhline(result['reversible']['static_loop_area'],color='#777',ls=':',label='Analytic slow-sweep limit')
    ax[1,0].set(title='C | Slow sweeps separate hysteresis from lag',xlabel='Time held at each stimulus level',ylabel='Area between up/down curves')
    ax[1,0].legend(frameon=False,fontsize=8)
    for initial,color in [(0.,colors[0]),(2.,colors[1])]:
        ax[1,1].plot(data[f'basin_{initial:g}_time'],data[f'basin_{initial:g}_x'],color=color,label=f'Start x={initial:g}')
    middle=result['reversible']['basin_roots'][1]['x']
    ax[1,1].axhline(middle,color='#777',ls=':',label='Unstable separator')
    ax[1,1].set(title=f'D | Same stimulus, two stable outcomes (s={result["reversible"]["basin_stimulus"]:.3f})',xlabel='Scaled time τ',ylabel='Activity x',xlim=(0,40))
    ax[1,1].legend(frameon=False,fontsize=8)
    for a in ax.flat:a.spines[['top','right']].set_visible(False)
    fig.savefig(HERE/'switching_results.png',dpi=170)
    fig.savefig(HERE/'switching_results.svg')
    plt.close(fig)

    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    fig.suptitle('Rotation, bottlenecks, and phase locking',fontsize=18,fontweight='bold')
    delta=pd['scaling_delta']
    ax[0,0].loglog(delta,pd['scaling_period'],color=colors[0],lw=2,label='Exact period')
    ax[0,0].loglog(delta,np.pi*np.sqrt(2)/np.sqrt(delta),'--',color=colors[1],label='Near-threshold approximation')
    ax[0,0].set(title='A | Period diverges with inverse square-root scaling',xlabel='Distance to threshold δ = 1 − a (ω=1)',ylabel='Rotation period T')
    ax[0,0].legend(frameon=False,fontsize=8)
    for a,color in [(0.,colors[2]),(.9,colors[0]),(.9999,colors[1])]:
        key=f'rotation_a{a:g}'
        time=pd[key+'_time']
        ax[0,1].plot(time/time[-1],pd[key+'_theta']/(2*np.pi),color=color,label=f'a={a:g}')
    ax[0,1].set(title='B | Nearly all the cycle is spent at the bottleneck',xlabel='Fraction of one period t/T',ylabel='Unwrapped phase / 2π')
    ax[0,1].legend(frameon=False,fontsize=8)
    for mu,color in zip([.5,.99,1.01,1.1],colors):
        ax[1,0].plot(pd[f'mu{mu:g}_time'],pd[f'mu{mu:g}_phi']/(2*np.pi),color=color,label=f'μ={mu:g}: '+('locked' if mu<1 else 'slipping'))
    ax[1,0].set(title='C | Phase difference settles or repeatedly slips',xlabel='Scaled time τ',ylabel='Unwrapped phase difference / 2π',xlim=(0,60),ylim=(-.1,4.5))
    ax[1,0].legend(frameon=False,fontsize=8)
    for mu,color in [(.8,colors[0]),(1.2,colors[1])]:
        rows=[r for r in phase['pendulum'] if r['mu']==mu]
        ax[1,1].loglog([r['epsilon'] for r in rows],[r['max_unwrapped_angle_difference'] for r in rows],'-o',color=color,label=f'μ={mu:g}')
    ax[1,1].set(title='D | Full pendulum approaches the overdamped model',xlabel='Inertia parameter ε',ylabel='Maximum angle difference over τ=0–40 (rad)')
    ax[1,1].legend(frameon=False,fontsize=8)
    for a in ax.flat:a.spines[['top','right']].set_visible(False)
    fig.savefig(HERE/'phase_results.png',dpi=170)
    fig.savefig(HERE/'phase_results.svg')
    plt.close(fig)
    (D/'audit.json').write_text(json.dumps({'all_saved_arrays_finite':True,'sweep_areas_reconstructed':True,
          'pulse_endpoints_reconstructed':True,'plot_count':2,'passed':True},indent=2),encoding='utf-8')
    print('Saved-array audit and both figures completed.')


if __name__=='__main__':main()
