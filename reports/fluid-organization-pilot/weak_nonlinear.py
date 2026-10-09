"""Weak Van der Pol attraction versus conservative hardening Duffing motion."""
import numpy as np
from scipy.integrate import solve_ivp,quad
from core import ROOT,save_json


def vdp(t,z,eps):
    x,v=z
    return [v,-x-eps*(x*x-1)*v]


def duffing(t,z,eps):
    x,v=z
    return [v,-x-eps*x**3]


def duffing_energy(z,eps):
    a=np.asarray(z);x=a[...,0];v=a[...,1]
    return .5*(v*v+x*x)+eps*x**4/4


def duffing_period(amplitude,eps):
    return 4*quad(lambda u:1/np.sqrt(1+eps*amplitude**2/2*(1+np.sin(u)**2)),0,np.pi/2,epsabs=1e-12)[0]


def downcross(t,z):return z[0]
downcross.direction=-1


def main():
    out=ROOT/'results';arrays={};eps=.1
    times=np.linspace(0,200,10001);sol=solve_ivp(lambda t,z:vdp(t,z,eps),(0,200),[.1,0.],t_eval=times,events=downcross,
                        method='DOP853',rtol=1e-10,atol=1e-12,max_step=.1)
    z=sol.y.T;radius=np.linalg.norm(z,axis=1);envelope=2/np.sqrt(1+399*np.exp(-eps*times))
    arrays.update(vdp_time=times,vdp_state=z,vdp_radius=radius,averaged_envelope=envelope)
    vv={'eps':eps,'initial_state':[.1,0.],'late_min_radius':float(radius[times>160].min()),'late_max_radius':float(radius[times>160].max()),
        'late_max_abs_x':float(abs(z[times>160,0]).max()),'late_period':float(np.diff(sol.t_events[0][-5:]).mean()),
        'averaged_envelope_formula':'r(t)=2/sqrt(1+(4/r0^2-1)*exp(-eps*t)); leading-order averaging, not exact'}
    records=[]
    for a in [.5,1.,2.]:
        T=duffing_period(a,eps);t=np.linspace(0,5*T,3001)
        s=solve_ivp(lambda t,z:duffing(t,z,eps),(0,5*T),[a,0.],t_eval=t,events=downcross,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05)
        state=s.y.T;E=duffing_energy(state,eps);measured=float(np.diff(s.t_events[0]).mean())
        records.append({'amplitude':a,'epsilon_a_squared':eps*a*a,'integral_period':float(T),'measured_period':measured,
                        'leading_frequency':1+3*eps*a*a/8,'leading_period':float(2*np.pi/(1+3*eps*a*a/8)),
                        'max_energy_drift':float(abs(E-E[0]).max())})
        arrays[f'duffing_time_{a}']=t;arrays[f'duffing_state_{a}']=state
    np.savez_compressed(out/'weak_nonlinear.npz',**arrays)
    save_json(out/'weak_nonlinear.json',{'vdp':vv,'duffing':records,'scope':'Autonomous scalar oscillator controls. Duffing epsilon=0.1 >0 is hardening, undamped and unforced. No chaos or fluid constitutive-law claim.'})
    print('Weakly nonlinear Van der Pol and Duffing controls complete',flush=True)


def supplement():
    from analyze import RESULTS,plt,read,figsave,table,img
    r=read(RESULTS/'weak_nonlinear.json');z=np.load(RESULTS/'weak_nonlinear.npz')
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    q=z['vdp_state'];axs[0].plot(q[:,0],q[:,1],lw=.6);angle=np.linspace(0,2*np.pi,400)
    axs[0].plot(2*np.cos(angle),2*np.sin(angle),'k--',lw=.8,label='Leading radius 2')
    axs[0].set(title='Weak Van der Pol: attracting cycle',xlabel='x',ylabel='ẋ',aspect='equal');axs[0].legend(fontsize=8)
    axs[1].plot(z['vdp_time'],z['vdp_radius'],lw=.7,label='Numerical phase-plane radius')
    axs[1].plot(z['vdp_time'],z['averaged_envelope'],'--',label='Leading averaged envelope')
    axs[1].set(title='Slow amplitude growth, ε=0.1',xlabel='Time',ylabel='Radius');axs[1].legend(fontsize=8)
    for a in [.5,1.,2.]:
        q=z[f'duffing_state_{a}'];axs[2].plot(q[:,0],q[:,1],label=f'Amplitude {a}')
    axs[2].set(title='Duffing: conserved closed orbits',xlabel='x',ylabel='ẋ',aspect='equal');axs[2].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'weak_nonlinear')
    return '\n'.join([
        '<h2>7e. Added weakly nonlinear controls: attraction versus conservation</h2><p>The new excerpt considers ẍ+x+εh(x,ẋ)=0 with small ε. Its two examples have different energy mechanisms, even though both can look organized in a phase portrait.</p>',
        '<div class="eq">Van der Pol: ẍ+x+ε(x²−1)ẋ=0.<br>Duffing: ẍ+x+εx³=0.</div>',
        '<p>At ε=0.1, Van der Pol started at (x,ẋ)=(0.1,0) grows slowly toward a nearly circular attracting limit cycle. Leading-order averaging predicts ṙ=(εr/2)(1−r²/4), with stable radius 2. The analytic envelope in the plot is an approximation; the actual cycle has small radius variations. This is the weak-nonlinearity regime of the same active Van der Pol oscillator used in the preceding strong fast–slow comparison.</p>',
        img('weak_nonlinear','The weak Van der Pol oscillator attracts trajectories; autonomous unforced Duffing retains distinct energy orbits.'),
        f'<p>The completed Van der Pol run has late radius range [{r["vdp"]["late_min_radius"]:.5f}, {r["vdp"]["late_max_radius"]:.5f}] and period {r["vdp"]["late_period"]:.6f}. A radius close to 2 is a leading-order prediction, not an exact circular orbit.</p>',
        '<p>For positive ε, the Duffing spring is hardening. Its energy H=½ẋ²+½x²+εx⁴/4 is conserved; there is no attracting cycle in this undamped, unforced example. The leading frequency is 1+3εa²/8, valid when εa² is small. Periods measured from repeated zero crossings agree with the full energy integral, while the approximation becomes less accurate at larger amplitudes.</p>',
        table(['Duffing amplitude','εa²','Measured period','Energy-integral period','Leading-order period','Max energy drift'],[[x['amplitude'],f'{x["epsilon_a_squared"]:.3f}',f'{x["measured_period"]:.6f}',f'{x["integral_period"]:.6f}',f'{x["leading_period"]:.6f}',f'{x["max_energy_drift"]:.2g}'] for x in r['duffing']]),
        '<p>The comparison prevents three different meanings of persistence from being conflated: an attracting oscillation, a conservative family of closed orbits, and a slowly evolving coherent fluid pattern. Autonomous conservative Duffing here is not the driven/damped Duffing system used to study chaos. None of these oscillator results alone establishes chaotic fluid behavior.</p>',
        '<details><summary>Supplied weakly nonlinear excerpt</summary><img src="references/weak_nonlinear_reference.png" alt="User-supplied weak Van der Pol and Duffing equations"></details>'
    ])


if __name__=='__main__':main()
