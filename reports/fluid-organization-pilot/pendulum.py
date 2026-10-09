"""Conservative nonlinear pendulum control from the newly supplied excerpt.

theta'' + sin(theta) = 0, dimensionless time sqrt(g/length)*physical time.
Energy H = omega^2/2 + 1 - cos(theta); libration amplitude alpha < pi.
SciPy ellipk takes the parameter m=sin(alpha/2)^2, not the modulus itself.
"""
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import ellipk
from core import ROOT,save_json


def rhs(t,state):
    theta,omega=state
    return [omega,-np.sin(theta)]


def hamiltonian(state):
    a=np.asarray(state)
    return .5*a[...,1]**2+1-np.cos(a[...,0])


def period(alpha):
    if not 0<=alpha<np.pi:raise ValueError('Libration amplitude must be in [0, pi)')
    return float(4*ellipk(np.sin(alpha/2)**2))


def main():
    records=[];arrays={}
    def crossing(t,y):return y[0]
    crossing.direction=-1
    for alpha in [.05,.5,1.,2.,2.8,3.05,3.13]:
        T=period(alpha);times=np.linspace(0,4*T,2001)
        sol=solve_ivp(rhs,(0,4*T),[alpha,0.],t_eval=times,events=crossing,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.1)
        values=np.diff(sol.t_events[0]);measured=float(values.mean())
        m=np.sin(alpha/2)**2
        integral=4*quad(lambda u:1/np.sqrt(1-m*np.sin(u)**2),0,np.pi/2,epsabs=1e-11)[0]
        state=sol.y.T;E=hamiltonian(state)
        records.append({'amplitude_rad':alpha,'energy':float(E[0]),'elliptic_period':T,'quadrature_period':float(integral),
                        'measured_period':measured,'relative_period_error':abs(measured/T-1),
                        'measured_cycle_period_range':float(np.ptp(values)),'max_absolute_energy_error':float(abs(E-E[0]).max())})
        arrays[f'time_{alpha}']=times;arrays[f'state_{alpha}']=state
    # A localized time impulse: same pendulum state at t=0, then one velocity kick.
    # For alpha=2.8, the critical initial speed to reach the separatrix is small.
    # Compare below and above this known threshold, using energy to classify motion.
    alpha=2.8;crit=float(np.sqrt(2*(1+np.cos(alpha))))
    kickrows=[];times=np.linspace(0,30,3001)
    for kick in [0.,.1,.3,.5]:
        sol=solve_ivp(rhs,(0,30),[alpha,kick],t_eval=times,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05)
        state=sol.y.T;E=hamiltonian(state);arrays[f'kick_{kick}']=state
        kickrows.append({'velocity_kick':kick,'initial_energy':float(E[0]),'added_energy':float(.5*kick*kick),
                         'motion':'rotation' if E[0]>2 else 'libration','max_energy_drift':float(abs(E-E[0]).max())})
    arrays['kick_time']=times
    out=ROOT/'results';out.mkdir(exist_ok=True)
    np.savez_compressed(out/'pendulum.npz',**arrays)
    save_json(out/'pendulum.json',{'equation':'theta_ddot + sin(theta) = 0','period_formula':'4 K(m), m=sin(alpha/2)^2',
              'records':records,'paired_kick_initial_angle':alpha,'rotation_threshold_speed':crit,'kick_controls':kickrows,
              'scope':'Separate conservative one-degree-of-freedom control. No damping, stochastic noise, fluid coupling or material phase transition.'})
    print('Nonlinear pendulum: periods, energy checks and paired impulse controls completed',flush=True)


def supplement():
    from analyze import RESULTS,FIG,plt,COLORS,read,figsave,table,img
    r=read(RESULTS/'pendulum.json');z=np.load(RESULTS/'pendulum.npz')
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    a=np.linspace(.001,np.pi-.003,500);axs[0].plot(a,[period(x) for x in a],label='4K(sin²(α/2))')
    axs[0].scatter([x['amplitude_rad'] for x in r['records']],[x['measured_period'] for x in r['records']],color=COLORS[1],label='Measured cycles')
    axs[0].axhline(2*np.pi,color='.4',ls=':',label='Small-amplitude 2π');axs[0].set(title='Period depends on amplitude',xlabel='Amplitude α (rad)',ylabel='Dimensionless period');axs[0].legend(fontsize=8)
    for alpha in [.5,2.,2.8,3.13]:
        q=z[f'state_{alpha}'];axs[1].plot(q[:,0],q[:,1],label=f'α={alpha}')
    theta=np.linspace(-np.pi,np.pi,400);sep=np.sqrt(2*(1+np.cos(theta)))
    axs[1].plot(theta,sep,'k:',lw=1);axs[1].plot(theta,-sep,'k:',lw=1)
    axs[1].set(title='Closed orbits approach the separatrix',xlabel='θ (rad)',ylabel='Angular velocity');axs[1].legend(fontsize=8)
    for k in [0.,.1,.3,.5]:axs[2].plot(z['kick_time'],z[f'kick_{k}'][:,0],label=f'Velocity kick {k}')
    axs[2].set(title='A controlled kick can change motion type',xlabel='Time after impulse',ylabel='Unwrapped angle θ');axs[2].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'pendulum')
    return '\n'.join([
        '<h2>7c. Added nonlinear control: the pendulum period</h2><p>The latest excerpt supplies θ̈+sinθ=0 and a period integral for libration amplitude α. Conservation of H=½θ̇²+1−cosθ gives θ̇²=2(cosθ−cosα), so a quarter-cycle integral determines the period. The substitution sin(θ/2)=sin(α/2)sinψ removes its endpoint singularity:</p>',
        '<div class="eq">T(α)=4∫₀ᵅ [2(cosθ−cosα)]⁻¹ᐟ² dθ<br>=4∫₀<sup>π/2</sup> [1−sin²(α/2)sin²ψ]⁻¹ᐟ² dψ<br>=4K(m), &nbsp; m=sin²(α/2), &nbsp; 0≤α&lt;π.</div>',
        '<p>Here K uses the elliptic parameter m. Small amplitudes recover T→2π; for a dimensional simple pendulum multiply by √(length/g). As α approaches π from below, the trajectory spends increasing time near the unstable upright state and the period diverges logarithmically. This slowing is produced by the conservative phase-space geometry, without weakened coupling or energy dissipation.</p>',
        img('pendulum','Completed nonlinear controls: numerical cycle times match the energy-integral prediction; a sufficiently large kick crosses the libration/rotation separatrix.'),
        table(['Amplitude α','Analytic period','Measured period','Relative period error','Max energy error'],[[f'{x["amplitude_rad"]:.2f}',f'{x["elliptic_period"]:.6f}',f'{x["measured_period"]:.6f}',f'{x["relative_period_error"]:.2g}',f'{x["max_absolute_energy_error"]:.2g}'] for x in r['records']]),
        f'<p>Each of seven amplitudes was followed for four predicted cycles; successive same-direction zero crossings independently measure the period. Numerical quadrature supplies a second check. Paired impulse examples start at θ=2.8, θ̇=0 and add only a velocity kick to one copy. The critical speed is {r["rotation_threshold_speed"]:.6f}; kicks 0.1 and 0.3 remain in libration, whereas 0.5 produces rotation. The extra energy is exactly ½(kick)² and is subsequently conserved to numerical accuracy. The baseline is the zero-kick trajectory.</p>',
        '<p>This connects the excerpts on centers, saddles, coupling and organization: the downward pendulum equilibrium is a center, the upright equilibrium is a saddle, and crossing its separatrix changes the qualitative motion. It does not demonstrate fluid turbulence, a thermodynamic phase transition or energy amplification. A deterministic present state (θ, θ̇) already specifies this model’s future; a history of θ can help only when angular velocity is unobserved or observations are noisy.</p>',
        '<details><summary>Supplied pendulum excerpt</summary><img src="references/pendulum_reference.png" alt="User-supplied nonlinear pendulum period exercise"></details>'
    ])


if __name__=='__main__':main()
