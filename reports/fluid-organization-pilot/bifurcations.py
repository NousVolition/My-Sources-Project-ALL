"""Local normal-form controls: saddle-node bottleneck and supercritical Hopf.

These finite-dimensional examples do not establish a bifurcation in the fluid.
Transcritical/pitchfork branch diagrams below are analytic illustrative forms,
not equations inferred from the cropped Figure 8.1.7.
"""
import numpy as np
from scipy.integrate import solve_ivp
from core import ROOT,save_json


def bottleneck_time(epsilon,a=1.):
    return 2/np.sqrt(epsilon)*np.arctan(a/np.sqrt(epsilon))


def hopf(t,z,mu,omega=1.):
    x,y=z;r2=x*x+y*y
    return [(mu-r2)*x-omega*y,omega*x+(mu-r2)*y]


def main():
    rows=[];arrays={}
    def exit_event(t,z):return z[0]+1
    exit_event.terminal=True;exit_event.direction=-1
    for eps in [.25,.0625,.015625,.00390625]:
        expected=bottleneck_time(eps)
        sol=solve_ivp(lambda t,z:[-eps-z[0]**2,-z[1]],(0,1.1*expected),[1.,.5],events=exit_event,dense_output=True,
                       method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('Bottleneck crossing failed')
        measured=float(sol.t_events[0][0]);t=np.linspace(0,measured,1001)
        arrays[f'time_{eps}']=t;arrays[f'path_{eps}']=sol.sol(t).T
        rows.append({'epsilon_minus_mu':eps,'analytic_passage_time':float(expected),'measured_passage_time':measured,
                     'asymptotic_pi_over_sqrt_epsilon':float(np.pi/np.sqrt(eps)),'relative_error':float(abs(measured/expected-1))})
    hrows=[];t=np.linspace(0,60,6001);arrays['hopf_time']=t
    for mu in [-.2,0.,.2]:
        for r0 in [.05,.8]:
            sol=solve_ivp(lambda t,z:hopf(t,z,mu),(0,60),[r0,0],t_eval=t,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05)
            z=sol.y.T;r=np.linalg.norm(z,axis=1);key=f'hopf_{mu}_{r0}';arrays[key]=z
            if mu==0:exact=r0/np.sqrt(1+2*r0*r0*t)
            else:exact=np.sqrt(r0*r0*np.exp(2*mu*t)/(1+(r0*r0/mu)*np.expm1(2*mu*t)))
            hrows.append({'mu':mu,'initial_radius':r0,'final_radius':float(r[-1]),'max_radial_solution_error':float(abs(r-exact).max())})
    out=ROOT/'results';np.savez_compressed(out/'bifurcations.npz',**arrays)
    save_json(out/'bifurcations.json',{'saddle_node_equations':'x_dot=mu-x^2; y_dot=-y',
              'bottleneck_crossing':'x=+1 to x=-1 at mu=-epsilon<0','bottleneck':rows,
              'hopf_equation':'z_dot=(mu+i)z-|z|^2 z','hopf':hrows,
              'scope':'Local representative normal forms. No fluid bifurcation, global attractor classification or hidden historical force is inferred.'})
    print('Saddle-node ghost passage and Hopf radial controls complete',flush=True)


def supplement():
    from analyze import RESULTS,plt,read,figsave,table,img
    r=read(RESULTS/'bifurcations.json');z=np.load(RESULTS/'bifurcations.npz')
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained');mus=np.linspace(0,.5,300)
    axs[0].plot(mus,np.sqrt(mus),label='Stable node');axs[0].plot(mus,-np.sqrt(mus),'--',label='Saddle')
    axs[0].axvline(0,color='.5',lw=.8);axs[0].set_xlim(-.15,.5);axs[0].set(title='Saddle-node: equilibria meet at μ=0',xlabel='μ',ylabel='Equilibrium x');axs[0].legend(fontsize=8)
    e=np.logspace(-3,-.4,200);axs[1].loglog(e,[bottleneck_time(a) for a in e],label='Exact finite passage')
    axs[1].loglog(e,np.pi/np.sqrt(e),'--',label='Asymptotic π/√ε')
    axs[1].scatter([x['epsilon_minus_mu'] for x in r['bottleneck']],[x['measured_passage_time'] for x in r['bottleneck']],label='ODE measurements',s=20)
    axs[1].set(title='Slow bottleneck after equilibria vanish',xlabel='ε = −μ > 0',ylabel='Time from x=1 to x=−1');axs[1].legend(fontsize=8)
    for mu,color in zip([-.2,0.,.2],['#16718a','#894cab','#c34f5a']):
        for r0 in [.05,.8]:
            q=z[f'hopf_{mu}_{r0}'];axs[2].plot(z['hopf_time'],np.linalg.norm(q,axis=1),label=f'μ={mu:g}' if r0==.05 else None,color=color,ls='-' if r0==.05 else '--')
    axs[2].axhline(np.sqrt(.2),color='.5',ls=':',lw=.8);axs[2].set(title='Supercritical Hopf: birth of a cycle',xlabel='Time (solid: r₀=0.05; dashed: r₀=0.8)',ylabel='Radius');axs[2].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'bifurcations')
    return '\n'.join([
        '<h2>7g. Added bifurcation controls: bottlenecks and emerging oscillation</h2><p>The latest excerpts distinguish the creation or loss of equilibria from the onset of an oscillation, and stress that a local portrait need not describe distant trajectories. The explicit saddle-node model supplied in the image is:</p>',
        '<div class="eq">ẋ=μ−x², &nbsp; ẏ=−y.</div>',
        '<p>For μ&gt;0, (+√μ,0) is a stable node with eigenvalues −2√μ and −1, while (−√μ,0) is a saddle. They meet at μ=0; neither equilibrium exists for μ&lt;0. Their disappearance does not remove the slow bottleneck. With μ=−ε and fixed endpoints x=±a, integration gives passage time 2 arctan(a/√ε)/√ε, approaching π/√ε as ε→0⁺. Four completed ODE runs measure this delay directly.</p>',
        img('bifurcations','The saddle-node “ghost” is a slow region of the present vector field. The Hopf control shows a different route to an organized oscillation.'),
        table(['ε=−μ','Measured passage time','Exact finite-interval time','π/√ε asymptote','Relative numerical error'],[[q['epsilon_minus_mu'],f'{q["measured_passage_time"]:.6f}',f'{q["analytic_passage_time"]:.6f}',f'{q["asymptotic_pi_over_sqrt_epsilon"]:.6f}',f'{q["relative_error"]:.2g}'] for q in r['bottleneck']]),
        '<p>The “ghost” is a description of current phase-space geometry, not a stored force or information emitted by an equilibrium that used to exist. It is relevant to this project because slow response does not, by itself, reveal weakened coupling or a history-dependent law.</p>',
        '<p>Saddle-node, transcritical and pitchfork transitions generically involve a zero eigenvalue. Representative scalar normal forms are μ−x², μx−x² and μx−x³, respectively; symmetry and nondegeneracy assumptions matter, and an arbitrary fluid cannot be assigned one from a suggestive plot. Figure 8.1.7 was supplied without its preceding system, so no specific normal form is attributed to that crop.</p>',
        '<p>A Hopf bifurcation instead involves a complex conjugate pair crossing the imaginary axis at nonzero frequency, with transversality and nonlinear nondegeneracy conditions. The separate supercritical control ż=(μ+i)z−|z|²z has ṙ=μr−r³: the origin attracts for μ&lt;0 and a stable cycle of radius √μ appears for μ&gt;0. Six numerical trajectories match the exact radial solution. At μ=0 the radius decays algebraically, illustrating why imaginary linear eigenvalues alone cannot settle nonlinear stability.</p>',
        '<p>These are established local mechanisms; see Yuri Kuznetsov’s expert accounts of the <a href="https://www.scholarpedia.org/article/Saddle-node">saddle-node bifurcation</a> and <a href="https://www.scholarpedia.org/article/Andronov-Hopf">Andronov–Hopf bifurcation</a>. To establish a corresponding fluid transition, continue an actual flow state in a physical control parameter, compute the relevant linearized fluid spectrum, track the emerging branch and verify domain/grid independence. Instantaneous local velocity-gradient eigenvalues or changing marker patterns do not supply that evidence.</p>',
        '<details><summary>Supplied bifurcation excerpts</summary><img src="references/saddle_node_reference.png" alt="User-supplied saddle-node normal form and ghost bottleneck"><img src="references/local_bifurcation_reference.png" alt="User-supplied local-validity warning and transition to Hopf bifurcation"></details>'
    ])


if __name__=='__main__':main()
