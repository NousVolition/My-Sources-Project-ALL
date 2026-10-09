"""Representative fast-slow Van der Pol control for the supplied cubic nullcline.

F(x)=x^3/3-x, x'=mu*(y-F(x)), y'=-x/mu, mu>0.
This is equivalent to x''-mu*(1-x^2)*x'+x=0. It is an active oscillator,
not a passive isolated fluid: negative damping supplies energy at small |x|.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from core import ROOT,save_json


def cubic(x):return x**3/3-x


def rhs(t,state,mu):
    x,y=state[:2];v=mu*(y-cubic(x))
    ans=[v,-x/mu]
    if len(state)==3:ans.append(mu*(1-x*x)*v*v)
    return ans


def energy(xy,mu):
    x,y=np.asarray(xy)[...,:2].T
    return .5*(x*x+(mu*(y-cubic(x)))**2)


def upcross(t,z):return z[0]
upcross.direction=1


def solve(mu,z0,end,dense=True):
    def jac(t,state):
        x,y=state[:2];v=mu*(y-cubic(x));J=np.zeros((len(state),len(state)))
        J[0,:2]=[mu*(1-x*x),mu];J[1,0]=-1/mu
        if len(state)==3:
            J[2,0]=mu*(-2*x*v*v+2*mu*(1-x*x)**2*v)
            J[2,1]=2*mu*mu*(1-x*x)*v
        return J
    sol=solve_ivp(lambda t,z:rhs(t,z,mu),(0,end),z0,method='Radau',rtol=1e-9,atol=1e-11,
                  dense_output=dense,events=upcross,max_step=.25,jac=jac)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def main():
    out=ROOT/'results';arrays={};rows=[]
    for mu in [1.,5.,10.,30.]:
        approx=mu*(3-2*np.log(2));end=10*(2*np.pi+approx)
        sol=solve(mu,[2.,0.,0.],end);cross=sol.t_events[0];periods=np.diff(cross[-5:]);T=float(periods.mean())
        ts=np.linspace(cross[-2],cross[-1],6001);xy=sol.sol(ts)[:2].T
        all_t=np.linspace(0,end,16001);full=sol.sol(all_t).T;E=energy(full[:,:2],mu)
        residual=E-E[0]-full[:,2]
        # Slow-branch departures: the cubic has folds at x=+/-1.
        # Fraction of time in the vicinity of the attracting outer branches.
        slow=float(np.mean((abs(xy[:,0])>1)&(abs(xy[:,1]-cubic(xy[:,0]))<.1)))
        arrays[f'cycle_time_{mu}']=ts-ts[0];arrays[f'cycle_{mu}']=xy
        rows.append({'mu':mu,'measured_period':T,'period_range_last_four':float(np.ptp(periods)),
                     'large_mu_period':float(approx),'relative_asymptotic_period_error':float(abs(approx/T-1)),
                     'slow_branch_time_fraction':slow,'max_energy_balance_absolute_residual':float(abs(residual).max())})
        if mu==10:
            anchor=sol.sol(cross[-2]+T/4)[:2];kick=anchor+np.array([.2,0.]);t=np.linspace(0,6*T,6001)
            b=solve(mu,anchor,6*T);p=solve(mu,kick,6*T);xb=b.sol(t).T;xp=p.sol(t).T
            # A well-resolved numerical reference orbit; distance is approximate.
            tree=cKDTree(xy);distance=tree.query(xp)[0]
            bt=b.t_events[0];pt=p.t_events[0];count=min(len(bt),len(pt));dt=pt[-count:]-bt[-count:]
            pair={'x_kick':.2,'anchor':anchor.tolist(),'final_same_time_separation':float(np.linalg.norm(xp[-1]-xb[-1])),
                  'final_distance_to_sampled_cycle':float(distance[-1]),'last_three_mean_crossing_time_shift':float(dt[-3:].mean()),
                  'reference_cycle_period':T,'reference_cycle_samples':len(xy)}
            arrays.update(pair_time=t,pair_base=xb,pair_perturbed=xp,pair_cycle_distance=distance)
    np.savez_compressed(out/'relaxation.npz',**arrays)
    save_json(out/'relaxation.json',{'model':'Representative Van der Pol / Lienard system; F=x^3/3-x. The supplied excerpt omits the preceding equation, so this explicit standard choice is not attributed as unseen source text.',
              'records':rows,'paired_mu10':pair,'interpretation':'Separate active nonlinear oscillator. Slow recovery and rapid jumps do not by themselves diagnose weakened interaction, fluid instability or loss of energy conservation in an isolated system.'})
    print('Fast-slow relaxation controls complete',flush=True)


def supplement():
    from analyze import RESULTS,plt,read,figsave,table,img
    r=read(RESULTS/'relaxation.json');z=np.load(RESULTS/'relaxation.npz')
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained');x=np.linspace(-2.3,2.3,500)
    axs[0].plot(x,cubic(x),'k--',lw=1,label='Cubic nullcline y=F(x)')
    for mu in [1.,10.,30.]:
        xy=z[f'cycle_{mu}'];axs[0].plot(xy[:,0],xy[:,1],label=f'μ={mu:g}')
    axs[0].scatter([-1,1],[2/3,-2/3],color='black',s=20)
    ymax=max(float(abs(z[f'cycle_{mu}'][:,1]).max()) for mu in [1.,10.,30.])
    axs[0].set(title='Slow branches and rapid jumps',xlabel='x',ylabel='y',ylim=(-1.08*ymax,1.08*ymax));axs[0].legend(fontsize=8)
    for mu in [1.,10.,30.]:
        t=z[f'cycle_time_{mu}'];axs[1].plot(t/t[-1],z[f'cycle_{mu}'][:,0],label=f'μ={mu:g}')
    axs[1].set(title='Waveform over one measured cycle',xlabel='Fraction of cycle',ylabel='x');axs[1].legend(fontsize=8)
    t=z['pair_time'];sep=np.linalg.norm(z['pair_perturbed']-z['pair_base'],axis=1)
    axs[2].semilogy(t,np.maximum(sep,1e-12),label='Same-time separation')
    axs[2].semilogy(t,np.maximum(z['pair_cycle_distance'],1e-12),label='Distance to sampled cycle')
    axs[2].set(title='Cycle recovery can retain a phase shift',xlabel='Time after x kick (μ=10)',ylabel='Phase-plane distance');axs[2].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'relaxation')
    q=r['paired_mu10']
    return '\n'.join([
        '<h2>7d. Added fast–slow control: relaxation oscillations</h2><p>The newest excerpt shows slow motion along a cubic nullcline and fast jumps between its outer branches. A standard representative model consistent with this picture is the following Van der Pol system. The preceding equation is absent from the supplied crop, so the cubic below is our explicit model choice.</p>',
        '<div class="eq">F(x)=x³/3−x; &nbsp; ẋ=μ[y−F(x)], &nbsp; ẏ=−x/μ, &nbsp; μ&gt;0.<br>Equivalently: ẍ−μ(1−x²)ẋ+x=0.</div>',
        '<p>At large μ, displacement away from y=F(x) drives rapid horizontal motion. On the outer branches, the trajectory evolves slowly until it reaches a fold near x=±1, where that slow branch ends and a fast jump follows. This change of speed occurs under fixed equations and fixed parameters. Integrating the outer slow branches gives T≈μ(3−2 ln2) as μ→∞; finite-μ corrections remain measurable.</p>',
        img('relaxation','Completed stiff ODE controls distinguish return toward a stable cycle from return to the same phase on that cycle.'),
        table(['μ','Measured period','Large-μ estimate','Relative estimate error','Slow-branch time fraction'],[[f'{x["mu"]:.0f}',f'{x["measured_period"]:.6f}',f'{x["large_mu_period"]:.6f}',f'{100*x["relative_asymptotic_period_error"]:.2f}%',f'{100*x["slow_branch_time_fraction"]:.1f}%'] for x in r['records']]),
        f'<p>After settling onto the μ=10 cycle, a copied state receives only an x kick of 0.2. After six cycles its same-time distance from the baseline is {q["final_same_time_separation"]:.5g}, while its distance to the sampled reference cycle is {q["final_distance_to_sampled_cycle"]:.5g}. The last three crossings retain a mean timing shift of {q["last_three_mean_crossing_time_shift"]:.5g}. The sampled-cycle distance has a finite discretization floor and is not an exact transverse exponent. An attracting cycle can therefore recover its shape while retaining a phase shift.</p>',
        '<p>This oscillator is active: for E=½(x²+ẋ²), Ė=μ(1−x²)ẋ². Small-amplitude motion receives energy through negative damping; larger-amplitude motion dissipates it. The numerical integration checks this balance. It must not be interpreted as a passive unforced fluid generating energy. To claim corresponding fast–slow fluid behavior would require identifying the relevant state variables, separated time scales and a validated reduced model.</p>',
        '<p>For the history question, x alone is an incomplete observation: y controls ẋ at a given x. Past x can help estimate that missing state. Once both x and y are observed, the first-order deterministic model needs no additional physical memory.</p>',
        '<details><summary>Supplied fast–slow excerpt</summary><img src="references/relaxation_reference.png" alt="User-supplied cubic nullcline and fast-slow relaxation oscillation excerpt"></details>'
    ])


if __name__=='__main__':main()
