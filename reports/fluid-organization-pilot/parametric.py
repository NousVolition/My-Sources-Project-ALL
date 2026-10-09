"""Parametrically driven pendulum/swing control from the supplied excerpt.

x'' + [1 + eps*gamma + eps*cos(2t)]*sin(x) = 0.
The modulation is an explicit external energy source. Exact zero state is
invariant; instability amplifies a nonzero seed, not an exactly absent motion.
"""
import numpy as np
from scipy.integrate import solve_ivp
from core import ROOT,save_json


def monodromy(eps,gamma):
    def rhs(t,h):
        A=np.array([[0.,1.],[-(1+eps*gamma+eps*np.cos(2*t)),0.]])
        return (A@h.reshape(2,2)).ravel()
    sol=solve_ivp(rhs,(0,np.pi),np.eye(2).ravel(),method='DOP853',rtol=1e-12,atol=1e-14,max_step=.02)
    M=sol.y[:,-1].reshape(2,2);multipliers=np.linalg.eigvals(M)
    return M,multipliers


def swing(eps,gamma,kick,times,modulation=True):
    def rhs(t,z):
        x,v,W=z;k=1+eps*gamma+(eps*np.cos(2*t) if modulation else 0.)
        work=-2*eps*np.sin(2*t)*(1-np.cos(x)) if modulation else 0.
        return [v,-k*np.sin(x),work]
    sol=solve_ivp(rhs,(times[0],times[-1]),[kick,0.,0.],t_eval=times,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y.T


def main():
    eps=.1;times=np.linspace(0,200,10001);rows=[];arrays={'time':times}
    for gamma in [-1.,-.75,-.25,0.,.25,.75,1.]:
        M,eig=monodromy(eps,gamma);rate=float(np.log(max(abs(eig)))/np.pi)
        leading=eps/4*np.sqrt(max(1-4*gamma**2,0))
        rows.append({'gamma':gamma,'floquet_growth_rate':rate,'leading_growth_rate':float(leading),'determinant':float(np.linalg.det(M)),
                     'multiplier_real_imag':[[float(e.real),float(e.imag)] for e in eig]})
    checks=[]
    for gamma,kick,mod in [(0.,0.,True),(0.,.001,True),(-.75,.001,True),(.75,.001,True),(0.,.001,False),(0.,.01,True)]:
        z=swing(eps,gamma,kick,times,mod);key=f'g{gamma}_x{kick}_mod{int(mod)}';arrays[key]=z
        k=1+eps*gamma+(eps*np.cos(2*times) if mod else 0.)
        H=.5*z[:,1]**2+k*(1-np.cos(z[:,0]));res=H-H[0]-z[:,2]
        checks.append({'gamma':gamma,'kick':kick,'modulation':mod,'max_abs_angle':float(abs(z[:,0]).max()),
                       'max_energy_work_residual':float(abs(res).max())})
    # Leading averaged equations use slow time T=eps*t; phi0=0 for x0>0,v0=0.
    def avg(T,z):
        r,phi=z;return [.25*r*np.sin(2*phi),.25*np.cos(2*phi)]
    av=solve_ivp(avg,(0,eps*times[-1]),[.001,0.],t_eval=eps*times,rtol=1e-11,atol=1e-13,method='DOP853')
    arrays['averaged_radius']=av.y[0];arrays['averaged_phase']=av.y[1]
    out=ROOT/'results';np.savez_compressed(out/'parametric.npz',**arrays)
    save_json(out/'parametric.json',{'epsilon':eps,'periodic_coefficient_period':float(np.pi),'floquet':rows,'trajectory_controls':checks,
              'leading_instability_band':'|gamma| < 1/2; finite-epsilon Floquet results determine the actual linear behavior',
              'scope':'Prescribed external modulation; no damping, noise or forcing independent of x. Exact rest remains exact rest.'})
    print('Parametric swing: rest, seeded amplification, detuning, forcing-removal and energy-work controls complete',flush=True)


def supplement():
    from analyze import RESULTS,plt,read,figsave,table,img
    r=read(RESULTS/'parametric.json');z=np.load(RESULTS/'parametric.npz');t=z['time']
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    rr=r['floquet'];axs[0].plot([q['gamma'] for q in rr],[q['floquet_growth_rate'] for q in rr],'o-',label='Numerical Floquet rate')
    g=np.linspace(-1,1,400);axs[0].plot(g,.025*np.sqrt(np.maximum(1-4*g*g,0)),'--',label='Leading averaged rate')
    axs[0].set(title='Detuning controls amplification',xlabel='γ',ylabel='Linear growth rate per time');axs[0].legend(fontsize=8)
    for key,label in [('g0.0_x0.001_mod1','Near resonance; small seed'),('g0.75_x0.001_mod1','Detuned'),('g0.0_x0.001_mod0','Modulation removed')]:
        state=z[key];rad=np.linalg.norm(state[:,:2],axis=1);axs[1].semilogy(t,rad,label=label,lw=.8)
    axs[1].semilogy(t,z['averaged_radius'],'k--',label='Leading averaged envelope')
    axs[1].set(title='A seed is required in this model',xlabel='Time',ylabel='Phase-plane radius');axs[1].legend(fontsize=7)
    state=z['g0.0_x0.01_mod1'];H=.5*state[:,1]**2+(1+.1*np.cos(2*t))*(1-np.cos(state[:,0]))
    axs[2].plot(t,H-H[0],label='Energy change');axs[2].plot(t,state[:,2],'--',label='Integrated modulation work')
    axs[2].set(title='Growth has an explicit energy source',xlabel='Time',ylabel='Energy / work');axs[2].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'parametric_swing')
    return '\n'.join([
        '<h2>7f. Added averaging, amplitude and parametric-forcing controls</h2><h3>Pendulum frequency correction</h3><p>Expanding the pendulum result gives ω≈1−α²/16 for small amplitude α, equivalently T≈2π(1+α²/16). This agrees with the exact elliptic-integral limit already checked. The crop begins with an incomplete preceding exercise; no missing equation from that exercise is inferred.</p>',
        '<h3>Why the weak Van der Pol radius is approximately 2</h3><p>For (ẋ,v̇)=(v,−x−ε(x²−1)v), the planar divergence is ε(1−x²). The net flux through a periodic orbit is zero because the vector field is tangent to it. Approximating the weakly nonlinear orbit by a circle of radius a, Green’s theorem gives 0=ε∬(1−x²)dA=επa²(1−a²/4), so the nonzero leading radius is a≈2. The circular approximation, rather than the divergence theorem, is the approximation in this argument.</p>',
        '<h3>A swing driven by periodic parameter modulation</h3><div class="eq">ẍ+[1+εγ+εcos(2t)]sinx=0.<br>For small x, slow time T=εt, x≈r(T)cos[t+φ(T)]:<br>r′=¼r sin(2φ), &nbsp; φ′=½[γ+½cos(2φ)].</div>',
        '<p>The model can amplify a small displacement near parametric resonance. Its leading small-ε instability band is |γ|&lt;1/2, with growth rate (ε/4)√(1−4γ²) per ordinary time. The completed ε=0.1 experiment computes the exact linear monodromy over the forcing period π to check this approximation; finite-ε boundaries need not equal the leading ones.</p>',
        img('parametric_swing','Exact rest, a controlled nonzero seed, detuning and removal of modulation distinguish invariance from parametric instability.'),
        table(['γ','Numerical Floquet rate','Leading averaged rate','det(monodromy)'],[[q['gamma'],f'{q["floquet_growth_rate"]:.7f}',f'{q["leading_growth_rate"]:.7f}',f'{q["determinant"]:.9f}'] for q in rr]),
        '<p>Starting exactly at x=0, ẋ=0 leaves the state exactly at rest. The modulation multiplies sinx and supplies no additive kick at the origin; this idealized swing needs a push or another nonzero perturbation. A 0.001 displacement grows near resonance, stays bounded in the checked detuned cases, and does not grow when modulation is removed. A larger 0.01 seed shows nonlinear departure from the leading small-amplitude envelope.</p>',
        '<p>The instantaneous energy H(t)=½ẋ²+[1+εγ+εcos(2t)](1−cosx) obeys Ḣ=−2εsin(2t)(1−cosx). Numerical energy change matches integrated work from the prescribed modulation. Amplification therefore has an identified external source. In a fluid experiment, comparable claims would require measuring periodic forcing work and the relevant response mode; this oscillator control does not establish fluid resonance.</p>',
        '<details><summary>Supplied averaging and swing excerpt</summary><img src="references/parametric_reference.png" alt="User-supplied pendulum averaging, Green theorem and parametric swing exercises"></details>'
    ])


if __name__=='__main__':main()
