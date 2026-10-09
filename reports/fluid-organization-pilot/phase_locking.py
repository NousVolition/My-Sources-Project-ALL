"""Adler phase-locking control motivated by the supplied Figure 4.5.1 excerpt.

phi = stimulus phase - oscillator phase, dphi/dt = detuning - A*sin(phi).
This is a separate reduced oscillator, not a fluid or particle interaction law.
"""
import numpy as np
from scipy.integrate import solve_ivp
from core import ROOT,save_json


def evolve(detuning,A,initial,times,rtol=1e-10):
    return solve_ivp(lambda t,p:detuning-A*np.sin(p),(times[0],times[-1]),np.atleast_1d(initial).astype(float),t_eval=times,rtol=rtol,atol=rtol/100,method='DOP853').y.T


def main():
    out=ROOT/'results';out.mkdir(exist_ok=True)
    # Long independent initial phases estimate the mean slipping frequency.
    times=np.linspace(0,400,4001);inits=np.linspace(-2.7,2.7,8)
    records=[];tracks={}
    for mu in [0.,.5,.9,1.,1.2,1.5]:
        x=evolve(mu,1.,inits,times);fine=evolve(mu,1.,inits,times,1e-12)
        measured=(x[-1]-x[2000])/200
        theory=0. if abs(mu)<=1 else np.sqrt(mu*mu-1)
        # Period from complete 2pi crossings avoids finite-window endpoint bias.
        periods=[]
        if mu>1:
            for q in range(8):
                levels=np.arange(np.ceil(x[2000,q]/(2*np.pi)),np.floor(x[-1,q]/(2*np.pi))+1)*2*np.pi
                crosses=np.interp(levels,x[:,q],times)
                periods.extend(np.diff(crosses).tolist())
        period_rate=float(2*np.pi/np.mean(periods)) if periods else 0.
        anchor=evolve(mu,1.,[.1],[0,20])[-1,0]
        short=np.linspace(0,20,401);paired=evolve(mu,1.,[anchor,anchor+.05],short)
        tracks[f'mu{mu}']=x[::4];tracks[f'pair{mu}']=paired
        records.append({'mu':mu,'stable_phase_rad':float(np.arcsin(mu)) if mu<1 else None,
                        'linear_recovery_rate':float(np.sqrt(1-mu*mu)) if mu<1 else 0.,
                        'predicted_mean_slip_rate':float(theory),'window_slip_mean':float(measured.mean()),
                        'window_slip_min':float(measured.min()),'window_slip_max':float(measured.max()),
                        'complete_cycle_slip_rate':period_rate,
                        'max_tolerance_refinement_phase_error':float(abs(x-fine).max())})
    # Same physical detuning, coupling removed/restored across four controls.
    coupling=[]
    for A in [0.,.5,1.,2.]:
        x=evolve(1.,A,inits,times)
        coupling.append({'detuning':1.,'A':A,'observed_window_drift':float(np.mean((x[-1]-x[2000])/200)),
                         'predicted_drift':float(np.sqrt(1-A*A)) if A<1 else 0.})
    np.savez_compressed(out/'phase_locking.npz',time=times[::4],pair_time=short,**tracks)
    save_json(out/'phase_locking.json',{'equation':'dphi/dt = Delta - A sin(phi); mu=Delta/A, tau=A t when A>0',
              'records':records,'coupling_controls':coupling,'initial_phases':inits.tolist(),
              'interpretation':'Strict |Delta|<A gives robust stable phase lock except the unstable fixed point. Equality is the nonhyperbolic saddle-node threshold, not robust locking. Above threshold phases slip with nonuniform instantaneous rate. No fluid phase was identified or fitted.'})
    print('Phase-locking controls complete',flush=True)


if __name__=='__main__':main()
