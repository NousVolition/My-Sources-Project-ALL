"""The user's rotating-hoop example: retain inertia, then compare overdamping.

The book's tangential force convention is m*r*phi_tt+b*phi_t=m*g*f(phi).
Time s=t*sqrt(g/r), q=r*omega^2/g, beta=b/(m*sqrt(g*r)).
phi_ss + beta*phi_s = sin(phi)*(q*cos(phi)-1).
In slow time tau=s/beta: beta^-2*phi_tautau + phi_tau = f(phi).
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from damped_leapfrog import leapfrog


def main():
    out = Path(__file__).resolve().parent/'data'
    out.mkdir(exist_ok=True)
    arrays, rows = {}, []
    for q in [.5, 1., 2.]:
        f = lambda phi: np.sin(phi)*(q*np.cos(phi)-1)
        for beta in [2., 20.]:
            # Same displacement and zero physical initial velocity.
            s, phi, velocity = leapfrog(f, .2, 0., beta, .005, round(12*beta/.005))
            tau = s/beta
            first = solve_ivp(lambda t, y: f(y), [0, 12], [.2], t_eval=tau,
                              rtol=1e-10, atol=1e-12).y[0]
            key = f'q{q:g}_b{beta:g}'
            keep = np.unique(np.r_[np.arange(0, len(tau), max(1, len(tau)//1200)), len(tau)-1])
            arrays[key+'_tau'], arrays[key+'_phi'], arrays[key+'_overdamped'] = tau[keep], phi[keep], first[keep]
            rows.append({'q': q, 'beta': beta, 'max_angle_error_vs_overdamped': float(np.max(abs(phi-first))),
                         'final_phi': float(phi[-1]), 'overdamped_final_phi': float(first[-1])})
    # Meaningful approximation check, not an assertion of monotonicity for all starts.
    for q in [.5, 1., 2.]:
        matches = [r for r in rows if r['q'] == q]
        assert matches[1]['max_angle_error_vs_overdamped'] < matches[0]['max_angle_error_vs_overdamped']/30
    # Phase portrait in the book's slow variables: Omega=dphi/dtau.
    phase_rows = []
    for beta in [2., 20.]:
        f = lambda phi: np.sin(phi)*(2*np.cos(phi)-1)
        for index, phi0 in enumerate([-5., -3., -.25, .25, 3., 5.]):
            omega0 = 2.5 if index % 2 == 0 else -2.5
            s, phi, velocity = leapfrog(f, phi0, omega0/beta, beta, .005, round(5*beta/.005))
            tau = s/beta
            omega = beta*velocity
            keep = np.unique(np.r_[np.arange(min(100, len(tau))), np.arange(100, len(tau), 10), len(tau)-1])
            key = f'phase_b{beta:g}_i{index}'
            arrays[key+'_phi'], arrays[key+'_omega'], arrays[key+'_tau'] = phi[keep], omega[keep], tau[keep]
            after = tau >= .1
            phase_rows.append({'beta': beta, 'phi0': phi0, 'Omega0': omega0,
                               'max_abs_Omega_minus_f_after_tau_0_1': float(np.max(abs(omega[after]-f(phi[after]))))})
    f = lambda phi: np.sin(phi)*(2*np.cos(phi)-1)
    _, coarse, _ = leapfrog(f, .2, 0., 20., .005, 48000)
    _, fine, _ = leapfrog(f, .2, 0., 20., .0025, 96000)
    convergence = float(np.max(abs(coarse-fine[::2])))
    assert convergence < 1e-5
    result = {'rows': rows, 'phase_rows': phase_rows, 'half_step_angle_difference_q2_beta20': convergence, 'initial_phi': .2, 'initial_physical_velocity': 0,
              'ds': .005, 'slow_time_end': 12,
              'equilibria': 'phi=0,pi always; phi=+/-acos(1/q) when q>1 (mod 2*pi).',
              'stability': 'For beta>0: bottom stable for q<1 and unstable for q>1; top unstable; off-center branches stable. At q=1 bottom is nonhyperbolic but nonlinearly stable.',
              'scaling': 'Using book force convention mr*phi_tt+b*phi_t=mg*f(phi): T=b/(mg), tau=t/T, beta=b/(m*sqrt(gr)), epsilon=1/beta^2=m^2*g*r/b^2.',
              'phase_equations': 'phi_tau=Omega; epsilon*Omega_tau=f(phi)-Omega. C: Omega=f(phi) is a nullcline, not exactly invariant for finite epsilon.',
              'qualification': 'Large beta at fixed q and after the initial velocity relaxation supports the overdamped reduction. Damping alone is not a proof that inertia can be dropped.'}
    (out/'hoop.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    np.savez_compressed(out/'hoop.npz', **arrays)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
