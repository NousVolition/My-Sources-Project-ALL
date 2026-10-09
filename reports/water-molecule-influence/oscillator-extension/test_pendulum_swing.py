"""Additional controls from the supplied weak-nonlinearity exercises."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import ellipk
from test_oscillators import checked_solve, crossing

HERE=Path(__file__).resolve().parent


def run():
    out=HERE/'data';out.mkdir(parents=True,exist_ok=True)
    arrays={};pendulum=[];floquet=[]
    for amp in [.1,.2,.4]:
        f=lambda t,z: [z[1],-np.sin(z[0])]
        t=np.linspace(0,50,5001)
        sol=checked_solve(f,[0,50],[amp,0.],method='DOP853',rtol=1e-12,atol=1e-14,
                          max_step=.1,t_eval=t,events=crossing(1,direction=-1))
        period=float(np.mean(np.diff(sol.t_events[0])))
        exact_period=float(4*ellipk(np.sin(amp/2)**2))
        exact_omega=2*np.pi/exact_period
        approximate_omega=1-amp*amp/16
        relative_period_error=abs(period-exact_period)/exact_period
        assert relative_period_error<1e-9
        pendulum.append({'amplitude':amp,'period':period,'elliptic_integral_period':exact_period,
                         'relative_period_error':relative_period_error,'exact_frequency':exact_omega,
                         'averaged_frequency':approximate_omega,
                         'averaging_error':float(abs(approximate_omega-exact_omega))})
        arrays[f'pendulum_{amp:g}_t'],arrays[f'pendulum_{amp:g}_z']=t,sol.y.T
    orders=np.log2(np.array([r['averaging_error'] for r in pendulum[1:]])/
                   np.array([r['averaging_error'] for r in pendulum[:-1]]))
    assert np.all((orders>3.8)&(orders<4.2))

    # The supplied epsilon=2 exercise has x_dot^3 (velocity damping), not x^3.
    damping=[]
    for eps in [.1,.2,2.]:
        amp=1.;t=np.linspace(0,50,10001)
        def f(t,z):
            x,v,loss=z
            return [v,-x-eps*v**3,eps*v**4]
        sol=checked_solve(f,[0,50],[amp,0.,0.],method='DOP853',rtol=1e-12,atol=1e-14,
                          max_step=.05,t_eval=t)
        radius=amp/np.sqrt(1+3*eps*amp*amp*t/4)
        approximation=radius*np.cos(t)
        energy=(sol.y[0]**2+sol.y[1]**2)/2
        balance=float(np.max(abs(energy+sol.y[2]-.5)))
        assert balance<1e-9 and np.max(np.diff(energy))<1e-10
        damping.append({'epsilon':eps,'initial':[amp,0.],'end_time':50.,
                         'maximum_position_error':float(np.max(abs(sol.y[0]-approximation))),
                         'RMS_position_error':float(np.sqrt(np.mean((sol.y[0]-approximation)**2))),
                         'maximum_energy_balance_residual':balance,
                         'final_energy':float(energy[-1]),'initial_energy':.5})
        arrays[f'cubic_damping_{eps:g}_t']=t
        arrays[f'cubic_damping_{eps:g}_z']=sol.y[:2].T
        arrays[f'cubic_damping_{eps:g}_approximation']=approximation

    # Linearized swing over one forcing period pi: M maps initial perturbations to final ones.
    for eps in [.1,.05,.025]:
        for gamma in [-.75,-.25,0.,.25,.75]:
            def matrix_equation(t,q):
                M=q.reshape(2,2)
                A=np.array([[0.,1.],[-(1+eps*gamma+eps*np.cos(2*t)),0.]])
                return (A@M).ravel()
            sol=checked_solve(matrix_equation,[0,np.pi],np.eye(2).ravel(),
                              method='DOP853',rtol=2e-13,atol=1e-15,max_step=.05)
            M=sol.y[:,-1].reshape(2,2);eig=np.linalg.eigvals(M)
            exponent=float(np.log(np.max(abs(eig)))/np.pi)
            prediction=float(eps*np.sqrt(max(0.,1-4*gamma*gamma))/4)
            assert abs(np.linalg.det(M)-1)<1e-10
            if abs(gamma)<.5:
                assert exponent>0 and abs(exponent-prediction)/prediction<.06
            else:
                assert abs(exponent)<1e-10
            floquet.append({'epsilon':eps,'gamma':gamma,'growth_rate_per_time':exponent,
                            'averaged_growth_rate':prediction,'determinant':float(np.linalg.det(M)),
                            'multipliers':[[float(q.real),float(q.imag)] for q in eig]})
    # Exact rest is invariant. A nonzero seed in the resonance band grows.
    eps=.1;gamma=0.;t=np.linspace(0,160,8001)
    def nonlinear_swing(t,z):
        return [z[1],-(1+eps*gamma+eps*np.cos(2*t))*np.sin(z[0])]
    trials=[]
    for seed in [0.,1e-4]:
        sol=checked_solve(nonlinear_swing,[0,160],[seed,0.],method='DOP853',
                          rtol=1e-11,atol=1e-14,max_step=.1,t_eval=t)
        radius=np.linalg.norm(sol.y,axis=0)
        trials.append({'initial_x':seed,'initial_v':0.,'final_radius':float(radius[-1]),
                        'maximum_radius':float(np.max(radius))})
        if seed==0:
            assert np.max(radius)==0.
        else:
            assert radius[-1]>10*seed and np.max(radius)<.02
        arrays[f'swing_{seed:g}_t'],arrays[f'swing_{seed:g}_z']=t,sol.y.T
    result={'pendulum':pendulum,'pendulum_averaging_error_orders':orders.tolist(),
            'cubic_velocity_damping':{'equation':'x_ddot+epsilon*x_dot^3+x=0',
              'averaged_amplitude':'dr/dt=-3*epsilon*r^3/8; r(t)=a/sqrt(1+3*epsilon*a^2*t/4)',
              'position_approximation':'x(t)=r(t)*cos(t)',
              'energy_identity':'E=(x^2+v^2)/2; dE/dt=-epsilon*v^4',
              'solver':'DOP853; rtol=1e-12, atol=1e-14, max_step=0.05; integrate the loss integral alongside x,v.',
              'rows':damping,'limitation':'The epsilon=2 case is outside a controlled small-parameter expansion. The leading approximation matches x(0)=a, but its envelope derivative gives an O(epsilon) mismatch to v(0)=0; no uniform high-accuracy claim is made.'},
            'swing':{'equation':'x_ddot+(1+epsilon*gamma+epsilon*cos(2*t))*sin(x)=0',
              'linearized_equation':'x_ddot+(1+epsilon*gamma+epsilon*cos(2*t))*x=0',
              'averaged_equations':'dr/dT=r*sin(2*phi)/4; dphi/dT=(gamma+cos(2*phi)/2)/2; T=epsilon*t',
              'phase_convention':'x=r*cos(t+phi), x_dot=-r*sin(t+phi)',
              'leading_instability_band':'abs(gamma)<1/2; finite-epsilon boundaries can shift',
              'floquet_method':'Integrate the 2x2 fundamental matrix over pi with DOP853, rtol=2e-13, atol=1e-15; growth=log(max(abs(eigenvalues)))/pi.',
              'floquet':floquet,'nonlinear_seed_test':trials,'seed_test_epsilon':.1,'seed_test_gamma':0.,'seed_test_end':160.},
            'passed':True}
    (out/'exercises.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(out/'exercises.npz',**arrays)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':run()
