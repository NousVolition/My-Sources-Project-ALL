"""Controlled weak/strong van der Pol and conservative Duffing benchmarks.

All quantities are dimensionless. No random inputs; no water or biological model.
Run directly to save trajectories, event measurements, and acceptance results.
"""
from pathlib import Path
import json, platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent


def checked_solve(f, interval, initial, **kwargs):
    sol = solve_ivp(f, interval, initial, **kwargs)
    assert sol.success, sol.message
    assert np.isfinite(sol.y).all()
    return sol


def crossing(component, value=0., direction=1):
    def event(t, z):
        return z[component] - value
    event.direction = direction
    return event


def rk4_final(f, initial, duration, h):
    z = np.array(initial, dtype=float)
    for _ in range(round(duration/h)):
        k1 = h*f(0, z)
        k2 = h*f(0, z+k1/2)
        k3 = h*f(0, z+k2/2)
        k4 = h*f(0, z+k3)
        z += (k1+2*k2+2*k3+k4)/6
    return z


def run():
    out = HERE/'data'
    out.mkdir(parents=True, exist_ok=True)
    arrays, weak, relaxation, duffing = {}, [], [], []
    protocol = {
        'date': '2026-10-09', 'python': platform.python_version(),
        'numpy': np.__version__, 'scipy': scipy.__version__,
        'van_der_pol': 'x_dot=v; v_dot=mu*(1-x^2)*v-x',
        'weak': {'mu': [.05, .1, .2], 'initial_x': [.2, 1., 3.],
                 'initial_v': 0., 'slow_time_end_mu_t': 16.,
                 'solver': 'DOP853', 'rtol': 1e-10, 'atol': 1e-12,
                 'max_step': .2, 'observable': 'positive maxima located by v=0, direction=-1'},
        'relaxation': {'mu': [5., 10., 20., 40.], 'initial_xy': [2., 0.],
                       'coordinates': 'F(x)=x^3/3-x; y=F(x)+v/mu; x_dot=mu*(y-F(x)); y_dot=-x/mu',
                       'time_end': '12*mu+40', 'solver': 'Radau',
                       'rtol': 1e-9, 'atol': 1e-11,
                       'refined_rtol': 1e-11, 'refined_atol': 1e-13,
                       'slow_interval': 'descending positive branch from x=1.8 to x=1.2',
                       'fast_interval': 'descending jump from x=0.5 to x=-0.5'},
        'duffing': {'equation': 'x_dot=v; v_dot=-x-epsilon*x^3',
                    'epsilon': .1, 'initial_amplitudes': [.5, 1., 2.],
                    'initial_v': 0., 'time_end': 200., 'solver': 'DOP853',
                    'rtol': 1e-11, 'atol': 1e-13, 'max_step': .1},
        'acceptance': {
            'weak_limit_cycle': 'last positive peaks within 0.01 of 2; initial-state spread below 0.002 at mu*t=16',
            'weak_averaging': 'maximum peak-envelope discrepancy at mu=0.05 below 0.03',
            'RK4_refinement': 'observed endpoint order 3.5 to 4.5 against DOP853',
            'relaxation_numerics': 'period changes by less than 1e-6 relative under tighter tolerance; last periods differ by less than 1e-6 relative',
            'relaxation_asymptotics': 'period relative error decreases with mu; slow time/mu within 0.02 and fast time*mu within 0.15 of asymptotic constants at mu=40',
            'duffing_energy': 'maximum relative energy error below 1e-7; period quadrature error below 1e-7 relative; peak amplitude error below 1e-7'
        }
    }
    (out/'protocol.json').write_text(json.dumps(protocol, indent=2), encoding='utf-8')

    for mu in protocol['weak']['mu']:
        f = lambda t,z: [z[1], mu*(1-z[0]**2)*z[1]-z[0]]
        for r0 in protocol['weak']['initial_x']:
            end = 16/mu
            sol = checked_solve(f, [0, end], [r0, 0.], method='DOP853',
                                rtol=1e-10, atol=1e-12, max_step=.2,
                                events=crossing(1, direction=-1), dense_output=True)
            peaks_t, peaks = sol.t_events[0], sol.y_events[0][:, 0]
            predicted = 2/np.sqrt(1+(4/r0**2-1)*np.exp(-mu*peaks_t))
            error = float(np.max(abs(peaks-predicted)))
            tail_period = float(np.mean(np.diff(peaks_t[-5:])))
            assert abs(peaks[-1]-2) < .01
            if mu == .05:
                assert error < .03
            key = f'weak_{mu:g}_{r0:g}'
            t = np.linspace(0, end, 5001)
            arrays[key+'_t'], arrays[key+'_z'] = t, sol.sol(t).T
            arrays[key+'_peak_t'], arrays[key+'_peaks'] = peaks_t, peaks
            weak.append({'mu': mu, 'initial_amplitude': r0, 'last_peak': float(peaks[-1]),
                         'period_last_four': tail_period,
                         'max_peak_envelope_error': error,
                         'last_peak_slow_time': float(mu*peaks_t[-1])})
        values = [row['last_peak'] for row in weak if row['mu']==mu]
        assert max(values)-min(values) < .002

    # Independent implementation control for the explicitly shown RK4 formula.
    f = lambda t,z: np.array([z[1], .1*(1-z[0]**2)*z[1]-z[0]])
    reference = checked_solve(f, [0,20], [.2,0.], method='DOP853', rtol=2e-13, atol=1e-15)
    steps = [.1, .05, .025]
    errors = [float(np.linalg.norm(rk4_final(f,[.2,0.],20,h)-reference.y[:,-1])) for h in steps]
    orders = np.log2(np.array(errors[:-1])/errors[1:])
    assert np.all((orders>3.5)&(orders<4.5))

    F = lambda x: x**3/3-x
    slow_constant = .9-np.log(1.5)
    fast_constant = quad(lambda x: 1/(2/3+F(x)), -.5, .5, epsabs=1e-13)[0]
    period_constant = 3-2*np.log(2)
    for mu in protocol['relaxation']['mu']:
        f = lambda t,z: np.array([mu*(z[1]-F(z[0])), -z[0]/mu])
        end = 12*mu+40
        periods = []
        for tol in [1e-9, 1e-11]:
            sol = checked_solve(f, [0,end], [2.,0.], method='Radau',
                                rtol=tol, atol=tol*.01, dense_output=True,
                                events=crossing(0, direction=1))
            events = sol.t_events[0]
            assert len(events)>=5
            periods.append(float(events[-1]-events[-2]))
        period = periods[-1]
        convergence = float(abs(periods[-1]-periods[-2])/period)
        cycle_change = float(abs((events[-2]-events[-3])-period)/period)
        assert convergence < 1e-6 and cycle_change < 1e-6
        # The last full crossing-to-crossing cycle contains each descending threshold once.
        start, stop = events[-2], events[-1]
        probe = np.linspace(start, stop, 20001)
        xprobe = sol.sol(probe)[0]
        def descending_time(level):
            idx = np.flatnonzero((xprobe[:-1]>level)&(xprobe[1:]<=level))
            assert len(idx)==1
            i = idx[0]
            return brentq(lambda t: sol.sol(t)[0]-level, probe[i], probe[i+1], xtol=1e-12)
        t18, t12, tp, tm = [descending_time(level) for level in [1.8,1.2,.5,-.5]]
        slow_time, fast_time = t12-t18, tm-tp
        slow_t = np.linspace(t18,t12,2001)
        slow_x, slow_y = sol.sol(slow_t)
        residual = float(np.median(abs(slow_y-F(slow_x))))
        t = np.linspace(start,stop,12001)
        z = sol.sol(t).T
        arrays[f'relax_{mu:g}_t'], arrays[f'relax_{mu:g}_z'] = t-start,z
        # Save a dense trace specifically through the short jump for plotting and audit.
        jt = np.linspace(tp,tm,1001)
        arrays[f'jump_{mu:g}_t'], arrays[f'jump_{mu:g}_z'] = jt-tp,sol.sol(jt).T
        arrays[f'slow_{mu:g}_t'], arrays[f'slow_{mu:g}_z'] = slow_t-t18,sol.sol(slow_t).T
        row = {'mu':mu, 'period':period, 'period_asymptote':float(period_constant*mu),
               'period_relative_asymptotic_error':float(abs(period-period_constant*mu)/(period_constant*mu)),
               'tolerance_period_relative_difference':convergence,
               'last_cycle_relative_difference':cycle_change,
               'slow_interval_duration':float(slow_time), 'fast_interval_duration':float(fast_time),
               'slow_over_mu':float(slow_time/mu), 'fast_times_mu':float(fast_time*mu),
               'slow_fast_ratio':float(slow_time/fast_time),
               'median_slow_nullcline_residual':residual,
               'mu_squared_nullcline_residual':float(mu*mu*residual)}
        if mu==10.:
            # Different solver family; compare an entire late cycle, not just an endpoint.
            independent = checked_solve(f,[0,end],[2.,0.],method='DOP853',rtol=1e-11,atol=1e-13,
                                        dense_output=True,events=crossing(0,direction=1))
            other_period = float(np.diff(independent.t_events[0])[-1])
            row['DOP853_period_relative_difference'] = abs(other_period-period)/period
            row['DOP853_cycle_max_distance'] = float(np.max(np.linalg.norm(independent.sol(t).T-z,axis=1)))
            assert row['DOP853_period_relative_difference']<1e-7
            assert row['DOP853_cycle_max_distance']<1e-5
        relaxation.append(row)
    p_errors = [row['period_relative_asymptotic_error'] for row in relaxation]
    assert np.all(np.diff(p_errors)<0)
    assert abs(relaxation[-1]['slow_over_mu']-slow_constant)<.02
    assert abs(relaxation[-1]['fast_times_mu']-fast_constant)<.15

    eps = .1
    for amp in [.5,1.,2.]:
        f = lambda t,z: [z[1], -z[0]-eps*z[0]**3]
        t = np.linspace(0,200,10001)
        sol = checked_solve(f,[0,200],[amp,0.],method='DOP853',rtol=1e-11,atol=1e-13,
                            max_step=.1,t_eval=t,events=crossing(1,direction=-1))
        x,v = sol.y
        energy = v*v/2+x*x/2+eps*x**4/4
        energy_error = float(np.max(abs(energy-energy[0]))/energy[0])
        # x=A*sin(theta) removes the endpoint singularity in the energy quadrature.
        exact_period = 4*quad(lambda theta: 1/np.sqrt(1+eps*amp**2/2*(1+np.sin(theta)**2)),
                              0,np.pi/2,epsabs=1e-13,epsrel=1e-13)[0]
        period = float(np.mean(np.diff(sol.t_events[0][-6:])))
        peak_error = float(np.max(abs(sol.y_events[0][:,0]-amp)))
        period_error = float(abs(period-exact_period)/exact_period)
        assert energy_error < 1e-7 and period_error < 1e-7 and peak_error < 1e-7
        duffing.append({'initial_amplitude':amp,'epsilon':eps,'period':period,
                        'quadrature_period':float(exact_period),'period_relative_error':period_error,
                        'averaged_frequency':float(1+3*eps*amp**2/8),
                        'measured_frequency':float(2*np.pi/period),
                        'maximum_relative_energy_error':energy_error,'maximum_peak_amplitude_error':peak_error})
        arrays[f'duffing_{amp:g}_t'],arrays[f'duffing_{amp:g}_z'] = t,sol.y.T

    result = {'weak':weak,'rk4':{'steps':steps,'errors':errors,'orders':orders.tolist()},
              'relaxation':relaxation,'asymptotic_constants':{'period_over_mu':float(period_constant),
               'slow_over_mu':float(slow_constant),'fast_times_mu':float(fast_constant)},
              'duffing':duffing,'passed':True}
    np.savez_compressed(out/'trajectories.npz',**arrays)
    (out/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2), flush=True)


if __name__=='__main__':
    run()
