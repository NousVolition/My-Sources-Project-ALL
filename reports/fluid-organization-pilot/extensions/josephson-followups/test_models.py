import numpy as np
from scipy.integrate import solve_ivp
from circuit import rhs,to_state,stable_velocity,stable_rhs,synchronized_voltage
from noise import mobility,heun,drift
from long_cycles import cycles,continued_voltage
from capacitance import segment
from calibrate import parameters

def test_log_coordinates_match_direct_physics():
    rng=np.random.default_rng(22);phi=rng.uniform(-2,2,(50,2));s,sign=to_state(phi)
    np.testing.assert_allclose(stable_velocity(s,sign,1.2,.5),rhs(phi,1.2,.5),atol=1e-14)
def test_load_mobility_and_fdt_positive():
    for alpha in (0,.5,2):
        m=mobility(alpha);np.testing.assert_allclose(m@(np.eye(2)+alpha*np.ones((2,2))),np.eye(2),atol=1e-14)
        assert np.linalg.eigvalsh(m).min()>0
        chol=np.linalg.cholesky(m);np.testing.assert_allclose(chol@chol.T,m)
def test_noise_free_drift_exact():
    phi=np.array([.2,-.3]);inc=np.array([.01,.02]);zero=np.zeros(2)
    np.testing.assert_allclose(heun(phi,1.1,.5,zero,.01,inc),phi+.01*drift(phi,1.1,.5,zero)+inc)
def test_mismatch_current_balance():
    p=np.array([.7,-.2]);critical=np.array([.98,1.02]);v=drift(p,1.1,.5,critical)
    np.testing.assert_allclose(v+.5*v.sum()+critical*np.sin(p),[1.1,1.1])
def test_cycle_exact_uncoupled_voltage_and_return():
    s,sign=to_state(np.array([[.1,.4],[-1.,.5]]));sol,t,ret,perr=cycles(s,[1.02,1.2],0)
    np.testing.assert_allclose(4*np.pi/t,2*np.sqrt(np.array([1.02,1.2])**2-1),rtol=1e-8)
    assert ret.max()<1e-8 and perr.max()<1e-7
def test_periodic_continuation_matches_direct_short():
    s,_=to_state(np.array([[.1,1.2]]));sol,t,_,_=cycles(s,1.2,.5)
    v=continued_voltage(sol,t,[50])[0,0]
    q=solve_ivp(lambda t,y:stable_rhs(y,1.2,.5),(0,50),s[0],rtol=1e-11,atol=1e-12,method='DOP853')
    assert abs(v-2*(q.y[0,-1]-s[0,0])/50)<1e-8
def test_capacitance_static_equilibrium():
    s=np.array([[np.arcsin(.4),0.]])
    q,v,r=segment(s,.4,2,.5,5,5)
    np.testing.assert_allclose(q[-1,0,:2],s[0],atol=1e-12);assert abs(v[0])<1e-12
def test_capacitance_energy_balance():
    q,v,r=segment(np.array([[.2,.3]]),1.2,2,.5,10,10)
    assert abs(r[0])<1e-9
def test_calibration_dimensions_and_mismatch():
    p=parameters(1e-6,1.02e-6,100,200,.1,1e-12)
    assert p['alpha']==.5 and p['beta']>0 and p['Gamma']>0
    np.testing.assert_allclose(sum(p['critical_current_coefficients']),2)
    from scipy.constants import hbar,elementary_charge
    assert abs(p['time_unit_seconds']*p['voltage_unit_V']-hbar/(2*elementary_charge))<1e-29
