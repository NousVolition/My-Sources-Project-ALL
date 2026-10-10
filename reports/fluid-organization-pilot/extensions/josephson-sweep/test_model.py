import numpy as np
import pytest
from scipy.integrate import solve_ivp
from model import rhs, integrate, voltage, potential, dissipation, power_residual, synchronized_voltage
from model import to_state, from_state, stable_rhs, stable_velocity, integrate_state


def test_transformed_equations_match_original():
    p=np.array([[.3,-1.2],[2.1,.1],[-.2,-1.],[-2.,2.]])
    s,sign=to_state(p)
    np.testing.assert_allclose(from_state(s,sign),p,atol=1e-15)
    np.testing.assert_allclose(stable_velocity(s,sign,1.005,2.),rhs(p,1.005,2.),atol=1e-15)


def test_transformed_short_trajectory_matches_direct_reference():
    p0=np.array([1.28,2.09]);s,sign=to_state(p0)
    _,states,_=integrate_state(s,sign,1.005,2.,30,.025,stride=1.)
    ref=solve_ivp(lambda t,y:rhs(y,1.005,2.),(0,30),p0,method='DOP853',rtol=1e-12,atol=1e-13,t_eval=np.arange(31),max_step=.1)
    np.testing.assert_allclose(from_state(states,sign),ref.y.T,atol=1e-8,rtol=0)


def test_tiny_separation_is_retained_and_can_regrow():
    state=np.array([np.pi,-100.])
    # Physical output rounds to equal phases, but internal separation survives.
    assert from_state(state,1.)[0]==from_state(state,1.)[1]
    assert stable_rhs(state,0.,2.)[1]>0
    _,s,_=integrate_state(state,1.,0.,2.,1,.025,stride=1.)
    assert s[-1,1]>state[1] and np.isfinite(s).all()


def test_transformed_synchronized_branch():
    s,sign=to_state(np.zeros(2))
    _,p,_=integrate_state(s,sign,1.2,.5,10,.025,stride=1.)
    physical=from_state(p,sign)
    np.testing.assert_array_equal(physical[...,0],physical[...,1])
    _,direct,_=integrate(np.zeros(2),1.2,.5,10,.025)
    np.testing.assert_allclose(physical,direct,atol=1e-12)


def test_kirchhoff_current_balance():
    p=np.array([[.2,2.1],[-.7,.4],[7.2,-2.3]])
    a=np.array([0,.5,2]); b=np.array([.9,1.01,2])
    v=rhs(p,b,a)
    np.testing.assert_allclose(v+a[:,None]*v.sum(axis=-1)[:,None]+np.sin(p),np.repeat(b[:,None],2,axis=1),atol=2e-15)


def test_passivity_identity():
    p=np.array([[.6,-1.2],[2,5.]])
    v=rhs(p,1.1,.5)
    np.testing.assert_allclose(1.1*v.sum(-1)-(np.sin(p)*v).sum(-1),dissipation(v,.5),atol=1e-14)
    assert np.all(dissipation(v,.5)>=0)


def test_static_solution_below_threshold():
    p=np.full((4,2),np.arcsin(.8))
    _,q,_=integrate(p,.8,.5,10,.05)
    np.testing.assert_allclose(q[-1],p,atol=1e-13)


def test_unloaded_independent_junctions():
    p=np.array([.1,2.2]); v=rhs(p,1.2,0)
    np.testing.assert_allclose(v,1.2-np.sin(p))
    assert rhs([.1,-8.],1.2,0)[0]==v[0]


def test_exchange_and_periodicity():
    p=np.array([.2,2.])
    np.testing.assert_allclose(rhs(p[::-1],1.3,.5),rhs(p,1.3,.5)[::-1])
    np.testing.assert_allclose(rhs(p+2*np.pi*np.array([2,-3]),1.3,.5),rhs(p,1.3,.5),atol=1e-14)


def test_unwrapped_voltage():
    np.testing.assert_allclose(voltage([0,0],[10*np.pi,8*np.pi],20),[np.pi/2,2*np.pi/5])


def test_rk4_refinement_and_reference():
    p0=np.array([.1,2.]); end=10
    ref=solve_ivp(lambda t,p:rhs(p,1.2,.5),(0,end),p0,method='DOP853',rtol=1e-12,atol=1e-13,max_step=.1).y[:,-1]
    errs=[]
    for dt in [.2,.1,.05]:
        _,p,_=integrate(p0,1.2,.5,end,dt)
        errs.append(np.linalg.norm(p[-1]-ref))
    assert errs[0]/errs[1]>12 and errs[1]/errs[2]>12
    assert errs[-1]<1e-6


def test_integrated_power_budget():
    _,p,d=integrate(np.array([[.2,2.],[1.,-.3]]),1.1,2.,40,.05)
    assert np.max(abs(power_residual(p,d,1.1)))<2e-6
    assert np.all(np.diff(d,axis=0)>=0)


def test_synchronized_full_cycle():
    b=1.1; a=.5; period=2*np.pi*(1+2*a)/np.sqrt(b*b-1)
    # Choose dt to end at exactly one analytical cycle.
    _,p,_=integrate(np.zeros(2),b,a,period,period/4000,stride=period/100)
    np.testing.assert_allclose(p[-1],2*np.pi,atol=1e-8)
    np.testing.assert_allclose(voltage(p[0],p[-1],period).sum(),synchronized_voltage(b,a),rtol=1e-8)


def test_full_state_replay_has_no_hidden_history():
    _,p,_=integrate(np.array([.7,-.9]),1.2,.5,10,.05)
    _,a,_=integrate(p[-1],.99,.5,10,.05)
    _,b,_=integrate(p[-1].copy(),.99,.5,10,.05)
    np.testing.assert_array_equal(a,b)


@pytest.mark.parametrize('dt,stride,duration',[(0,1,10),(.3,1,10),(.1,.25,10),(.1,1,10.1)])
def test_reject_invalid_time_grid(dt,stride,duration):
    with pytest.raises(ValueError):integrate(np.zeros(2),1,.5,duration,dt,stride)
