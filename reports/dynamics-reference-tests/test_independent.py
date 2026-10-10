"""Independent checks in addition to the study's recorded acceptance checks."""
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.linalg import expm
from core import integrate, lorenz, jacobian, delayed_branch, wheel_rhs, wheel_project


@pytest.mark.parametrize("method,expected", [("euler",1),("heun",2),("rk4",4)])
def test_vector_problem_against_matrix_exponential(method, expected):
    matrix = np.array([[-1., .3],[-.7,-2.]])
    initial = np.array([.4,-.2])
    exact = expm(matrix) @ initial
    errors=[]
    for step in (.04,.02,.01):
        _,state=integrate(lambda y:matrix@y,initial,1.,step,method)
        errors.append(np.linalg.norm(state[-1]-exact))
    assert abs(np.log2(errors[-2]/errors[-1])-expected)<.12


def test_jacobian_by_finite_differences():
    state=np.array([3.,-5.,14.])
    step=1e-5
    finite=np.column_stack([(lorenz(state+step*d)-lorenz(state-step*d))/(2*step) for d in np.eye(3)])
    np.testing.assert_allclose(jacobian(state),finite,rtol=1e-9,atol=1e-8)


def test_logistic_reciprocal_change_of_variable():
    # For x'=x(1-x), u=1/x satisfies u'=1-u while x is nonzero.
    t,x=integrate(lambda y:y*(1-y),[.2],4.,.005)
    _,u=integrate(lambda y:1-y,[5.],4.,.005)
    np.testing.assert_allclose(u[:,0],1/x[:,0],rtol=0,atol=2e-11)


@pytest.mark.parametrize("sign", [-1,1])
def test_delayed_branch_derivative_independently(sign):
    t=np.array([.1,.5,1.1,1.8,3.])
    h=1e-5
    right=delayed_branch(t+h,1.,sign)[0]
    left=delayed_branch(t-h,1.,sign)[0]
    x=delayed_branch(t,1.,sign)[0]
    np.testing.assert_allclose((right-left)/(2*h),np.cbrt(x),atol=2e-9,rtol=1e-8)
    # At the joining point the derivative approaches zero from both sides.
    for h in (1e-4,1e-6,1e-8):
        assert abs(delayed_branch(np.array([1.+h]),1.,sign)[0][0]/h) < np.sqrt(h)


def test_waterwheel_instantaneous_mass_budget_for_nonuniform_state():
    theta,rhs=wheel_rhs(64)
    density=30+np.random.default_rng(12).uniform(-2,2,64)
    state=np.r_[density,3.5]
    rate=rhs(state)
    np.testing.assert_allclose(2*np.pi*np.mean(rate[:-1]),2*np.pi*60-2*np.pi*np.mean(density),atol=1e-11)


def test_periodic_transport_against_translated_profile():
    theta,rhs=wheel_rhs(32,q0=0.,q1=0.,leak=0.,drag=0.,torque=0.)
    profile=lambda a:1+.1*np.cos(3*a)+.2*np.sin(2*a)
    _,state=integrate(rhs,np.r_[profile(theta),1.7],.8,.0005)
    np.testing.assert_allclose(state[-1,:-1],profile(theta-1.7*.8),rtol=0,atol=2e-12)


def test_no_torque_ablation_has_exact_drag_decay():
    theta,rhs=wheel_rhs(16,torque=0.)
    t,state=integrate(rhs,np.r_[np.full(16,60.),1.],1.,.001)
    np.testing.assert_allclose(state[:,-1],np.exp(-10*t),atol=4e-11,rtol=0)


def test_symmetric_water_no_initial_rotation_stays_still():
    theta,rhs=wheel_rhs(16,q1=0.)
    t,state=integrate(rhs,np.r_[np.full(16,45.),0.],1.,.001)
    assert np.max(abs(state[:,-1])) < 1e-12
    np.testing.assert_allclose(np.mean(state[:,:-1],axis=1),60-15*np.exp(-t),atol=1e-11,rtol=0)


def test_higher_waterwheel_mode_has_no_torque_feedback():
    theta,rhs=wheel_rhs(32)
    base=45+np.sin(theta)+27*np.cos(theta)
    _,a=integrate(rhs,np.r_[base,1.],1.,.001)
    _,b=integrate(rhs,np.r_[base+2*np.cos(3*theta),1.],1.,.001)
    np.testing.assert_allclose(wheel_project(a,theta),wheel_project(b,theta),rtol=0,atol=1e-10)


def test_frozen_study_outputs_pass_all_acceptance_checks():
    checks=json.loads((Path(__file__).parent/"data"/"checks.json").read_text())
    assert len(checks)>=50
    assert all(item["passed"] for item in checks)
