import numpy as np
import pytest
from core import PilotFlow,initial,energy,L
from recovery import rk4,strength
from network import matrix,integrate,rhs_log,switches,fit_coefficients,waiting
from common import arrivals,P

def test_residence_requires_full_observed_window():
    t=np.arange(6.)
    assert arrivals(t,[2,2,2,2,0,0],1,2)==(4.,None)
    assert arrivals(t,[2,0,2,0,0,0],1,2)==(1.,3.)

def test_forcing_controls():
    assert strength(0,'abrupt',2)==0
    assert strength(1,'gradual',2)==.5
    assert strength(3,'gradual',2)==0
    assert strength(100,'unchanged',2)==1

def test_beltrami_decay_and_energy_budget():
    f=PilotFlow(12,.2);h=initial(f,7,kind='beltrami');h0=h.copy();e0=energy(f,h);budget=np.zeros(2)
    for j in range(50):h,b,_=rk4(f,h,j*.02,.02,np.zeros_like(h),'abrupt',2);budget+=b
    expected=h0*np.exp(-.2*(2*np.pi/L)**2)
    assert np.sqrt(f.inner(h-expected,h-expected)/f.inner(expected,expected))<1e-10
    assert abs(energy(f,h)-e0+budget[1])/e0<1e-9

def test_forced_beltrami_exact_solution():
    f=PilotFlow(12,.2);shape=initial(f,9,kind='beltrami');h=.4*shape;rate=.2*(2*np.pi/L)**2
    for j in range(50):h,_,_=rk4(f,h,j*.02,.02,.06*shape,'unchanged',2)
    expected=(.4*np.exp(-rate)+.06/rate*(1-np.exp(-rate)))*shape
    assert np.sqrt(f.inner(h-expected,h-expected)/f.inner(expected,expected))<1e-10

def test_independent_clones_and_divergence():
    f=PilotFlow(12,.2);h=initial(f,2);a=h.copy();b=h.copy()
    a,_,_=rk4(f,a,0,.01,np.zeros_like(a),'abrupt',2)
    assert np.array_equal(b,h)
    assert np.max(abs(f.real(1j*sum(k*x for k,x in zip(f.k,a)))))<1e-12

def test_network_saddles_and_invariant_planes():
    C=matrix()
    for j in range(3):
        x=np.eye(3)[j];f=x*(1-C@x)
        J=np.diag(1-C@x)-np.diag(x)@C
        assert np.max(abs(f))==0
        assert np.allclose(np.sort(np.linalg.eigvals(J)),[-1,-.6,.4])
    x=np.array([.5,.7,0]);assert (x*(1-C@x))[2]==0

def test_uncoupled_logistic_solution():
    x0=np.array([.1,.2,.7]);t,y=integrate(np.log(x0),matrix(0),5)
    exact=1/(1+(1/x0-1)*np.exp(-t[:,None]))
    assert np.max(abs(np.exp(y)-exact))<1e-8

def test_cycle_order_and_refinement():
    t,y=integrate(np.log([.6,.2,.1]),matrix(),100)
    labels,ix=switches(t,y)
    assert len(ix)>=5
    assert np.all((labels[ix]-labels[ix-1])%3==1)
    _,fine=integrate(np.log([.6,.2,.1]),matrix(),100,tight=True)
    assert np.max(abs(np.exp(y)-np.exp(fine)))<1e-7

def test_boundary_connection_and_coupling_regimes():
    from scipy.integrate import solve_ivp
    C=matrix();x0=np.array([1-1e-6,1e-6,0.])
    s=solve_ivp(lambda t,x:x*(1-C@x),(0,110),x0,rtol=1e-10,atol=1e-12)
    assert np.max(abs(s.y[:,-1]-[0,1,0]))<1e-7
    for coupling,stable in ((.8,True),(1.,False)):
        C=matrix(coupling);coexist=np.linalg.solve(C,np.ones(3));ev=np.linalg.eigvals(-np.diag(coexist)@C)
        assert bool(np.max(ev.real)<0)==stable

def test_known_parameter_recovery_from_observations():
    data=[]
    for x in ([.6,.2,.1],[.15,.75,.08],[.5,.1,.9]):
        _,y=integrate(np.log(x),matrix(),30);data.append({'base':y})
    estimate=fit_coefficients(data)
    assert np.max(abs(estimate-matrix()))<.003

def test_censored_wait_is_restricted_not_dropped():
    t=np.arange(0,11.);y=np.c_[np.ones(11),np.zeros(11),np.zeros(11)]
    assert waiting(t,y,0,5)==(5,True)
    y[4:,1]=2
    assert waiting(t,y,0,5)==(4.,False)

def test_coexistence_is_not_a_saddle_visit():
    from network import saddle_visits
    assert len(saddle_visits(np.log(np.full((10,3),.35))))==0
    assert saddle_visits(np.log([[.95,.02,.02],[.94,.02,.02],[.02,.95,.02]])).tolist()==[0,1]

def test_run_splits_do_not_overlap():
    for system in ('fluid','network'):
        sets=[set(P[system][s+'_seeds']) for s in ('train','validation','test')]
        assert all(not (sets[i]&sets[j]) for i,j in ((0,1),(0,2),(1,2)))

def test_recurrence_requires_exit_then_sustained_return():
    from recurrence import return_events
    t=np.arange(6.)
    assert return_events(t,np.zeros(6),residence=1)==[]
    assert return_events(t,np.array([0,.3,.05,.15,.05,.05]),residence=1)==[4.]
