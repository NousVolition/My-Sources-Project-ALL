import json
import numpy as np
from core import PilotFlow,initial,impulse,energy,advance,advance_pair,sample,wrap,L


def test_projection_and_pressure():
    f=PilotFlow(26,.03)
    r=np.random.default_rng(12).normal(size=(3,26,26,26))
    h=f.project(f.hat(r))
    assert np.max(abs(sum(k*c for k,c in zip(f.k,h))))<1e-10
    assert np.allclose(f.project(h),h,atol=1e-11)


def test_beltrami_exact_decay():
    f=PilotFlow(26,.07);h=initial(f,11,kind="beltrami");h0=h.copy()
    x=np.array([[.13,.27,.31]])
    for _ in range(10):h,x,*_=advance(f,h,x,.02)
    expected=h0*np.exp(-.07*(2*np.pi/L)**2*.2)
    assert np.linalg.norm(h-expected)/np.linalg.norm(expected)<1e-10


def test_nonlin_energy_and_enstrophy_balance():
    f=PilotFlow(26,.03);h=initial(f,32);r=f.rhs_only(h)
    w=f.curl(h)
    assert abs(f.inner(h,r)+f.nu*f.inner(w,w))<1e-10
    wr=f.real(w);G=f.gradient(h)
    prod=np.mean(np.einsum('iabc,jabc,ijabc->abc',wr,wr,G))*L**3
    assert abs(f.inner(w,f.curl(r))-prod+f.nu*f.inner(w,f.k2*w))<1e-9


def test_difference_energy_production_identity():
    f=PilotFlow(26,.03);b=initial(f,5);p=b+impulse(f,.05);d=p-b
    u=f.real(d);G=f.gradient(b)
    prod=-np.mean(np.einsum('iabc,jabc,ijabc->abc',u,u,G))
    loss=f.nu*f.inner(f.curl(d),f.curl(d))/L**3
    rate=f.inner(d,f.rhs_only(p)-f.rhs_only(b))/L**3
    assert abs(rate-prod+loss)<1e-13


def test_identical_clones_and_marker_label_independence():
    f=PilotFlow(26,.03);h=initial(f,8);x=np.random.default_rng(5).uniform(-3,3,(12,3))
    b,p,xb,xp,*_=advance_pair(f,h,h.copy(),x,x.copy(),.01)
    assert np.array_equal(b,p) and np.array_equal(xb,xp)
    _,xx,*_=advance(f,h,x[::-1],.01)
    assert np.allclose(xx[::-1],xb,atol=1e-15)


def test_grid_independent_start_and_impulse():
    f=PilotFlow(26,.03);g=PilotFlow(38,.03)
    for h,k in [(initial(f,5),initial(g,5)),(impulse(f,.05),impulse(g,.05))]:
        for i in range(-5,6):
            for j in range(-5,6):
                assert np.allclose(h[:,i%26,j%26,:6]/26**3,k[:,i%38,j%38,:6]/38**3,atol=2e-15)


def test_impulse_strength_and_localization():
    f=PilotFlow(26,.03);h=impulse(f,.05)
    assert np.isclose(energy(f,h),.5*.05**2,rtol=1e-13)
    a=np.arange(26)*L/26-L/2
    rr=sum(x*x for x in np.meshgrid(a,a,a,indexing='ij',sparse=True))
    e=np.sum(f.real(h)**2,axis=0)
    assert e[rr<1.5**2].sum()/e.sum()>.8


def test_constant_velocity_marker_and_periodic_wrap():
    f=PilotFlow(26,.03);u=np.zeros((3,26,26,26));u[0]=1
    h=f.hat(u);x=np.array([[2.99,0.,0.]])
    _,xx,*_=advance(f,h,x,.1)
    assert np.allclose(xx,[[3.09,0,0]],atol=1e-14)
    assert np.allclose(wrap(xx),[[-2.91,0,0]])


def test_spline_refinement():
    x=np.random.default_rng(1).uniform(-3,3,(50,3));errors=[]
    for n in [26,38,50]:
        f=PilotFlow(n,.03);a=np.arange(n)*f.dx-L/2
        u=np.broadcast_to(np.sin(2*np.pi*a[:,None,None]/L),(n,n,n))
        errors.append(np.max(abs(sample(u,x)-np.sin(2*np.pi*x[:,0]/L))))
    assert errors[1]<errors[0] and errors[2]<errors[1]
    assert errors[-1]<1e-6


def test_stokes_null_exact_modes():
    f=PilotFlow(26,.03,nonlinear=False);h=initial(f,8);h0=h.copy();x=np.zeros((1,3))
    for _ in range(10):h,x,*_=advance(f,h,x,.01)
    assert np.linalg.norm(h-h0*np.exp(-f.nu*f.k2*.1))/np.linalg.norm(h)<1e-10


def test_seed_groups_disjoint():
    from core import ROOT
    p=json.loads((ROOT/'protocol.json').read_text())
    groups=[set(p[k]) for k in ['train_seeds','validation_seeds','test_seeds','regime_seeds']]
    assert sum(map(len,groups))==len(set.union(*groups))


def test_particle_drag_exact_and_step_refinement():
    from particles import particle_step,coefficients
    u=np.zeros((3,12,12,12));u[0]=1.;coeff=coefficients(u)
    errors=[]
    for ratio in [.5,.25]:
        x,v=particle_step(np.zeros((1,3)),np.zeros((1,3)),np.array([.02]),np.zeros((1,3)),coeff,coeff,.1,ratio)
        exactx=.1-.02*(1-np.exp(-5));exactv=1-np.exp(-5)
        errors.append(abs(v[0,0]-exactv)+abs(x[0,0]-exactx))
    assert errors[1]<errors[0]/10 and errors[1]<2e-6


def test_history_features_relabel_and_past_only():
    from prediction import features
    rng=np.random.default_rng(82);x=rng.uniform(-2,2,(32,3));ts=np.array([0,.1,.2,.3,.4])
    pos=np.array([x*(1+t*.02) for t in ts]);v=rng.normal(size=x.shape);G=rng.normal(size=(32,3,3));z=np.zeros_like(v);Z=np.zeros_like(G)
    a,b,nn,*_=features(ts,pos,v,z,G,Z)
    perm=rng.permutation(32);inv=np.argsort(perm)
    aa,bb,*_=features(ts,pos[:,perm],v[perm],z[perm],G[perm],Z[perm])
    assert np.allclose(a,aa[inv]) and np.allclose(b,bb[inv])
    assert np.isfinite(b).all()


def test_phase_locking_fixed_point_and_uncoupled():
    from phase_locking import evolve
    ts=np.linspace(0,20,101)
    x=evolve(.5,1.,[np.arcsin(.5)],ts)
    assert np.max(abs(x-np.arcsin(.5)))<1e-12
    y=evolve(.7,0.,[.3],ts)
    assert np.allclose(y[:,0],.3+.7*ts,atol=1e-10)


def test_triangle_lock_and_equal_peak():
    from dynamics import triangle,tri_evolve
    a=np.linspace(-10,10,200)
    assert np.allclose(triangle(a),triangle(a+2*np.pi))
    x=tri_evolve(.5,[.5,.55],np.linspace(0,2,21))
    assert np.max(abs(x[:,0]-.5))<1e-12
    assert np.isclose(x[-1,1]-.5,.05*np.exp(-2),rtol=1e-8)


def test_nonnormal_gain_and_trace_classification():
    from scipy.linalg import expm
    from dynamics import classify
    A=np.array([[-1.,8.],[0.,-2.]])
    assert classify(A)=='stable node'
    assert np.max(np.linalg.svd(expm(A*.6),compute_uv=False))**2>4
    assert classify(np.diag([-1.,1.]))=='saddle'


def test_conditional_circuit_decoupling_and_symmetry():
    from dynamics import circuit_rhs
    x=np.array([.2,.8]);bias=1.2
    assert np.allclose(circuit_rhs(x,bias,0),bias-np.sin(x))
    assert np.allclose(circuit_rhs(x,bias,.5)[::-1],circuit_rhs(x[::-1],bias,.5))
    y=circuit_rhs(np.array([.3,.3]),bias,.5)
    assert np.allclose(y,(bias-np.sin(.3))/2)


def test_pendulum_period_and_conservation():
    from scipy.integrate import solve_ivp
    from pendulum import period,rhs,hamiltonian
    assert np.isclose(period(0),2*np.pi)
    assert period(2.8)>period(1)>period(.1)
    T=period(2.);sol=solve_ivp(rhs,(0,T),[2.,0.],rtol=1e-11,atol=1e-13,method='DOP853')
    assert np.linalg.norm(sol.y[:,-1]-[2.,0.])<1e-9
    H=hamiltonian(sol.y.T);assert np.max(abs(H-H[0]))<1e-9


def test_relaxation_energy_identity_and_origin():
    from relaxation import rhs,cubic
    assert np.allclose(rhs(0,[0.,0.],10),[0.,0.])
    x,y,mu=.6,.7,5.
    xd,yd=rhs(0,[x,y],mu);xdd=mu*(yd-(x*x-1)*xd)
    assert np.isclose(x*xd+xd*xdd,mu*(1-x*x)*xd*xd)
    assert np.isclose(cubic(1),-2/3) and np.isclose(cubic(-1),2/3)


def test_duffing_period_and_energy_identity():
    from weak_nonlinear import duffing,duffing_period
    assert np.isclose(duffing_period(1,0),2*np.pi)
    assert duffing_period(2,.1)<duffing_period(1,.1)<2*np.pi
    x,v,eps=.4,.8,.1;dx,dv=duffing(0,[x,v],eps)
    assert abs(v*dv+(x+eps*x**3)*dx)<1e-14


def test_parametric_exact_rest_and_instability():
    from parametric import swing,monodromy
    assert np.array_equal(swing(.1,0,0,np.linspace(0,10,11)),np.zeros((11,3)))
    M,eig=monodromy(.1,0);assert abs(np.linalg.det(M)-1)<1e-10
    assert max(abs(eig))>1.05
    _,stable=monodromy(.1,.75);assert max(abs(stable))<1+1e-9


def test_bifurcation_radial_identity_and_bottleneck():
    from bifurcations import bottleneck_time,hopf
    from scipy.integrate import quad
    eps=.04;numeric=quad(lambda x:1/(eps+x*x),-1,1)[0]
    assert np.isclose(numeric,bottleneck_time(eps),rtol=1e-12)
    x,y,mu=.3,.4,.2;dx,dy=hopf(0,[x,y],mu);r2=x*x+y*y
    assert np.isclose(x*dx+y*dy,mu*r2-r2*r2)
