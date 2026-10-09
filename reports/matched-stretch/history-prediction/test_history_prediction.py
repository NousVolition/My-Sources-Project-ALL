import json
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from simulate import Flow,L,initial,sample,ROOT
from features import present,history,future_target,neighbors,shuffle_blocks


def test_interpolant_gradient_matches_finite_differences():
    f=Flow(14,.02,workers=1)
    u=f.real(initial(f,19))
    x=np.random.default_rng(5).uniform(-2.9,2.9,(20,3))
    v,J=sample(u,x)
    for j in range(3):
        dx=np.eye(3)[j]*1e-6
        derivative=(sample(u,x+dx)[0]-sample(u,x-dx)[0])/2e-6
        np.testing.assert_allclose(J[:,:,j],derivative,rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(sample(u,x+L)[0],v,atol=1e-13)


def test_constant_velocity_and_uniform_affine_stretch():
    u=np.broadcast_to(np.array([1.,2.,-3.])[:,None,None,None],(3,8,8,8))
    v,J=sample(u,np.array([[2.99,0.,0.],[-3.01,0.,0.]]))
    np.testing.assert_allclose(v,[[1,2,-3],[1,2,-3]])
    np.testing.assert_allclose(J,0,atol=1e-14)
    t0,t1=.4,.7
    a=np.diag([2.,-1.,-1.])
    m0=np.diag(np.exp(np.diag(a)*t0))[None]
    m1=np.diag(np.exp(np.diag(a)*t1))[None]
    np.testing.assert_allclose(future_target(m0,m1,t1-t0),2)


def test_target_invariant_to_prior_tangent_basis():
    rng=np.random.default_rng(9)
    m0=np.eye(3)[None]+rng.normal(size=(10,3,3))*.1
    F=np.eye(3)[None]+rng.normal(size=(10,3,3))*.2
    m1=F@m0
    basis=np.eye(3)+rng.normal(size=(3,3))*.1
    np.testing.assert_allclose(future_target(m0,m1,.2),future_target(m0@basis,m1@basis,.2),atol=1e-12)


def fixture():
    rng=np.random.default_rng(8)
    x=rng.uniform(-2,2,(5,50,3))
    # Smooth, nondegenerate trajectories with no neighbor distance ties.
    x=x[:1]+np.arange(5)[:,None,None]*rng.normal(0,.03,(50,3))
    u=rng.normal(size=(50,3)); J=rng.normal(size=(50,3,3))
    return x,u,J


def test_particle_identity_permutation_equivariance():
    x,u,J=fixture()
    perm=np.random.default_rng(32).permutation(len(u))
    A=present(x[-1],u,J,8)[0]; H=history(x,u,J,8,.2)
    np.testing.assert_allclose(present(x[-1,perm],u[perm],J[perm],8)[0],A[perm],atol=1e-12)
    np.testing.assert_allclose(history(x[:,perm],u[perm],J[perm],8,.2),H[perm],atol=1e-10)


def test_static_geometry_history_zero():
    x,u,J=fixture(); x=np.repeat(x[:1],5,axis=0)
    np.testing.assert_allclose(history(x,u,J,8,.2),0,atol=1e-7)


def test_past_features_cannot_read_mutated_future():
    x,u,J=fixture(); full=np.concatenate([x,x+3])
    A=present(full[4],u,J,8)[0]; H=history(full[:5],u,J,8,.2)
    full[5:]=np.nan
    np.testing.assert_array_equal(A,present(full[4],u,J,8)[0])
    np.testing.assert_array_equal(H,history(full[:5],u,J,8,.2))


def test_independent_seed_splits_and_variants():
    p=json.loads((ROOT/"protocol.json").read_text())
    a,b,c=[set(p[k]) for k in ["train_seeds","validation_seeds","test_seeds"]]
    assert not a&b and not a&c and not b&c
    assert set(p["refinement_seeds"])<=c
    for t in p["anchors"]:
        assert t>=max(p["lookbacks"])
        assert t+p["gap"]+max(p["horizons"])<=p["base"]["end"]


def test_projected_solver_and_grid_consistent_initial_conditions():
    f=Flow(26,.02,workers=1); g=Flow(52,.02,workers=1)
    a=initial(f,1000); b=initial(g,1000)
    np.testing.assert_allclose(f.real(a),g.real(b)[:,::2,::2,::2],atol=1e-13)
    div=f.real(1j*sum(k*c for k,c in zip(f.k,a)))
    assert abs(div).max()<1e-13
    rh,_,_=f.rhs(a)
    assert f.inner(a,rh)<0


def test_shuffle_stays_in_run_and_anchor():
    values=np.arange(24); group=np.repeat([1,2],12); anchor=np.tile(np.repeat([.1,.2],6),2)
    shuffled=shuffle_blocks(values,group,anchor,17)
    assert not np.array_equal(values,shuffled)
    for g in [1,2]:
        for t in [.1,.2]:
            ix=(group==g)&(anchor==t)
            assert sorted(values[ix])==sorted(shuffled[ix])


def test_heun_navier_stokes_against_decaying_taylor_green():
    f=Flow(14,.1,workers=1)
    axis=np.arange(f.n)*f.dx
    x,y,z=np.meshgrid(axis,axis,axis,indexing="ij")
    k=2*np.pi/L
    u=np.stack([np.sin(k*x)*np.cos(k*y),-np.cos(k*x)*np.sin(k*y),np.zeros_like(x)])
    h0=f.project(f.hat(u))
    errors=[]
    for dt in [.04,.02]:
        h=h0.copy()
        for _ in range(round(.4/dt)): h,_,_=f.step(h,dt)
        exact=h0*np.exp(-2*.1*k*k*.4)
        errors.append(np.linalg.norm(h-exact)/np.linalg.norm(exact))
    assert errors[0]/errors[1]>3.8
    assert errors[-1]<1e-6


def test_rigid_rotation_has_zero_target_and_affine_history_recovers_stretch():
    x,u,J=fixture()
    r=np.array([[0,-1.,0],[1.,0,0],[0,0,1.]])
    np.testing.assert_allclose(future_target(np.eye(3)[None],r[None],.2),0,atol=1e-12)
    # Keep this manufactured affine test far from periodic boundaries.
    x=x[:1]*.1
    rates=np.array([.2,-.1,-.1])
    path=x*np.exp(np.linspace(0,.2,5)[:,None,None]*rates)
    H=history(path,np.zeros_like(u),np.zeros_like(J),8,.2)
    np.testing.assert_allclose(H[:,15:18],np.tile([.2,-.1,-.1],(len(u),1)),atol=1e-7)


def test_model_preprocessing_and_selection_do_not_fit_test_data():
    from evaluate import fit_pair
    rng=np.random.default_rng(7)
    X=rng.normal(size=(90,4)); y=X[:,0]+rng.normal(size=90)*.1
    masks=[np.arange(90)//30==i for i in range(3)]
    group=np.arange(90)//15
    event=(y>=np.median(y[masks[0]])).astype(int)
    _,model=fit_pair(X,y,event,masks,group)
    X2=X.copy();X2[masks[2]]+=10000
    y2=y.copy();y2[masks[2]]=-999
    event2=event.copy();event2[masks[2]]=1-event2[masks[2]]
    _,model2=fit_pair(X2,y2,event2,masks,group)
    for m,n in zip(model,model2):
        np.testing.assert_allclose(m.steps[0][1].mean_,X[masks[0]].mean(0))
        np.testing.assert_array_equal(m.steps[0][1].mean_,n.steps[0][1].mean_)
        np.testing.assert_array_equal(m.steps[-1][1].coef_,n.steps[-1][1].coef_)


def test_saved_forecaster_rejects_wrong_window():
    from predict import forecast
    x,u,J=fixture()
    observation={"time":np.linspace(0,.3,5),"positions":x,"velocity":u,"gradient":J}
    with pytest.raises(ValueError,match="lookback"):
        forecast(observation,{"lookback":.2})
