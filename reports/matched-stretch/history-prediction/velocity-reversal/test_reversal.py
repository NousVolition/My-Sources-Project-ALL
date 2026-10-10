from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_pairs import advance,directional,identity,Flow,initial,sample


def test_velocity_reversal_preserves_energy_and_changes_gradient_sign():
    f=Flow(14,.02,workers=1); h=initial(f,1024)
    x=np.random.default_rng(602).uniform(-2.9,2.9,(30,3))
    v,J=sample(f.real(h),x); vm,Jm=sample(f.real(-h),x)
    np.testing.assert_allclose(vm,-v,atol=1e-14)
    np.testing.assert_allclose(Jm,-J,atol=1e-14)
    assert f.inner(h,h)==f.inner(-h,-h)
    # Reversing u preserves the nonlinear RHS and reverses the viscous term.
    np.testing.assert_allclose(f.rhs(-h)[0]-f.rhs(h)[0],2*f.nu*f.k2*h,atol=1e-11)


def test_joint_step_matches_original_fluid_step_and_uniform_advection():
    f=Flow(14,.02,workers=1); h=initial(f,1024)
    x=np.random.default_rng(600).uniform(-2,2,(20,3))
    result=advance(f,h,x,identity(len(x)),.005)
    np.testing.assert_array_equal(result[0],f.step(h,.005)[0])
    u=np.broadcast_to(np.array([.3,-.2,.1])[:,None,None,None],(3,14,14,14)).copy()
    for sign in [-1,1]:
        _,xx,M,_,_=advance(f,f.hat(sign*u),x,identity(len(x)),.01)
        np.testing.assert_allclose(xx,x+sign*.01*np.array([.3,-.2,.1]),atol=1e-14)
        np.testing.assert_allclose(M,identity(len(x)),atol=1e-14)


def test_directional_rates_match_affine_separation_and_ignore_identity():
    rng=np.random.default_rng(603); x=rng.uniform(-.4,.4,(30,3))
    J=np.broadcast_to(np.diag([2.,-1.,-1.]),(len(x),3,3)); v=x@J[0].T
    rates=directional(x,v,J,8)
    np.testing.assert_allclose(rates[:,0],rates[:,1],atol=1e-13)
    np.testing.assert_allclose(directional(x,-v,-J,8),-rates,atol=1e-13)
    perm=rng.permutation(len(x))
    np.testing.assert_allclose(directional(x[perm],v[perm],J[perm],8),rates[perm],atol=1e-13)


def test_audit_accepts_only_newline_equivalent_source(tmp_path):
    import hashlib
    import importlib.util
    spec=importlib.util.spec_from_file_location('velocity_analysis',Path(__file__).with_name('analyze.py'))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    p=tmp_path/'code.py'; p.write_bytes(b'x=1\n')
    assert module.same_source(p,hashlib.sha256(b'x=1\r\n').hexdigest())
    assert not module.same_source(p,hashlib.sha256(b'x=2\n').hexdigest())
    comparison=module.compare(np.array([0.,1.]),np.array([1.,0.]),.5)
    assert comparison['same_label_RMSE']==1 and comparison['event_Jaccard']==0
