from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_matrix import Flow,R,mirror_field,prepare,cloud_target,advance,identity,sample
from simulate import initial


def test_whole_state_field_gradient_and_rhs_mirror():
    f=Flow(14,.02,workers=1); h=initial(f,1024); hm=mirror_field(f,h)
    x=np.random.default_rng(700).uniform(-2.9,2.9,(30,3))
    v,J=sample(f.real(h),x); vm,Jm=sample(f.real(hm),x@R)
    np.testing.assert_allclose(vm,v@R,atol=1e-13)
    np.testing.assert_allclose(Jm,R@J@R,atol=1e-13)
    np.testing.assert_allclose(f.rhs(hm)[0],mirror_field(f,f.rhs(h)[0]),atol=1e-11)


def test_anchored_geometry_keeps_centers_distances_and_labels():
    x=np.random.default_rng(701).uniform(-2,2,(30,3)); idx,r,local=prepare(x,8)
    np.testing.assert_allclose(local-x[:,None],r@R,atol=1e-14)
    np.testing.assert_allclose(np.linalg.norm(r,axis=-1),np.linalg.norm(r@R,axis=-1),atol=1e-14)
    perm=np.random.default_rng(702).permutation(len(x))
    _,rp,lp=prepare(x[perm],8)
    np.testing.assert_allclose(rp,r[perm],atol=1e-14)
    np.testing.assert_allclose(lp,local[perm],atol=1e-14)


def test_probe_addition_cannot_change_flow_or_center_targets():
    f=Flow(14,.02,workers=1); h=initial(f,1024)
    x=np.random.default_rng(703).uniform(-2,2,(30,3)); _,_,local=prepare(x,8)
    original=advance(f,h,x,identity(len(x)),.005)
    combined=np.concatenate([x,x@R,local.reshape(-1,3)])
    result=advance(f,h,combined,identity(len(combined)),.005)
    np.testing.assert_array_equal(result[0],original[0])
    np.testing.assert_array_equal(result[1][:len(x)],original[1])
    np.testing.assert_array_equal(result[2][:len(x)],original[2])


def test_finite_cloud_target_has_expected_affine_and_reflection_behavior():
    r=np.random.default_rng(704).normal(0,.1,(20,8,3))
    F=np.diag(np.exp([.4,-.2,-.2])); y=r@F.T
    np.testing.assert_allclose(cloud_target(r,y,.2),2,atol=1e-8)
    np.testing.assert_allclose(cloud_target(r@R,y@R,.2),cloud_target(r,y,.2),atol=1e-12)
    np.testing.assert_allclose(cloud_target(r@R,(r@R)@F.T,.2),2,atol=1e-8)


def test_four_cell_contrasts_separate_additive_and_interaction_cases():
    import importlib.util
    spec=importlib.util.spec_from_file_location('mirror_direction_analysis',Path(__file__).with_name('analyze.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cells={k:np.array([v,v]) for k,v in zip(['original_normal','original_reversed','mirror_normal','mirror_reversed'],[1.,3.,4.,6.])}
    result=module.factorial(cells)
    assert result['arrangement_average']['mean']==3
    assert result['direction_average']['mean']==2
    assert result['interaction_difference_in_differences']['RMS']==0
    cells['mirror_reversed']+=1
    assert module.factorial(cells)['interaction_difference_in_differences']['mean']==1
