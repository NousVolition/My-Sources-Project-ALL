"""Independent analytic checks for the new observation method and statistics."""
import importlib.util
from pathlib import Path
import numpy as np
import sys,json

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
spec=importlib.util.spec_from_file_location('fresh_sampler',HERE/'sampler.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
sample=module.sample


def field(n):
    axis=np.arange(n)*6/n-3
    x,y,z=np.meshgrid(axis,axis,axis,indexing='ij')
    k=2*np.pi/6
    return np.stack([np.sin(k*y),np.sin(k*z),np.sin(k*x)])


def test_exact_interpolant_jacobian_and_periodicity():
    x=np.random.default_rng(61).uniform(-2.9,2.9,(80,3));u=field(22)
    velocity,J=sample(u,x)
    for j in range(3):
        dx=np.zeros_like(x);dx[:,j]=1e-6
        derivative=(sample(u,x+dx)[0]-sample(u,x-dx)[0])/2e-6
        np.testing.assert_allclose(J[:,:,j],derivative,atol=2e-9,rtol=2e-8)
    vp,Jp=sample(u,x+np.array([12,-18,6]))
    np.testing.assert_allclose(vp,velocity,atol=2e-14,rtol=0)
    np.testing.assert_allclose(Jp,J,atol=2e-14,rtol=0)


def test_analytic_velocity_and_gradient_refinement():
    x=np.random.default_rng(62).uniform(-3,3,(500,3));k=2*np.pi/6
    exact=np.column_stack([np.sin(k*x[:,1]),np.sin(k*x[:,2]),np.sin(k*x[:,0])])
    J=np.zeros((len(x),3,3))
    for i,j in [(0,1),(1,2),(2,0)]:J[:,i,j]=k*np.cos(k*x[:,j])
    errors=[]
    for n in [14,22,46]:
        v,g=sample(field(n),x)
        errors.append([np.sqrt(np.mean((v-exact)**2)),np.sqrt(np.mean((g-J)**2))])
    errors=np.array(errors)
    assert np.all(errors[1:]<errors[:-1]/2.5)
    assert errors[-1,0]<2e-5 and errors[-1,1]<1e-4


def test_consistent_reflection_and_particle_relabeling():
    n=22;u=field(n);x=np.random.default_rng(63).uniform(-3,3,(80,3))
    R=np.diag([-1,1,1]);ur=np.take(u,(-np.arange(n))%n,axis=1).copy();ur[0]*=-1
    v,J=sample(u,x);vr,Jr=sample(ur,x@R)
    np.testing.assert_allclose(vr,v@R,atol=1e-14,rtol=0)
    np.testing.assert_allclose(Jr,R@J@R,atol=1e-14,rtol=0)
    perm=np.random.default_rng(64).permutation(len(x));vp,Jp=sample(u,x[perm])
    np.testing.assert_array_equal(vp,v[perm]);np.testing.assert_array_equal(Jp,J[perm])


def fixture_archive(path,permutation=None,future_change=False):
    time=np.arange(49)*.025;rng=np.random.default_rng(69);x=rng.uniform(-2.5,2.5,(32,3))
    positions=np.array([x+np.sin(t)*.03*x[:,[1,2,0]] for t in time])
    velocity=np.array([np.cos(t)*.03*x[:,[1,2,0]] for t in time])
    gradient=np.zeros((49,32,3,3));gradient[:,:,0,0]=.1;gradient[:,:,1,1]=-.07;gradient[:,:,2,2]=-.03
    tangent=np.array([np.broadcast_to(np.diag(np.exp([.1*t,-.07*t,-.03*t])),(32,3,3)) for t in time])
    after=time>.3+1e-10
    if future_change:positions[after]+=rng.normal(size=positions[after].shape);velocity[after]*=17;gradient[after]+=100
    if permutation is not None:
        positions=positions[:,permutation];velocity=velocity[:,permutation];gradient=gradient[:,permutation];tangent=tangent[:,permutation]
    np.savez_compressed(path,time=time,positions=positions,velocity=velocity,gradient=gradient,tangent=tangent)
    path.with_suffix('.json').write_text(json.dumps({'seed':69}))


def test_features_are_causal_and_follow_arrangement_under_relabeling(tmp_path):
    from fresh_features import make_dataset
    p={'anchors':[.3],'gap':.025,'primary_lookback':.2,'primary_horizon':.2,'neighbors':8,'base':{'markers':32}}
    paths=[tmp_path/(n+'.npz') for n in ['base','future','permuted']];perm=np.random.default_rng(70).permutation(32)
    fixture_archive(paths[0]);fixture_archive(paths[1],future_change=True);fixture_archive(paths[2],permutation=perm)
    a,b,c=[make_dataset([path],p) for path in paths]
    for key in ['present','history','reverse','irrelevant','mirrored','mirrored_reverse','target']:
        np.testing.assert_array_equal(a[key],b[key])
        np.testing.assert_allclose(a[key][perm],c[key],atol=2e-8,rtol=2e-8)


def test_spatial_metric_compares_places_at_same_time_and_bootstraps_64_runs():
    from fresh_statistics import spatial_scores,paired_values,classify
    event=np.array([1,0,0,0,0,1]);prob=np.array([.9,.8,.7,.3,.2,.1])
    group=np.zeros(6);anchor=np.repeat([.3,.4],3)
    r=spatial_scores(event,prob,group,anchor)
    assert np.isclose(r['0']['mean_AP'],(1+1/3)/2)
    p={'bootstrap':{'seed':84017,'replicates':10000}}
    b=paired_values(np.full(64,.006),p)
    np.testing.assert_allclose(b['ci95'],[.006,.006],atol=1e-16)
    assert classify(b['ci95'],.005,True).startswith('useful benefit supported')
    assert classify(b['ci95'],.005,False).startswith('numerically unresolved')
    assert 'ruled out' in classify([-.001,.003],.005,True)


def test_new_recorder_preserves_eulerian_heun_and_cannot_start_test_without_lock(tmp_path,monkeypatch):
    import pytest
    from record_history import run,Flow,initial,validate
    c={'n':10,'nu':.02,'dt':.0025,'save_dt':.025,'end':.05,'markers':24}
    path=run(771,c,tmp_path/'data');meta=validate(path,c)
    f=Flow(10,.02,workers=1);h=initial(f,771)
    expected=[]
    for step in range(21):
        if step%10==0:expected.append(f.inner(h,h)/2)
        if step<20:h=f.step(h,.0025)[0]
    np.testing.assert_array_equal([d['energy'] for d in meta['diagnostics']],expected)
    import record_history
    (tmp_path/'protocol.json').write_bytes((HERE/'protocol.json').read_bytes())
    monkeypatch.setattr(record_history,'HERE',tmp_path)
    with pytest.raises(FileNotFoundError):run(53000,c,tmp_path/'forbidden-data')
    assert not (tmp_path/'forbidden-data').exists()


def test_fit_split_scaler_and_model_round_trip(tmp_path):
    import joblib
    from fit_models import fit_case
    rng=np.random.default_rng(71);group=np.repeat([1,2,3,4],40);anchor=np.tile(np.repeat([.3,.4],20),4)
    A=rng.normal(size=(160,31));H=rng.normal(size=(160,23));A[group>=3]+=5
    y=A[:,0]+.2*H[:,0]+rng.normal(size=160)*.3
    d={'present':A,'history':H,'reverse':H[:,::-1],'irrelevant':H[::-1],'mirrored':H,'mirrored_reverse':H,
       'target':y,'group':group,'anchor':anchor,'particle':np.tile(np.arange(40),4),'position':A[:,:3]}
    p={'train_seeds':[1,2],'validation_seeds':[3,4]}
    case={'name':'fixture','quantile':.5,'options':{}}
    result=fit_case(d,p,case,tmp_path)
    for name in ['A_full','B_history']:
        info=result['models'][name];models=joblib.load(tmp_path/info['file']);X=A if name=='A_full' else np.column_stack([A,H])
        for model in models:np.testing.assert_allclose(model.steps[0][1].mean_,X[group<3].mean(0),atol=1e-15)
        assert info['input_dimensions']==X.shape[1]
        assert result['train_event_cutoff']==float(np.quantile(y[group<3],.5))
    from numeric_audit import model_screens
    protocol=json.loads((HERE/'protocol.json').read_text())
    models={'primary':{name:joblib.load(tmp_path/info['file']) for name,info in result['models'].items()}}
    screens=model_screens(d,d,protocol,{'cases':{'primary':result}},models,'identical fixture')
    assert all(s['pass'] for s in screens) and all(s['value']==0 for s in screens)


def test_entire_run_splits_bounded_budget_and_pretest_claim_guard(tmp_path,monkeypatch):
    import pytest
    import fit_models
    p=json.loads((HERE/'protocol.json').read_text())
    groups=[set(p[k]) for k in ['train_seeds','validation_seeds','test_seeds']]+[set(p['pilot']['seeds'])]
    assert all(not a&b for i,a in enumerate(groups) for b in groups[i+1:])
    assert [len(x) for x in groups]==[32,16,64,4]
    assert p['max_new_simulations']==22+96+64+4+8
    monkeypatch.setattr(fit_models,'HERE',tmp_path)
    fit_models.assert_test_absent(p)
    (tmp_path/'claims').mkdir();(tmp_path/'claims'/'s53000_unstarted.json').write_text('{}')
    with pytest.raises(ValueError,match='fitting forbidden'):fit_models.assert_test_absent(p)
