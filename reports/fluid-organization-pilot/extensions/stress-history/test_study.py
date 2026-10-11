import json
import numpy as np
import pytest
from common import ROOT
from physics import PilotFlow,initial,impulse,energy,L
from simulate import stage,advance,observe,markers,geometry,marker_sample
from features import inputs,shuffle_groups,derange
from models import fit,predict
def test_rotational_advection_conserves_energy():
    f=PilotFlow(24,.04);h=initial(f,83);u=f.real(h);w=f.real(f.curl(h));a=f.project(f.hat(np.cross(u,w,axisa=0,axisb=0,axisc=0)))
    assert abs(f.inner(h,a))<1e-12*f.inner(h,h)
def test_incompressibility_and_strict_dealiasing():
    f=PilotFlow(24,.04);h=initial(f,9)+impulse(f,.75,1)
    assert abs(f.real(1j*sum(k*a for k,a in zip(f.k,h)))).max()<1e-12
    assert np.max(abs(h[:,~f.keep]))==0
def test_beltrami_decay_against_exact_solution():
    f=PilotFlow(24,.2);h=initial(f,3,kind='beltrami');x=markers(3,2,.18);old=h.copy();dt=.02
    for _ in range(10):hs,xs,*_=advance(f,[h],[x],dt);h,x=hs[0],xs[0]
    exact=old*np.exp(-f.nu*(2*np.pi/L)**2*.2)
    assert np.sqrt(f.inner(h-exact,h-exact)/f.inner(exact,exact))<1e-9
def test_identical_clones_do_not_diverge():
    f=PilotFlow(24,.04);h=initial(f,7);x=markers(7,2,.18)
    hs,xs,b,pd,_=advance(f,[h,h.copy()],[x,x.copy()],.005)
    assert np.array_equal(hs[0],hs[1]);assert np.array_equal(xs[0],xs[1]);assert np.array_equal(pd,np.zeros((1,2)))
def test_difference_production_matches_strain_identity():
    f=PilotFlow(24,.04);h=initial(f,83);d=impulse(f,.5,1);hp=h+d
    diss=f.nu*f.inner(f.curl(d),f.curl(d))/L**3
    production=f.inner(d,f.rhs_only(hp)-f.rhs_only(h))/L**3+diss
    du=f.real(d);g=f.gradient(h);direct=-np.mean(np.einsum('iabc,jabc,ijabc->abc',du,du,g))
    assert abs(production-direct)<1e-12
def test_budget_closes_over_saved_interval():
    f=PilotFlow(24,.04);h=initial(f,4);d=impulse(f,.75,1);xs=[markers(4,2,.18)]*2;hs=[h,h+d];e0=energy(f,d);eb=np.array([energy(f,q) for q in hs]);b=np.zeros((2,2));pd=np.zeros((1,2))
    for _ in range(20):hs,xs,db,dp,_=advance(f,hs,xs,.0025);b+=db;pd+=dp
    assert np.max(abs(np.array([energy(f,q) for q in hs])-eb-b[:,0]+b[:,1]))<1e-8
    assert abs(energy(f,hs[1]-hs[0])-e0-pd[0,0]+pd[0,1])<1e-8
def test_local_strain_and_vorticity_measurement():
    f=PilotFlow(24,.04);h=initial(f,83);obs,_=observe(f,h,markers(83,2,.18),.18);g=f.gradient(h)[:,:,12,12,12];w=f.real(f.curl(h))[:,12,12,12]
    assert np.isclose(obs[14],w@((g+g.T)/2)@w,atol=1e-12)
    assert np.allclose(obs[-9:].reshape(3,3),g,atol=1e-12);assert np.allclose(obs[-12:-9],f.real(h)[:,12,12,12],atol=1e-12)
def test_same_initial_field_on_different_grids():
    a=PilotFlow(24,.04);b=PilotFlow(48,.04);ha=initial(a,97);hb=initial(b,97)
    assert np.allclose(a.real(ha),b.real(hb)[:,::2,::2,::2],atol=1e-13)
    assert np.allclose(a.real(impulse(a,.75,1)),b.real(impulse(b,.75,1))[:,::2,::2,::2],atol=1e-13)
def test_geometry_is_translation_and_rotation_invariant():
    x=markers(8,5,.18);q=np.linalg.qr(np.random.default_rng(5).normal(size=(3,3)))[0]
    assert np.allclose(geometry(x,.18),geometry(x@q+3,.18),atol=1e-12)
def test_inputs_ignore_future_recordings():
    rng=np.random.default_rng(2);z={'history_time':np.linspace(-2,0,41),'history':rng.normal(size=(41,21)),'relative_amplitudes':np.array([.25,.75]),'initial_rates':rng.normal(size=(2,2)),'initial_perturbation_energy':np.ones(2),'positions':np.zeros((31,3,16,4,3)),'perturbation_energy':np.ones((31,2))}
    a=inputs(z,[.25,.5,1]);z['positions'][:]=999;z['perturbation_energy'][:]=999
    assert all(np.array_equal(x,y) for x,y in zip(a,inputs(z,[.25,.5,1])))
def test_shuffling_keeps_strengths_together_and_has_no_fixed_groups():
    h=np.arange(60).reshape(10,6);perm=derange(5,3);out=shuffle_groups(h,perm)
    assert np.all(perm!=np.arange(5));assert np.array_equal(out.reshape(5,2,6),h.reshape(5,2,6)[perm])
def test_model_scaling_uses_training_only():
    x=np.arange(40.).reshape(20,2);y=np.sin(x[:,0]);m=fit(x,y,x[:4]+100,y[:4],{'ridge':[.1],'bandwidth_factors':[1]})
    assert np.allclose(m['mean'],x.mean(0));assert np.allclose(m['scale'],x.std(0));assert np.isfinite(predict(m,x)).all()
def test_locked_run_splits_do_not_overlap():
    path=ROOT/'protocol.json'
    if not path.exists():return
    p=json.loads(path.read_text(encoding='utf-8'));a,b,c=[set(p['splits'][k]) for k in ['train','validation','test']]
    assert not(a&b or a&c or b&c);assert not((a|b|c)&set(p['pilot_seeds']));assert set(p['refinement_seeds'])<=c
def test_future_history_is_rejected():
    z={'history_time':np.array([-.5,0,.25]),'history':np.ones((3,21))}
    with pytest.raises(ValueError,match='past observations'):inputs(z,[.25,.5,1])
def test_cubic_marker_velocity_accuracy_and_periodicity():
    from scipy.ndimage import map_coordinates
    n=24;dx=L/n;a=np.arange(n)*dx-L/2;xg,yg,zg=np.meshgrid(a,a,a,indexing='ij');u=np.stack([np.sin(2*2*np.pi*xg/L),np.cos(3*2*np.pi*yg/L),np.zeros_like(zg)])
    x=np.random.default_rng(41).uniform(-L/2,L/2,(8,4,3));exact=np.stack([np.sin(2*2*np.pi*x[...,0]/L),np.cos(3*2*np.pi*x[...,1]/L),np.zeros(x.shape[:-1])],axis=-1)
    cubic=marker_sample(u,x,dx);coords=(((x.reshape(-1,3)+L/2)%L)/dx).T;linear=np.stack([map_coordinates(v,coords,order=1,mode='grid-wrap',prefilter=False) for v in u],axis=-1).reshape(x.shape)
    assert np.max(abs(cubic-exact))<.002;assert np.max(abs(cubic-exact))<.1*np.max(abs(linear-exact));assert np.allclose(cubic,marker_sample(u,x+L,dx),atol=1e-12)
def test_diagnostic_optimization_preserves_screened_physics():
    import importlib.util
    spec=importlib.util.spec_from_file_location('screened',ROOT/'screen_cubic_simulate.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    f=PilotFlow(24,.04);h=initial(f,72);x=markers(72,2,.18);a,da=observe(f,h,x,.18);b,db=old.observe(f,h,x,.18)
    assert np.allclose(a[:21],b,rtol=1e-13,atol=1e-13)
    assert np.allclose(da['spectrum'],db['spectrum'],rtol=1e-12,atol=1e-13);assert np.isclose(da['integral_length'],db['integral_length'],rtol=1e-13)
    hs=[h,h+impulse(f,.75,1)];xs=[x.copy(),x.copy()];current=advance(f,hs,xs,.005);screened=old.advance(f,hs,xs,.005)
    for aa,bb in zip(current[:4],screened[:4]):assert np.array_equal(np.asarray(aa),np.asarray(bb))
def test_null_bootstrap_excludes_same_field_in_both_roles():
    from forecast import bootstrap_null
    matrix=np.ones((8,8));np.fill_diagonal(matrix,10000.)
    result=bootstrap_null(np.zeros(8),matrix,{'bootstrap_replicates':500,'interval_quantiles':[.0125,.9875]},51)
    assert result['null_mean_field_rmse']==1.;assert np.allclose(result['rmse_difference_interval'],[-1.,-1.]);assert result['relative_improvement_percent']==100.
