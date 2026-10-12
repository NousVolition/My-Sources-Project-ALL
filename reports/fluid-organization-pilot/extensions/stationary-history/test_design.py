import numpy as np
from features import inputs,target,shuffle,derange,NAMES
from models import fit,predict
from study_statistics import donor_contrast
from solver import Flow,initial
from simulate import response

def fixture():
 rng=np.random.default_rng(4)
 return dict(history=rng.normal(size=(4,4,len(NAMES))),strengths=np.array([.25,.75]),events=np.array([32,36,40,44]),initial_rates=rng.normal(size=(4,2,3)),response_time=np.broadcast_to(np.linspace(0,1.5,31),(4,31)),energy_ratio=np.ones((4,31,2)))
def test_future_measurements_do_not_change_inputs():
 z=fixture();a=inputs(z);z['energy_ratio'][:]=999
 assert all(np.array_equal(x,y) for x,y in zip(a,inputs(z)))
def test_normalized_energy_integral():
 z=fixture();assert np.allclose(target(z),np.full(8,1.5))
def test_whole_run_shuffle():
 h=np.arange(5*8*3).reshape(40,3);permutation=derange(5,9);out=shuffle(h,5,9)
 assert np.array_equal(out.reshape(5,8,3),h.reshape(5,8,3)[permutation]);assert np.all(permutation!=np.arange(5))
def test_train_only_scaling_and_prediction():
 x=np.arange(120.).reshape(40,3);y=np.sin(x[:,0]);m=fit(x,y,x[:8]+1000,y[:8],dict(bandwidth_factors=[1],ridge=[.1]),8)
 assert np.allclose(m['mean'],x.mean(0));assert np.allclose(m['scale'],x.std(0));assert np.isfinite(predict(m,x)).all()
def test_donor_bootstrap_excludes_original_self_donors():
 m=np.ones((8,8));np.fill_diagonal(m,100000);r=donor_contrast(np.zeros(8),m,500)
 assert r['null_mean_run_rmse']==1.;assert np.allclose(r['difference_interval'],[-1,-1])
def test_response_shapes_and_budget():
 f=Flow(24);h=initial(f,4);c=dict(strengths=[.25,.75],impulse_width=1.,response_dt=.005,horizon=.05,save_dt=.025)
 z=response(f,h,c)
 assert z['energy_ratio'].shape==(3,2);assert np.allclose(z['energy_ratio'][0],1.)
 assert np.max(z['response_checks'][:,2:])<1e-7
 assert np.allclose(z['energy_ratio'][-1]-1,z['difference_budget'][-1,:,0]-z['difference_budget'][-1,:,1]+z['difference_budget'][-1,:,2],atol=1e-7)

def test_independent_splits_and_declared_controls():
 from common import protocol
 p=protocol();sets=[set(p['splits'][s]) for s in ['train','validation','test']]
 assert not(sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2]);assert set(p['refinement_seeds'])<=sets[2]
 assert not set.union(*sets)&set(range(300000,300100))
 for t in p['events']:
  assert all(lag>=0 and t-lag>=0 for lag in p['history_lags'])

def test_null_shapes_keep_all_times_and_strengths_together():
 from forecast import null_errors
 x=np.arange(48.).reshape(16,3);y=np.sin(x[:,0]);m=fit(x,y,x[:8],y[:8],dict(bandwidth_factors=[1],ridge=[.1]),8)
 a=x[:,:2];h=x[:,2:];matrix=null_errors({'shuffled_0':m},a,h,y,8)
 assert matrix.shape==(1,2,2);assert np.isfinite(matrix).all()

def stationary_fixture():
 rows=[dict(t=float(t),energy=1.,dissipation=.2,integral_length=1.5,re_lambda=30.,flux_k4=.05,component_energies=[2/3]*3,turnover=1.5) for t in np.arange(0,44.5,.5)]
 import copy
 return [dict(baseline=copy.deepcopy(rows),window_spectra=[[0,.6,.3,.1],[0,.6,.3,.1]]) for _ in range(12)]

def test_stationarity_accepts_constant_statistics():
 from common import protocol
 from study_statistics import stationarity
 s=stationarity(stationary_fixture(),protocol());assert s['passed'];assert s['spectrum']['total_variation']==0

def test_stationarity_rejects_energy_or_structural_drift():
 from common import protocol
 from study_statistics import stationarity
 records=stationary_fixture()
 for r in records:
  for row in r['baseline']:
   if row['t']>=32:row['energy']=1.5
 assert not stationarity(records,protocol())['passed']
 records=stationary_fixture()
 for r in records:r['window_spectra'][1]=[0,.3,.3,.4]
 s=stationarity(records,protocol());assert not s['passed'];assert not s['spectrum']['pass_check']
 records=stationary_fixture()
 for r in records:
  for row in r['baseline']:row['flux_k4']=1e-14
 assert not stationarity(records,protocol())['passed']
