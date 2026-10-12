"""Verify delivered bytes and reconstruct the frozen forecasts without rerunning DNS."""
import json,datetime
import numpy as np
from threadpoolctl import threadpool_limits
from common import ROOT,digest,read,protocol
from forecast import verify_freeze,load,families,null_errors
from features import inputs,target
from models import predict,run_errors
from study_statistics import contrast,donor_contrast,stationarity

def close(a,b):
 if isinstance(a,dict):
  assert set(a)==set(b)
  for k in a:close(a[k],b[k])
 elif isinstance(a,(list,tuple)):
  assert len(a)==len(b)
  for x,y in zip(a,b):close(x,y)
 elif isinstance(a,(int,float)) and not isinstance(a,bool):assert np.isclose(a,b,rtol=1e-9,atol=1e-11),(a,b)
 else:assert a==b,(a,b)

def audit():
 manifest=read('manifest.json')
 for f in manifest['files']:
  path=ROOT/f['path'];assert path.stat().st_size==f['size'] and digest(path)==f['sha256'],f['path']
 p=protocol();m=verify_freeze();c,a,h,y,_,meta=load('test');n=len(p['splits']['test']);rows=len(p['events'])*len(p['strengths']);x=families(c,a,h,n,'test')
 frozen_at=datetime.datetime.fromisoformat(m['created_utc'])
 for record in meta:assert datetime.datetime.fromisoformat(record['started_utc'])>frozen_at
 registration_path=ROOT/'results/public_protocol_receipt.json'
 if registration_path.exists():
  registration=read('results/public_protocol_receipt.json');assert registration['frozen_models_sha256']==digest(ROOT/'results/frozen_models.json')
  published_at=datetime.datetime.fromisoformat(registration['published_utc'])
  for record in meta:assert datetime.datetime.fromisoformat(record['started_utc'])>published_at
 pred={name:predict(model,x[name]) for name,model in m['models'].items()};pred['training_prior']=np.tile(m['prior_by_strength'],n*len(p['events']));pred['event_prior']=np.tile(m['prior_by_event_strength'],n)
 matrix=null_errors(m['models'],a,h,y,rows)
 with np.load(ROOT/'data/predictions.npz') as z:
  assert np.allclose(y,z['truth'],rtol=1e-12,atol=1e-12);assert np.allclose(matrix,z['null_errors'],rtol=1e-9,atol=1e-11)
  for key,value in pred.items():assert np.allclose(value,z[key],rtol=1e-9,atol=1e-11),key
 errors={key:run_errors(v,y,rows) for key,v in pred.items()};r=read('results/prediction_summary.json')
 for j,name in enumerate(['global_present','present','training_prior','event_prior']):close(contrast(errors['history'],errors[name],p['bootstrap_replicates'],901+j),r['contrasts'][name])
 close(donor_contrast(errors['history'],matrix.mean(0),p['bootstrap_replicates'],950),r['contrasts']['wrong_history']);close(stationarity(meta,p),read('results/test_stationarity.json'))
 # Each saved history is reconstructed from genuinely pre-intervention baseline data.
 for split in p['splits']:
  for seed in p['splits'][split]:
   with np.load(ROOT/'data/runs'/f'base_{seed}.npz') as z:
    for j,t in enumerate(p['events']):
     times=np.array([t]+[t-lag for lag in p['history_lags']]);idx=[int(np.argmin(abs(z['baseline_time']-tt))) for tt in times]
     assert np.allclose(z['baseline_time'][idx],times,atol=1e-10);assert np.array_equal(z['history'][j],z['baseline_features'][idx])
 from checks import verify_all
 numerical=verify_all(write=False);close(numerical,read('results/numerical_checks.json'));assert numerical['all_passed']
 print(json.dumps(dict(manifest_files_verified=len(manifest['files']),independent_test_runs=n,frozen_forecasts_rebuilt=True,all_history_times_causal=True)))
if __name__=='__main__':
 with threadpool_limits(limits=1):audit()
