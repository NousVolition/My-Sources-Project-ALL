import argparse,json
import numpy as np
from threadpoolctl import threadpool_limits
from common import ROOT,digest,save,read,protocol,now,sources
from features import inputs,target,shuffle
from models import fit,predict,run_errors
from study_statistics import contrast,donor_contrast,stationarity

LOCKED=['protocol.json','solver.py','kernels.py','simulate.py','features.py','models.py','study_statistics.py','forecast.py','common.py','checks.py']
def load(split):
 p=protocol();c=[];a=[];h=[];y=[];hashes={};meta=[]
 from checks import numerical
 for seed in p['splits'][split]:
  path=ROOT/'data/runs'/f'base_{seed}.npz';q=read(f'data/runs/base_{seed}.json');assert digest(path)==q['sha256'];assert q['source_sha256']==sources();checks=numerical(path,q,p);assert all(checks.values()),(seed,checks)
  with np.load(path) as z:
   cc,aa,hh=inputs(z);c.extend(cc);a.extend(aa);h.extend(hh);y.extend(target(z))
   window_spectra=[]
   for left,right in p['stationarity_windows']:
    sp=z['spectra'][(z['baseline_time']>=left)&(z['baseline_time']<right)];window_spectra.append((sp/sp.sum(1)[:,None]).mean(0).tolist())
   q['window_spectra']=window_spectra
  hashes[path.name]=digest(path);hashes[path.with_suffix('.json').name]=digest(path.with_suffix('.json'));meta.append(q)
 return *(np.asarray(v) for v in [c,a,h,y]),hashes,meta
def families(c,a,h,n,split):
 p=protocol();x=dict(global_present=c,present=a,history=np.c_[a,h])
 for j in range(p['shuffled_models']):x[f'shuffled_{j}']=np.c_[a,shuffle(h,n,801+100*j+{'train':0,'validation':1,'test':2}[split])]
 return x
def freeze():
 p=protocol();path=ROOT/'results/frozen_models.json';assert not path.exists(),'Freeze exists; no refitting'
 from checks import screening
 screen=screening();save('results/pilot_response_checks.json',screen);assert screen['passed'],screen
 test_files=[q.name for q in (ROOT/'data/runs').glob('*') if any(f'_{s}' in q.name for s in p['splits']['test'])];assert not test_files
 c,a,h,y,hashes,meta=load('train');cv,av,hv,yv,hvhash,mv=load('validation');stat=stationarity(meta+mv,p);save('results/calibration_stationarity.json',stat);assert stat['passed'],stat
 n=len(p['splits']['train']);nv=len(p['splits']['validation']);rows=len(p['events'])*len(p['strengths']);x=families(c,a,h,n,'train');v=families(cv,av,hv,nv,'validation')
 m={name:fit(value,y,v[name],yv,p['model_grid'],rows) for name,value in x.items()}
 out=dict(created_utc=now(),source_sha256={s:digest(ROOT/s) for s in LOCKED},data_sha256={**hashes,**hvhash},pilot_data_sha256=screen['sha256'],models=m,prior_by_strength=y.reshape(n,len(p['events']),len(p['strengths'])).mean((0,1)).tolist(),prior_by_event_strength=y.reshape(n,rows).mean(0).tolist(),absent_test_files_at_freeze=test_files)
 save('results/frozen_models.json',out);save('results/freeze_receipt.json',dict(created_utc=out['created_utc'],sha256=digest(path)))
 print(json.dumps({k:{z:v[z] for z in ['kind','ridge','factor','validation_mean_run_rmse']} for k,v in m.items()}),flush=True)
def verify_freeze():
 m=read('results/frozen_models.json');r=read('results/freeze_receipt.json');assert digest(ROOT/'results/frozen_models.json')==r['sha256'];assert not m['absent_test_files_at_freeze']
 for path,sha in m['source_sha256'].items():assert digest(ROOT/path)==sha,path
 for path,sha in m['data_sha256'].items():assert digest(ROOT/'data/runs'/path)==sha,path
 for path,sha in m['pilot_data_sha256'].items():assert digest(ROOT/path)==sha,path
 return m
def null_errors(models,present,history,y,rows):
 n=len(y)//rows;a=present.reshape(n,rows,-1);h=history.reshape(n,rows,-1);yy=y.reshape(n,rows);result=[]
 for name,m in models.items():
  if not name.startswith('shuffled_'):continue
  matrix=np.empty((n,n))
  for i in range(n):
   x=np.concatenate([np.broadcast_to(a[i],(n,rows,a.shape[-1])),h],axis=2).reshape(n*rows,-1)
   pred=predict(m,x).reshape(n,rows);matrix[i]=np.sqrt(((pred-yy[i])**2).mean(1))
  result.append(matrix)
 return np.array(result)
def evaluate():
 p=protocol();m=verify_freeze();c,a,h,y,hashes,meta=load('test');n=len(p['splits']['test']);rows=len(p['events'])*len(p['strengths']);x=families(c,a,h,n,'test');pred={name:predict(model,x[name]) for name,model in m['models'].items()};pred['training_prior']=np.tile(m['prior_by_strength'],n*len(p['events']));pred['event_prior']=np.tile(m['prior_by_event_strength'],n);errors={name:run_errors(v,y,rows) for name,v in pred.items()};matrix=null_errors(m['models'],a,h,y,rows)
 comparisons={name:contrast(errors['history'],errors[name],p['bootstrap_replicates'],901+j) for j,name in enumerate(['global_present','present','training_prior','event_prior'])};comparisons['wrong_history']=donor_contrast(errors['history'],matrix.mean(0),p['bootstrap_replicates'],950)
 supported=all(q['difference_interval'][1]<0 for q in comparisons.values());stat=stationarity(meta,p)
 save('results/prediction_summary.json',dict(independent_test_runs=n,snapshots_per_run=len(p['events']),strengths=p['strengths'],target_mean=float(y.mean()),target_std=float(y.std()),mean_run_rmse={k:float(v.mean()) for k,v in errors.items()},contrasts=comparisons,history_supported=bool(supported),stationarity_pass=stat['passed'],test_data_sha256=hashes))
 save('results/test_stationarity.json',stat)
 dest=ROOT/'data/predictions.npz';np.savez_compressed(dest,truth=y,**pred,null_errors=matrix)
 import csv
 with (ROOT/'results/predictions.csv').open('w',newline='',encoding='utf-8') as fh:
  writer=csv.DictWriter(fh,fieldnames=['seed','event','strength','target',*pred]);writer.writeheader()
  for i,s in enumerate(p['splits']['test']):
   for j,t in enumerate(p['events']):
    for k,alpha in enumerate(p['strengths']):
     idx=i*rows+j*len(p['strengths'])+k;writer.writerow(dict(seed=s,event=t,strength=alpha,target=y[idx],**{name:v[idx] for name,v in pred.items()}))
 print(json.dumps({'history_supported':supported,'stationarity_pass':stat['passed'],'contrasts':comparisons}),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','evaluate']);a=ap.parse_args()
 with threadpool_limits(limits=1):(freeze if a.action=='freeze' else evaluate)()
