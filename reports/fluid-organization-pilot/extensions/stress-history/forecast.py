import argparse,json,datetime
import numpy as np
from threadpoolctl import threadpool_limits
from common import ROOT,save,table,digest,read_protocol,numerical_checks
from features import inputs,target,derange,shuffle_groups
from models import fit,predict
def load(regime,split):
    p=read_protocol();cc=[];pp=[];hh=[];yy=[];hashes={}
    for seed in p['splits'][split]:
        path=ROOT/'data/runs'/f'{regime}_base_{seed}.npz'
        meta=json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
        if not all(numerical_checks(meta,p['gates']).values()):raise ValueError('Numerical gate failed '+path.name)
        if digest(path)!=meta['sha256']:raise ValueError('Recording hash mismatch')
        with np.load(path) as z:
            c,a,h=inputs(z,p['history_lags']);y=target(z)
        cc.extend(c);pp.extend(a);hh.extend(h);yy.extend(y);hashes[path.name]=digest(path)
    return *(np.asarray(v) for v in (cc,pp,hh,yy)),hashes
def families(c,p,h,regime,split):
    proto=read_protocol();out={'conventional':c,'present':p,'history':np.c_[p,h]}
    for j in range(proto['shuffled_controls']):
        perm=derange(len(h)//2,proto['shuffle_seed']+100*j+10*list(proto['regimes']).index(regime)+{'train':0,'validation':1,'test':2}[split])
        out['shuffled_'+str(j)]=np.c_[p,shuffle_groups(h,perm)]
    return out
def freeze():
    if (ROOT/'results/frozen_models.json').exists():raise FileExistsError('Freeze already exists; never refit against test outcomes')
    p=read_protocol();out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':digest(ROOT/'protocol.json'),'source_sha256':{s:digest(ROOT/s) for s in ['features.py','feature_names.json','models.py','simulate.py','physics.py','vendor/numerics.py','forecast.py']},'data_sha256':{},'models':{},'prior_by_strength':{}}
    for regime in p['regimes']:
        c,a,h,y,hashes=load(regime,'train');cv,av,hv,yv,hashv=load(regime,'validation');out['data_sha256'].update(hashes);out['data_sha256'].update(hashv)
        train=families(c,a,h,regime,'train');val=families(cv,av,hv,regime,'validation')
        out['prior_by_strength'][regime]=y.reshape(-1,2).mean(0).tolist()
        out['models'][regime]={name:fit(x,y,val[name],yv,p['model_grid']) for name,x in train.items()}
        print(regime+' validation '+str({k:round(v['validation_rmse'],6) for k,v in out['models'][regime].items()}),flush=True)
    save('results/frozen_models.json',out);save('results/freeze_receipt.json',{'sha256':digest(ROOT/'results/frozen_models.json'),'created_utc':out['created_utc'],'protocol_sha256':out['protocol_sha256'],'test_simulations_present_at_freeze':[str(q.name) for q in (ROOT/'data/runs').glob('*') if any('_'+str(s)+'.' in q.name for s in p['splits']['test'])]})
def verify_freeze():
    q=json.loads((ROOT/'results/freeze_receipt.json').read_text(encoding='utf-8'));m=json.loads((ROOT/'results/frozen_models.json').read_text(encoding='utf-8'))
    assert not q['test_simulations_present_at_freeze'];assert digest(ROOT/'results/frozen_models.json')==q['sha256'];assert digest(ROOT/'protocol.json')==m['protocol_sha256']
    for s,h in m['source_sha256'].items():assert digest(ROOT/s)==h
    for s,h in m['data_sha256'].items():assert digest(ROOT/'data/runs'/s)==h
    return m
def bootstrap(a,b,p,seed):
    rng=np.random.default_rng(seed);idx=rng.integers(0,len(a),(p['bootstrap_replicates'],len(a)))
    diff=(a-b)[idx].mean(1);relative=100*(1-a[idx].mean(1)/b[idx].mean(1))
    q=p['interval_quantiles']
    return {'mean_rmse_difference':float(np.mean(a-b)),'rmse_difference_interval':np.quantile(diff,q).tolist(),'relative_improvement_percent':float(100*(1-a.mean()/b.mean())),'relative_improvement_interval':np.quantile(relative,q).tolist(),'fields_better':int(np.sum(a<b)),'independent_fields':len(a)}
def null_errors(models,present,history,y):
    """All wrong whole-field history donors; retain both strength rows together."""
    n=len(y)//2;mat=[];h=history.reshape(n,2,-1);a=present.reshape(n,2,-1);yy=y.reshape(n,2)
    for name,model in models.items():
        if not name.startswith('shuffled'):continue
        errors=np.empty((n,n))
        for i in range(n):
            x=np.concatenate([np.tile(a[i][None],(n,1,1)),h],axis=2).reshape(2*n,-1)
            forecast=predict(model,x).reshape(n,2);errors[i]=np.sqrt(np.mean((forecast-yy[i])**2,axis=1))
        mat.append(errors)
    return np.array(mat)
def bootstrap_null(history_error,matrix,p,seed):
    """Resample independent field nodes jointly in current and donor roles.

    The null is an off-diagonal U-statistic, not independent rows: every donor
    also belongs to the current-field sample. Duplicate copies of the same
    original field are excluded as self donors, not just matrix diagonals.
    """
    n=len(history_error);mask=~np.eye(n,dtype=bool);base=float(matrix[mask].mean());perfield=np.array([matrix[i,np.arange(n)!=i].mean() for i in range(n)]);rng=np.random.default_rng(seed);diff=[];relative=[]
    for offset in range(0,p['bootstrap_replicates'],500):
        size=min(500,p['bootstrap_replicates']-offset);idx=rng.integers(0,n,(size,n));valid=idx[:,:,None]!=idx[:,None,:];denom=valid.sum((1,2))
        if np.any(denom==0):raise ValueError('Degenerate node-bootstrap sample; larger field ensemble required')
        b=(matrix[idx[:,:,None],idx[:,None,:]]*valid).sum((1,2))/denom;a=history_error[idx].mean(1);diff.extend(a-b);relative.extend(100*(1-a/b))
    return {'mean_rmse_difference':float(history_error.mean()-base),'rmse_difference_interval':np.quantile(diff,p['interval_quantiles']).tolist(),'relative_improvement_percent':float(100*(1-history_error.mean()/base)),'relative_improvement_interval':np.quantile(relative,p['interval_quantiles']).tolist(),'fields_better':int(np.sum(history_error<perfield)),'independent_fields':n,'null_mean_field_rmse':base,'null_definition':'Mean RMSE over all nonidentity whole-run history donors and the three frozen shuffled-history models','uncertainty':'Field-node bootstrap resamples both current and history-donor roles jointly; excludes same original seed even when duplicated'}
def evaluate():
    p=read_protocol();m=verify_freeze();rows=[];summary={};null_record={}
    for i,regime in enumerate(p['regimes']):
        c,a,h,y,hashes=load(regime,'test');x=families(c,a,h,regime,'test');pred={k:predict(model,x[k]) for k,model in m['models'][regime].items()};pred['prior_by_strength']=np.tile(m['prior_by_strength'][regime],len(y)//2);errors={k:np.sqrt(((v-y)**2).reshape(-1,2).mean(1)) for k,v in pred.items()}
        matrix=null_errors(m['models'][regime],a,h,y);null_record[regime]=matrix
        for j,seed in enumerate(p['splits']['test']):
            for k,strength in enumerate(p['kick_rms_over_baseline_rms']):rows.append({'regime':regime,'seed':seed,'relative_strength':strength,'observed':float(y[2*j+k]),**{name:float(v[2*j+k]) for name,v in pred.items()}})
        contrasts={'history_vs_present':bootstrap(errors['history'],errors['present'],p,1200+i),'history_vs_shuffled_mean':bootstrap_null(errors['history'],matrix.mean(0),p,1300+i),'present_vs_conventional':bootstrap(errors['present'],errors['conventional'],p,1400+i),'history_vs_training_prior':bootstrap(errors['history'],errors['prior_by_strength'],p,1500+i)}
        summary[regime]={'mean_field_rmse':{k:float(v.mean()) for k,v in errors.items()},'pooled_rmse':{k:float(np.sqrt(np.mean((v-y)**2))) for k,v in pred.items()},'target_mean':float(y.mean()),'target_std':float(y.std()),'contrasts':contrasts,'history_supported':all(contrasts[name]['rmse_difference_interval'][1]<0 for name in ['history_vs_present','history_vs_shuffled_mean','history_vs_training_prior']),'test_data_sha256':hashes}
        print(regime+' '+str(summary[regime]['contrasts']['history_vs_present']),flush=True)
    folder=ROOT/'data/forecast';folder.mkdir(exist_ok=True);np.savez_compressed(folder/'shuffled_errors.npz',**null_record)
    table('predictions.csv',rows);save('results/prediction_summary.json',summary)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','evaluate']);a=ap.parse_args()
    with threadpool_limits(limits=1):(freeze if a.action=='freeze' else evaluate)()
