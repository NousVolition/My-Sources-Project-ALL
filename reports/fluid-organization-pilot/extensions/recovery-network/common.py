"""Small reproducible statistics and first-arrival/sustained-arrival helpers."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent
P = json.loads((ROOT/'protocol.json').read_text(encoding='utf-8'))

def save(path, obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def arrivals(t, values, threshold, residence):
    """Grid-observed crossing; never count an incomplete residence window."""
    t=np.asarray(t); hit=np.asarray(values)<=threshold
    first=float(t[np.flatnonzero(hit)[0]]) if hit.any() else None
    sustained=None
    for i in np.flatnonzero(hit):
        j=np.searchsorted(t,t[i]+residence-1e-9)
        if j<len(t) and hit[i:j+1].all():
            sustained=float(t[i]);break
    return first,sustained

def interval(values):
    x=np.asarray(values,float);rng=np.random.default_rng(P['inference']['bootstrap_seed'])
    samples=x[rng.integers(0,len(x),(P['inference']['bootstrap_reps'],len(x)))].mean(1)
    return {'mean':float(x.mean()),'ci95':np.quantile(samples,[.025,.975]).tolist(),'n_independent_runs':len(x)}

def ridge(train,val,test,yt,yv):
    mean=train.mean(0);scale=train.std(0);scale[scale<1e-10]=1
    a,b,c=[(v-mean)/scale for v in (train,val,test)];ym=yt.mean(0)
    u,s,vt=np.linalg.svd(a,full_matrices=False);best=None
    for alpha in P['inference']['alphas']:
        w=vt.T@((s/(s*s+alpha))[:,None]*(u.T@(yt-ym)))
        loss=float(np.mean((b@w+ym-yv)**2))
        if best is None or loss<best[0]:best=(loss,alpha,c@w+ym)
    return best[2],{'alpha':best[1],'validation_rmse':float(np.sqrt(best[0])),'features':train.shape[1]}

def prediction_comparison(parts, target_names):
    """Whole-run split, train-only transforms, paired run bootstrap."""
    out={}; saved={};rng=np.random.default_rng(441)
    transformed={}
    m=parts['train']['current'].mean(0);s=parts['train']['current'].std(0);s[s<1e-10]=1
    for split,d in parts.items():
        a=(d['current']-m)/s
        # Same fixed present-state quadratic map in every predictor.
        i,j=np.triu_indices(a.shape[1]);cur=np.c_[a,a[:,i]*a[:,j]]
        history=d['history'];shuffled=history.copy()
        # Permute whole histories between runs at the same intervention/observation slot.
        for slot in np.unique(d['slot']):
            ix=np.flatnonzero(d['slot']==slot);shuffled[ix]=history[rng.permutation(ix)]
        transformed[split]={'current':cur,'history':np.c_[cur,history],
                            'shuffled_history':np.c_[cur,shuffled]}
    for j,name in enumerate(target_names):
        scores={};out[name]={}
        for model in ('current','history','shuffled_history'):
            pred,meta=ridge(*(transformed[k][model] for k in ('train','validation','test')),
                            parts['train']['y'][:,j:j+1],parts['validation']['y'][:,j:j+1])
            truth=parts['test']['y'][:,j];pred=pred[:,0];groups=parts['test']['groups']
            err=np.array([np.sqrt(np.mean((pred[groups==g]-truth[groups==g])**2)) for g in np.unique(groups)])
            scores[model]=err;out[name][model]={**meta,**interval(err),'per_run_rmse':err.tolist()}
            saved[name+'_'+model]=pred
        out[name]['history_minus_current']=interval(scores['history']-scores['current'])
        out[name]['history_minus_shuffled']=interval(scores['history']-scores['shuffled_history'])
    return out,saved
