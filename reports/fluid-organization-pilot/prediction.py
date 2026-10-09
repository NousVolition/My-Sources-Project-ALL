"""Past-only features, run-separated fitting, and held-out paired uncertainty."""
import argparse
import json
from pathlib import Path
import numpy as np
from core import ROOT,wrap,save_json


def neighbors(x,k=8):
    d=wrap(x[None,:,:]-x[:,None,:]);r2=np.sum(d*d,axis=-1)
    np.fill_diagonal(r2,np.inf)
    return np.argsort(r2,axis=1,kind='stable')[:,:k]


def edges(x,nn):
    return wrap(x[nn]-x[:,None,:])


def geometry(e):
    r=np.linalg.norm(e,axis=-1);unit=e/np.maximum(r[:,:,None],1e-15)
    C=np.einsum('mki,mkj->mij',e,e)/e.shape[1]
    components=C[:,np.triu_indices(3)[0],np.triu_indices(3)[1]]
    eig=np.linalg.eigvalsh(C)
    return np.column_stack([r.mean(1),r.std(1),r.min(1),r.max(1),components,eig])


def features(past_times,past_positions,velocity_base,velocity_delta,gradient_base,gradient_delta,k=8,lags=(.1,.2,.3)):
    """API accepts ONLY observations through the kick, never future positions."""
    t=past_times[-1];x=past_positions[-1];nn=neighbors(x,k);e=edges(x,nn)
    r=np.linalg.norm(e,axis=-1);ub=velocity_base;dp=velocity_delta
    pairrate=lambda v:np.sum(e*(v[nn]-v[:,None,:]),axis=-1)/np.maximum(r*r,1e-30)
    rb=pairrate(ub);rd=pairrate(dp)
    G=gradient_base+gradient_delta;S=(G+G.swapaxes(1,2))/2
    current=np.column_stack([ub,dp,gradient_base.reshape(len(x),9),gradient_delta.reshape(len(x),9),
                             np.linalg.eigvalsh(S),geometry(e),rb.mean(1),rb.std(1),rd.mean(1),rd.std(1)])
    history=[]
    nowgeo=geometry(e)
    for lag in lags:
        ix=int(np.argmin(abs(past_times-(t-lag))))
        if abs(past_times[ix]-(t-lag))>1e-7:raise ValueError('Missing requested history time')
        old=edges(past_positions[ix],nn);ro=np.linalg.norm(old,axis=-1)
        log=np.log(r/np.maximum(ro,1e-15))/lag
        cosine=np.sum(e*old,axis=-1)/np.maximum(r*ro,1e-15)
        history.append(np.column_stack([log.mean(1),log.std(1),log.min(1),log.max(1),
                                        cosine.mean(1),cosine.std(1),(nowgeo-geometry(old))/lag]))
    return current,np.column_stack(history),nn,rb.mean(1)+rd.mean(1),rd.mean(1)


def dataset(paths,p):
    A=[];B=[];Y=[];groups=[];frozen=[]
    for path in paths:
        z=np.load(path);meta=json.loads(path.with_suffix('.json').read_text());t=z['time'];kick=meta['config']['kick_time']
        ix=int(np.argmin(abs(t-kick)));pos=z['positions_base'];post=z['positions_perturbed']
        aa,bb,nn,fr,dr=features(t[:ix+1],pos[:ix+1],z['velocity_base'],z['velocity_delta'],z['gradient_base'],z['gradient_delta'],p['history']['neighbors'],p['history']['lags'])
        r=np.linalg.norm(edges(pos[ix],nn),axis=-1)
        future=[np.log(np.linalg.norm(edges(q[-1],nn),axis=-1)/r).mean(1)/(t[-1]-t[ix]) for q in [pos,post]]
        A.append(aa);B.append(bb);Y.append(np.column_stack([future[1],future[1]-future[0]]));groups.extend([meta['seed']]*len(aa));frozen.append(np.column_stack([fr,dr]))
    return {'A':np.concatenate(A),'B':np.concatenate(B),'y':np.concatenate(Y),'groups':np.array(groups),'frozen':np.concatenate(frozen)}


def expansion(a):
    # A fixed quadratic map is a stronger present-only comparator than linear ridge.
    # All feature sets receive exactly the same map, not test-selected complexity.
    i,j=np.triu_indices(a.shape[1])
    return np.column_stack([a,a[:,i]*a[:,j]])


def fit_predict(train,val,test,yt,yv,alphas,quadratic=True):
    mean=train.mean(0);scale=train.std(0);scale[scale<1e-10]=1
    tr=(train-mean)/scale;va=(val-mean)/scale;te=(test-mean)/scale
    if quadratic:tr,va,te=map(expansion,[tr,va,te])
    # Avoid a cubic solve of a large feature covariance: low-rank SVD shared by alphas.
    from scipy.linalg import svd
    from threadpoolctl import threadpool_limits
    xmean=tr.mean(0);xstd=tr.std(0);xstd[xstd<1e-10]=1
    tr=(tr-xmean)/xstd;va=(va-xmean)/xstd;te=(te-xmean)/xstd
    ym=float(yt.mean())
    with threadpool_limits(limits=1):
        u,s,vt=svd(tr,full_matrices=False,check_finite=False,lapack_driver='gesdd')
    scores=[];coeff=[]
    uy=u.T@(yt-ym)
    for alpha in alphas:
        w=vt.T@((s/(s*s+alpha))*uy)
        scores.append(np.mean((va@w+ym-yv)**2));coeff.append(w)
    best=int(np.argmin(scores));pred=te@coeff[best]+ym
    return pred,{'alpha':alphas[best],'validation_rmse':float(np.sqrt(scores[best])),'input_features':train.shape[1],'expanded_features':tr.shape[1]}


def rmse_by_group(y,pred,groups):
    return np.array([np.sqrt(np.mean((y[groups==s]-pred[groups==s])**2)) for s in np.unique(groups)])


def paired(a,b,reps=4000):
    delta=b-a;rng=np.random.default_rng(410)
    bs=delta[rng.integers(0,len(delta),(reps,len(delta)))].mean(1)
    return {'mean_B_minus_A':float(delta.mean()),'ci95':np.quantile(bs,[.025,.975]).tolist(),'per_seed_delta':delta.tolist()}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');ap.add_argument('--output',type=Path,default=ROOT/'results');a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);p=json.loads((ROOT/'protocol.json').read_text())
    sets={key:dataset([a.data/'ensemble'/f'seed{s}.npz' for s in p[key+'_seeds']],p) for key in ['train','validation','test']}
    rng=np.random.default_rng(713);controls={}
    for key,d in sets.items():
        hist=d['B'].copy()
        for s in np.unique(d['groups']):
            ix=np.flatnonzero(d['groups']==s);hist[ix]=hist[rng.permutation(ix)]
        controls[key]=hist
    out={};predictions={'groups':sets['test']['groups'],'truth':sets['test']['y']}
    # Quadratic expansion of all ~100 history features would add excessive flexibility.
    # Primary A: quadratic present. Primary B: same A quadratic + linear history.
    # This design is fixed here before fitting; strength selected only on validation.
    transformed={}
    m=sets['train']['A'].mean(0);s=sets['train']['A'].std(0);s[s<1e-10]=1
    for key,d in sets.items():
        cur=expansion((d['A']-m)/s)
        h=d['B'];hrev=h.reshape(len(h),3,-1)[:,::-1,:].reshape(h.shape)
        transformed[key]={'current':cur,'history':np.column_stack([cur,h]),'shuffled_history':np.column_stack([cur,controls[key]]),'reversed_lags':np.column_stack([cur,hrev])}
    for col,target in enumerate(['future_deformation','disturbance_response']):
        out[target]={};scores={}
        for model in ['current','history','shuffled_history','reversed_lags','shuffled_targets']:
            family='history' if model=='shuffled_targets' else model
            yt=sets['train']['y'][:,col].copy()
            if model=='shuffled_targets':yt=rng.permutation(yt)
            pred,meta=fit_predict(transformed['train'][family],transformed['validation'][family],transformed['test'][family],yt,sets['validation']['y'][:,col],p['history']['alphas'],quadratic=False)
            scores[model]=rmse_by_group(sets['test']['y'][:,col],pred,sets['test']['groups'])
            out[target][model]={**meta,'mean_run_rmse':float(scores[model].mean()),'per_seed_rmse':scores[model].tolist()}
            predictions[target+'_'+model]=pred
            print(f'{target} / {model}: mean run RMSE {scores[model].mean():.6g}',flush=True)
        for model in ['frozen_rate','training_mean']:
            pred=sets['test']['frozen'][:,col] if model=='frozen_rate' else np.full(len(sets['test']['groups']),sets['train']['y'][:,col].mean())
            scores[model]=rmse_by_group(sets['test']['y'][:,col],pred,sets['test']['groups'])
            out[target][model]={'mean_run_rmse':float(scores[model].mean()),'per_seed_rmse':scores[model].tolist()}
        out[target]['paired_history_vs_current']=paired(scores['current'],scores['history'])
        out[target]['relative_rmse_reduction_percent']=float(100*(1-scores['history'].mean()/scores['current'].mean()))
    np.savez_compressed(a.output/'predictions.npz',**predictions)
    save_json(a.output/'prediction.json',out)
    # Numerical target/feature sensitivity uses same seed groups, never extra replicates.
    sensitivity={}
    for label in ['grid','dt']:
        ref=dataset([a.data/'refinement'/f'{label}_{s}.npz' for s in p['test_seeds']],p)
        sensitivity[label]={'target_rms_change':np.sqrt(np.mean((ref['y']-sets['test']['y'])**2,axis=0)).tolist()}
        cur=expansion((ref['A']-m)/s)
        for col,target in enumerate(['future_deformation','disturbance_response']):
            err=[]
            for model,te in [('current',cur),('history',np.column_stack([cur,ref['B']]))]:
                pred,_=fit_predict(transformed['train'][model],transformed['validation'][model],te,sets['train']['y'][:,col],sets['validation']['y'][:,col],p['history']['alphas'],False)
                err.append(rmse_by_group(ref['y'][:,col],pred,ref['groups']))
            sensitivity[label][target]={'current_rmse':float(err[0].mean()),'history_rmse':float(err[1].mean()),'paired':paired(*err)}
    save_json(a.output/'prediction_sensitivity.json',sensitivity)


if __name__=='__main__':main()
