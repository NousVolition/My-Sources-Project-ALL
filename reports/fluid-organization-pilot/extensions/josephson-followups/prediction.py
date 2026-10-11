"""Frozen nonlinear predictors, run-wise holdout, and shuffled-history controls."""
import numpy as np
import json,hashlib
from scipy.linalg import solve
from threadpoolctl import threadpool_limits
from circuit import to_state,stable_velocity
from support import ROOT,P,save,table,interval,evolve

def squared_distance(a,b):
    return np.maximum((a*a).sum(1)[:,None]+(b*b).sum(1)[None,:]-2*a@b.T,0)
def fit_kernel(train,val,y,yv):
    mean=train.mean(0);scale=train.std(0);scale[scale<1e-10]=1
    a=(train-mean)/scale;b=(val-mean)/scale;ym=float(y.mean());d=squared_distance(a,a)
    median=float(np.median(d[np.triu_indices(len(a),1)]));median=max(median,1e-8);best=None
    for factor in P['prediction']['bandwidth_factors']:
        width=median*factor**2;k=np.exp(-d/(2*width));kv=np.exp(-squared_distance(b,a)/(2*width))
        for reg in P['prediction']['regularization']:
            weights=solve(k+reg*np.eye(len(k)),y-ym,assume_a='pos')
            loss=float(np.mean((kv@weights+ym-yv)**2))
            if best is None or loss<best[0]:best=(loss,reg,width,weights)
    loss,reg,width,weights=best
    return {'mean':mean,'scale':scale,'train_scaled':a,'target_mean':np.array(ym),'width':np.array(width),'weights':weights,'regularization':np.array(reg),'validation_rmse':np.array(np.sqrt(loss))}
def predict(model,x):
    b=(x-model['mean'])/model['scale'];return np.exp(-squared_distance(b,model['train_scaled'])/(2*model['width']))@model['weights']+model['target_mean']
def shuffled_history(d,rng):
    """Permute complete runs within prior-bias strata, keeping observation/arm slots."""
    h=d['history'];out=h.copy();groups=np.unique(d['groups'])
    strata=d.get('stratum',np.zeros(len(h)))
    for stratum in np.unique(strata):
        gs=np.array([g for g in groups if strata[np.flatnonzero(d['groups']==g)[0]]==stratum]);order=rng.permutation(gs)
        for g,other in zip(order,np.roll(order,1)):
            out[d['groups']==g]=h[d['groups']==other]
    return out
def compare(parts,name):
    rng=np.random.default_rng(43120);transformed={}
    for split,d in parts.items():
        present=d['current'];simple=d.get('simple',present)
        transformed[split]={'simple_present':simple,'present':present,'history':np.c_[present,d['history']],'shuffled_history':np.c_[present,shuffled_history(d,rng)]}
    scores={};predictions={};summary={};truth=parts['test']['y'];groups=parts['test']['groups']
    for modelname in ('simple_present','present','history','shuffled_history'):
        with threadpool_limits(limits=1):
            path=ROOT/f'data/{name}_model_{modelname}.npz';lock=ROOT/'data/confirmation_frozen_hashes.json'
            locked=json.loads(lock.read_text())['files'] if lock.exists() else {}
            key=path.relative_to(ROOT).as_posix()
            if key in locked:
                if hashlib.sha256(path.read_bytes()).hexdigest()!=locked[key]:raise ValueError('Frozen model changed: '+key)
                with np.load(path) as z:model={k:z[k] for k in z.files}
            else:
                model=fit_kernel(transformed['train'][modelname],transformed['validation'][modelname],parts['train']['y'],parts['validation']['y'])
                # Save a frozen model before evaluating any test outcomes.
                np.savez_compressed(path,**model)
            pred=predict(model,transformed['test'][modelname])
        rmse=np.array([np.sqrt(np.mean((pred[groups==g]-truth[groups==g])**2)) for g in np.unique(groups)])
        scores[modelname]=rmse;predictions[modelname]=pred
        summary[modelname]={**interval(rmse),'per_run_rmse':rmse.tolist(),'regularization':float(model['regularization']),'kernel_width':float(model['width']),'validation_rmse':float(model['validation_rmse']),'features':transformed['train'][modelname].shape[1]}
    for baseline in ('simple_present','present','shuffled_history'):
        summary['history_minus_'+baseline]=interval(scores['history']-scores[baseline])
    np.savez_compressed(ROOT/f'data/{name}_predictions.npz',truth=truth,groups=groups,**predictions)
    save(f'results/{name}_prediction.json',summary)
    return summary

def circuit_dataset():
    c=P['prediction'];rng=np.random.default_rng(c['seed']);n=c['train_starts']+c['validation_starts']+c['test_starts']
    phi=rng.uniform(-np.pi,np.pi,(n,2));state,sign=to_state(phi);bias=np.resize(c['prior_biases'],n)
    times=np.unique([t-lag for t in c['snapshots'] for lag in [0]+c['lags']]);path=evolve(state,bias,c['alpha'],max(times),times)
    currents=[];simple=[];history=[];groups=[];stratum=[];states=[];bs=[];changes=[]
    for g in range(n):
        for t in c['snapshots']:
            s=path[np.flatnonzero(times==t)[0],g];v=stable_velocity(s,sign[g],bias[g],c['alpha'])
            hist=[stable_velocity(path[np.flatnonzero(times==t-lag)[0],g],sign[g],bias[g],c['alpha']).sum() for lag in c['lags']]
            for change in c['interventions']:
                currents.append(np.r_[v,bias[g],bias[g]+change]);simple.append([v.sum(),bias[g],bias[g]+change]);history.append(hist)
                groups.append(g);stratum.append(bias[g]);states.append(s);bs.append(bias[g]);changes.append(change)
    states=np.array(states);bs=np.array(bs);changes=np.array(changes)
    end0=evolve(states,bs,c['alpha'],c['horizon'],[c['horizon']])[-1]
    end1=evolve(states,bs+changes,c['alpha'],c['horizon'],[c['horizon']])[-1]
    y=2*(end1[:,0]-end0[:,0])/c['horizon'];groups=np.array(groups)
    data={'current':np.array(currents),'simple':np.array(simple),'history':np.array(history),'groups':groups,'stratum':np.array(stratum),'y':y}
    parts={};masks={}
    for split,lo,hi in [('train',0,c['train_starts']),('validation',c['train_starts'],c['train_starts']+c['validation_starts']),('test',c['train_starts']+c['validation_starts'],n)]:
        mask=(groups>=lo)&(groups<hi);masks[split]=mask;parts[split]={k:v[mask] for k,v in data.items()}
    # Complete-state physics benchmark, using tighter settings than target generation.
    mask=masks['test'];a=evolve(states[mask],bs[mask],c['alpha'],c['horizon'],[c['horizon']],max_step=.1,rtol=2e-12,atol=2e-13)[-1]
    b=evolve(states[mask],(bs+changes)[mask],c['alpha'],c['horizon'],[c['horizon']],max_step=.1,rtol=2e-12,atol=2e-13)[-1]
    oracle=2*(b[:,0]-a[:,0])/c['horizon']
    np.savez_compressed(ROOT/'data/circuit_prediction_inputs.npz',initial_phase=phi,observation_times=times,path=path,branch_states=states,prior_bias=bs,interventions=changes,end_unchanged=end0,end_changed=end1,oracle=oracle,**data)
    summary=compare(parts,'circuit');summary['full_state_physics_max_error']=float(abs(oracle-y[mask]).max())
    summary['full_state_physics_rmse']=float(np.sqrt(np.mean((oracle-y[mask])**2)))
    summary['split_counts']={k:len(np.unique(d['groups'])) for k,d in parts.items()}
    save('results/circuit_prediction.json',summary);print('circuit history comparison',summary,flush=True)

def fluid_dataset():
    z=np.load(ROOT/'data/fluid_training.npz');parts={s:{k:z[s+'_'+k] for k in ('current','history','y','groups')} for s in ('train','validation')}
    rows=[];meta=[];test={'current':[],'history':[],'y':[],'groups':[]};checks=[]
    for seed in P['fluid_replication']['new_test_seeds']:
        path=ROOT/f'fluid/data/fluid/replication_{seed}.npz';q=np.load(path)
        import json
        meta.append(json.loads(path.with_suffix('.json').read_text()));h=q['history'];ht=q['history_times'];t=q['time'];e=q['energy']/q['energy'][0,0];ix=t<=8+1e-9
        integ=np.trapezoid(e[ix],t[ix],axis=0)
        hist=np.concatenate([(h[-1]-h[np.argmin(abs(ht+lag))])/lag for lag in (.5,1)])
        for arm in (1,2):
            test['current'].append(np.r_[h[-1],arm==1,arm==2]);test['history'].append(hist);test['y'].append(integ[arm]-integ[0]);test['groups'].append(seed)
            rows.append({'seed':seed,'arm':'abrupt' if arm==1 else 'gradual','integrated_normalized_energy_difference':float(integ[arm]-integ[0])})
        if seed in P['fluid_replication']['new_test_seeds'][:2]:
            half=np.load(ROOT/f'fluid/data/fluid/half_{seed}.npz');fine=np.load(ROOT/f'fluid/data/fluid/fine_{seed}.npz')
            checks.append({'seed':seed,'time_curve_error_over_initial':float(abs(q['energy']-half['energy']).max()/q['energy'][0,0]),'grid_curve_error_over_initial':float(abs(fine['energy']-half['energy']).max()/half['energy'][0,0])})
    parts['test']={k:np.array(v) for k,v in test.items()};summary=compare(parts,'fluid')
    summary['validation']={'max_energy_budget_error':max(d['energy_budget_residual'] for d in meta),'max_burn_budget_error':max(d['burn_budget_residual'] for d in meta),'max_divergence':max(d['divergence_max'] for d in meta),'max_cfl':max(d['max_cfl'] for d in meta),'max_clone_error':max(d['clone_error'] for d in meta),'refinement':checks}
    summary['gates']={'energy_budget':summary['validation']['max_energy_budget_error']<1e-5,'divergence':summary['validation']['max_divergence']<1e-8,'cfl':summary['validation']['max_cfl']<.5,'identical_clones':summary['validation']['max_clone_error']==0,'refinement':all(max(d['time_curve_error_over_initial'],d['grid_curve_error_over_initial'])<.01 for d in checks)}
    summary['split_counts']={k:len(np.unique(d['groups'])) for k,d in parts.items()}
    save('results/fluid_prediction.json',summary);table('fluid_runs.csv',rows)
    np.savez_compressed(ROOT/'data/fluid_test_features.npz',**parts['test']);print('fluid history comparison',summary,flush=True)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--fluid',action='store_true');args=ap.parse_args()
    fluid_dataset() if args.fluid else circuit_dataset()
