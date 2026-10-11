"""Larger prespecified independent fluid confirmation with frozen predictors."""
import argparse,json,hashlib,os,sys,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from support import ROOT,save,table
from prediction import predict,shuffled_history
C=json.loads((ROOT/'fluid_confirmation_protocol.json').read_text())
def filehash(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def freeze():
    names=C['frozen_models']+['fluid/core.py','fluid/recovery.py','fluid/vendor/numerics.py','fluid/protocol.json','fluid_confirmation_protocol.json']
    hashes={p:filehash(ROOT/p) for p in names};path=ROOT/'data/confirmation_frozen_hashes.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['files']!=hashes:raise ValueError('Frozen sources/models changed; preserve prior confirmation before altering them')
    else:save('data/confirmation_frozen_hashes.json',{'declared_before_new_fields':True,'files':hashes})
    return hashes
def simulate(workers):
    freeze();jobs=[(s,'confirmation',None,None) for s in range(C['test_seed_start'],C['test_seed_stop_exclusive'])]
    for s in C['numerical_refinement_seeds']:jobs.extend([(s,'confirmation_half',None,C['fine_dt']),(s,'confirmation_fine',C['fine_n'],C['fine_dt'])])
    env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join(sys.path);env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
    def task(job):
        seed,label,n,dt=job;args=[sys.executable,str(ROOT/'run_fluid.py'),'--seed',str(seed),'--label',label]
        if n:args+=['--n',str(n)]
        if dt:args+=['--dt',str(dt)]
        result=subprocess.run(args,env=env,capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        return result.stdout.strip() or f'Cached {label} {seed}'
    start=time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(task,j) for j in jobs]
        for count,f in enumerate(as_completed(futures),1):
            print(f'{count}/{len(jobs)}: {f.result()}',flush=True)
    print(f'Completed independent confirmation simulations in {time.perf_counter()-start:.1f}s',flush=True)
def analyze():
    freeze();features={'current':[],'history':[],'groups':[]};targets=[];rows=[];metas=[];refinements=[]
    for seed in range(C['test_seed_start'],C['test_seed_stop_exclusive']):
        path=ROOT/f'fluid/data/fluid/confirmation_{seed}.npz';q=np.load(path);metas.append(json.loads(path.with_suffix('.json').read_text()))
        h=q['history'];ht=q['history_times'];t=q['time'];e0=q['energy'][0,0];ix=t<=8+1e-9
        integral=np.trapezoid(q['energy'][ix]/e0,t[ix],axis=0);hist=np.concatenate([(h[-1]-h[np.argmin(abs(ht+lag))])/lag for lag in (.5,1)])
        for a in (1,2):
            features['current'].append(np.r_[h[-1],a==1,a==2]);features['history'].append(hist);features['groups'].append(seed);targets.append(integral[a]-integral[0])
        if seed in C['numerical_refinement_seeds']:
            half=np.load(ROOT/f'fluid/data/fluid/confirmation_half_{seed}.npz');fine=np.load(ROOT/f'fluid/data/fluid/confirmation_fine_{seed}.npz')
            refinements.append({'seed':seed,'time_error_over_initial_energy':float(abs(q['energy']-half['energy']).max()/e0),'grid_error_over_initial_energy':float(abs(fine['energy']-half['energy']).max()/half['energy'][0,0])})
    d={k:np.array(v) for k,v in features.items()};y=np.array(targets);inputs={'present':d['current'],'history':np.c_[d['current'],d['history']],'shuffled_history':np.c_[d['current'],shuffled_history(d,np.random.default_rng(C['shuffle_seed']))]}
    forecasts={k:predict(np.load(ROOT/f'data/fluid_model_{k}.npz'),v) for k,v in inputs.items()};errors={k:np.sqrt(np.mean((v-y).reshape(-1,2)**2,axis=1)) for k,v in forecasts.items()}
    rng=np.random.default_rng(C['bootstrap_seed']);ix=rng.integers(0,C['independent_test_fields'],(C['bootstrap_resamples'],C['independent_test_fields']))
    def ci(x):return {'mean':float(x.mean()),'ci95':np.quantile(x[ix].mean(1),[.025,.975]).tolist(),'n':len(x)}
    delta=errors['history']-errors['present'];relative=1-errors['history'][ix].mean(1)/errors['present'][ix].mean(1)
    summary={k:ci(v) for k,v in errors.items()};summary['history_minus_present']=ci(delta);summary['history_minus_shuffled']=ci(errors['history']-errors['shuffled_history'])
    summary['relative_mean_rmse_reduction']={'value':float(1-errors['history'].mean()/errors['present'].mean()),'ci95':np.quantile(relative,[.025,.975]).tolist()}
    summary['median_paired_error_difference']=float(np.median(delta));summary['fields_improved']=int(np.sum(delta<0));summary['independent_fields']=len(delta)
    summary['primary_improvement_resolved']=summary['history_minus_present']['ci95'][1]<0
    summary['per_arm_mean_absolute_error']={k:dict(zip(['abrupt','gradual'],np.mean(abs(v-y).reshape(-1,2),axis=0).tolist())) for k,v in forecasts.items()}
    summary['numerical']={'max_energy_budget':max(m['energy_budget_residual'] for m in metas),'max_burn_budget':max(m['burn_budget_residual'] for m in metas),'max_divergence':max(m['divergence_max'] for m in metas),'max_cfl':max(m['max_cfl'] for m in metas),'max_clone_error':max(m['clone_error'] for m in metas),'refinements':refinements}
    v=summary['numerical'];g=C['gates'];summary['gates']={'energy':max(v['max_energy_budget'],v['max_burn_budget'])<g['normalized_energy_budget'],'divergence':v['max_divergence']<g['divergence'],'cfl':v['max_cfl']<g['max_cfl'],'clones':v['max_clone_error']==0,'refinements':all(max(r['time_error_over_initial_energy'],r['grid_error_over_initial_energy'])<g['time_or_grid_relative_energy_error'] for r in refinements)}
    summary['interpretation']=C['scope'];summary['model_hashes']=freeze();save('results/fluid_confirmation.json',summary)
    for j,seed in enumerate(range(C['test_seed_start'],C['test_seed_stop_exclusive'])):
        rows.append({'seed':seed,'true_abrupt':float(y[2*j]),'true_gradual':float(y[2*j+1]),**{k+'_rmse':float(v[j]) for k,v in errors.items()},'history_minus_present':float(delta[j])})
    table('fluid_confirmation.csv',rows);np.savez_compressed(ROOT/'data/fluid_confirmation.npz',current=d['current'],history_features=d['history'],groups=d['groups'],truth=y,**forecasts)
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--analyze',action='store_true');ap.add_argument('--workers',type=int,default=4);args=ap.parse_args()
    if args.analyze:analyze()
    else:simulate(args.workers);analyze()
