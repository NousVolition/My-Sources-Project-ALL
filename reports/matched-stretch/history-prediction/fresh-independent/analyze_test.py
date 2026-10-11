"""Score frozen models on untouched runs; every reported interval resamples runs."""
from pathlib import Path
import json,argparse
import numpy as np
from fit_models import load_models,paths_for,HERE,make_dataset,make_matrices,save_json,sha
from evaluate import summarize
from fresh_statistics import spatial_scores,compare,classify


def score(d,case,models):
    Xs=make_matrices(d);y=d['target'];event=(y>=case['train_event_cutoff']).astype(int)
    results={};predictions={k:d[k] for k in ['target','group','anchor','particle','position']}
    predictions['event']=event
    for name,pair in models.items():
        key=name[5:] if name.startswith('tree_') else name
        pred=pair[0].predict(Xs[key]);prob=pair[1].predict_proba(Xs[key])[:,1]
        info=summarize(y,event,pred,prob,d['group'],case['models'][name]['decision_threshold'])
        info['spatial']=spatial_scores(event,prob,d['group'],d['anchor'])
        results[name]=info;predictions[name+'_prediction']=pred;predictions[name+'_probability']=prob
    return results,predictions


def evaluate_case(paths,p,case,models):
    d=make_dataset(paths,p,case['case']['options']);results,pred=score(d,case,models)
    paired={name+'_vs_A_full':compare(results['A_full'],result,p) for name,result in results.items() if name!='A_full'}
    if 'A_matched_current' in results:paired['B_vs_matched_current']=compare(results['A_matched_current'],results['B_history'],p)
    if 'tree_A_full' in results:paired['tree_B_vs_tree_A']=compare(results['tree_A_full'],results['tree_B_history'],p)
    if 'A_weak' in results:paired['B_weak_vs_A_weak']=compare(results['A_weak'],results['B_weak'],p)
    return {'models':results,'paired':paired,'n_rows':len(d['target']),
            'n_independent_runs':len(np.unique(d['group'])),'cutoff':case['train_event_cutoff']},pred


def analyze():
    p=json.loads((HERE/'protocol.json').read_text());lock=json.loads((HERE/'test-lock.json').read_text())
    manifest_path=HERE/lock['manifest']
    if sha(manifest_path)!=lock['manifest_sha256']:raise ValueError('Test lock manifest changed')
    m,models=load_models(manifest_path);paths=paths_for(p['test_seeds'],m['grid'],p)
    result={'grid':m['grid'],'manifest_sha256':sha(manifest_path),'test_lock_sha256':sha(HERE/'test-lock.json'),
        'protocol_sha256':sha(HERE/'protocol.json'),'test_files':{x.name:sha(x) for x in paths},'cases':{},'viscosity':{}}
    for name,case in m['cases'].items():
        r,pred=evaluate_case(paths,p,case,models[name]);result['cases'][name]=r
        np.savez_compressed(HERE/(name+'-predictions.npz'),**pred)
    for nu in p['viscosity_controls']['nu']:
        vp=paths_for(p['viscosity_controls']['seeds'],m['grid'],p,nu=nu)
        subset={k:models['primary'][k] for k in ['A_full','B_history']}
        r,_=evaluate_case(vp,p,m['cases']['primary'],subset)
        r['dependent_variants_of_test_groups']=True;result['viscosity'][str(nu)]=r
    audit=json.loads((HERE/'test-numerical-audit.json').read_text())
    passed=audit['all_screens_pass'];primary=result['cases']['primary']['paired']['B_history_vs_A_full']
    result['numerically_audited']=passed
    result['primary_conclusion']=classify(primary['spatial_AP']['ci95'],p['primary_endpoint']['useful_absolute_AP_gain'],passed)
    result['secondary_continuous_conclusion']=classify(primary['relative_RMSE_reduction']['ci95'],p['secondary_continuous']['useful_gain'],passed)
    save_json(HERE/'results.json',result)
    print(result['primary_conclusion'],flush=True)


if __name__=='__main__':analyze()
