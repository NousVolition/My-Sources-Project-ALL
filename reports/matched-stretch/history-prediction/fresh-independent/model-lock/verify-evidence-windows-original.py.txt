"""Reproduce the frozen model lock; never simulate or overwrite study evidence.

Default checks saved models/scalers, training cutoffs, validation diagnostics and
the numerical audit. --refit independently repeats selection in memory only.
--pretest additionally requires all untouched test data and claims to be absent.
"""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import argparse, datetime, hashlib, json, sys, warnings
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_curve

HERE = Path(__file__).resolve().parent
F = HERE.parent
sys.path.insert(0, str(F))
import record_history as rh
import fit_models as fm
import numeric_audit as na
from evaluate import macro_loss, fit_pair, summarize
from fresh_statistics import spatial_scores
from features import shuffle_blocks


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(refit=False, pretest=False, exact=False):
    differences = []
    def same(a, b):
        if isinstance(a, dict):
            assert a.keys() == b.keys()
            for k in a: same(a[k], b[k])
        elif isinstance(a, list):
            assert len(a) == len(b)
            for x, y in zip(a, b): same(x, y)
        elif isinstance(a, (float, np.floating)):
            assert np.isfinite(a) and np.isfinite(b)
            differences.append(float(abs(a-b)))
            assert a == b if exact else abs(a-b) <= 1e-12
        else:
            assert a == b

    p = read(F/'protocol.json')
    if pretest: fm.assert_test_absent(p)
    lock = read(F/'test-lock.json')
    mp = F/lock['manifest']
    audit_path = F/'model-audit-112.json'
    assert lock['grid'] == 112 and lock['test_absent_at_lock']
    assert sha(mp) == lock['manifest_sha256']
    assert sha(audit_path) == lock['model_numeric_audit_sha256']
    assert sha(F/'protocol.json') == lock['protocol_sha256']
    m, models = fm.load_models(mp)
    assert m['training_only_scalers'] and m['no_train_validation_refit'] and m['test_absent_at_freeze']
    assert m['sources'] == lock['sources']
    inventory = read(F/'pilot-complete/source-inventory.json')
    for group, base in [('fresh_runtime', F.parent.parent), ('fresh_implementation', F), ('predecessor', F.parent)]:
        for rel, digest in inventory[group].items():
            assert rh.source_hash_matches((base/rel).read_bytes(), digest), rel
    paths = fm.paths_for(p['train_seeds']+p['validation_seeds'], 112, p)
    assert len(paths) == 48 and set(m['files']) == {x.name for x in paths}
    archives = []
    for path in paths:
        meta = rh.validate(path, {**p['base'], 'n':112})
        assert meta['study_sources'] == m['sources']
        assert meta['completed_UTC'] < m['models_frozen_at_UTC'] < lock['locked_at_UTC']
        archives.append({'path':str(path.relative_to(F)).replace('\\','/'), 'sha256':sha(path),
                         'metadata_sha256':sha(path.with_suffix('.json'))})
    assert len(list((F/'models-112').glob('*.joblib'))) == 39
    expected_cases = fm.cases(p)
    assert set(m['cases']) == {x['name'] for x in expected_cases}
    cases_checked = []
    for case in expected_cases:
        entry = m['cases'][case['name']]
        same(entry['case'], case)
        assert set(entry['models']) == set(fm.model_names(case['name']=='primary'))
        d = fm.make_dataset(paths, p, case['options'])
        tr = np.isin(d['group'],p['train_seeds']); va = np.isin(d['group'],p['validation_seeds'])
        assert np.all(tr|va) and not np.any(tr&va)
        assert int(tr.sum()) == entry['n_train'] and int(va.sum()) == entry['n_validation']
        cutoff = float(np.quantile(d['target'][tr],case['quantile']))
        same(cutoff, entry['train_event_cutoff'])
        event = (d['target'] >= cutoff).astype(int)
        matrices = fm.make_matrices(d)
        for name, info in entry['models'].items():
            tree = name.startswith('tree_'); key = name[5:] if tree else name
            X = matrices[key]; pair = models[case['name']][name]
            assert X.shape[1] == info['input_dimensions'] and len(pair) == 2
            predictions = []
            for classification, model, kind in zip([False,True],pair,['continuous','event']):
                prediction = model.predict_proba(X[va])[:,1] if classification else model.predict(X[va])
                assert np.isfinite(prediction).all()
                predictions.append(prediction)
                target = event if classification else d['target']
                loss = macro_loss(target[va],prediction,d['group'][va],classification)
                same(loss,info[kind]['validation_loss'])
                candidates = info[kind]['candidates']
                expected = [1.,10.] if tree else p['models']['logistic_C' if classification else 'ridge_alphas']
                same([x['regularization'] for x in candidates],expected)
                winner = min(candidates,key=lambda x:x['validation_loss'])
                same(winner['regularization'],info[kind]['selected_regularization'])
                same(winner['validation_loss'],loss)
                if not tree:
                    scaler = model.steps[0][1]; recomputed = StandardScaler().fit(X[tr])
                    for attr in ['mean_','var_','scale_']:
                        aa=getattr(scaler,attr);bb=getattr(recomputed,attr)
                        same(aa.tolist(),bb.tolist())
                    assert int(scaler.n_samples_seen_) == int(tr.sum())
                    est = model.steps[-1][1]
                    same(est.coef_.reshape(-1).tolist(),info[kind]['standardized_coefficients'])
                    same(getattr(est,'C' if classification else 'alpha'),info[kind]['selected_regularization'])
                else:
                    assert model.max_iter==80 and model.max_leaf_nodes==7 and model.max_depth==3
                    assert model.learning_rate==.08 and model.min_samples_leaf==40
                    assert model.early_stopping is False and model.random_state==26
                    same(model.l2_regularization,info[kind]['selected_regularization'])
            pred,prob = predictions
            precision,recall,thresholds = precision_recall_curve(event[va],prob)
            f1=2*precision[:-1]*recall[:-1]/np.maximum(precision[:-1]+recall[:-1],1e-12)
            same(float(thresholds[np.argmax(f1)]),info['decision_threshold'])
            validation=summarize(d['target'][va],event[va],pred,prob,d['group'][va],info['decision_threshold'])
            validation['spatial']=spatial_scores(event[va],prob,d['group'][va],d['anchor'][va])
            same(validation,info['validation'])
            if refit:
                perm=None
                if name=='B_permuted_labels':
                    perm={False:shuffle_blocks(d['target'],d['group'],d['anchor'],61),True:shuffle_blocks(event,d['group'],d['anchor'],61)}
                repeated,repeated_models=fit_pair(X,d['target'],event,[tr,va,va],d['group'],family='tree' if tree else 'linear',permuted=perm)
                for kind in ['continuous','event','decision_threshold']:same(repeated[kind],info[kind])
                same(repeated['prediction'].tolist(),pred.tolist());same(repeated['probability'].tolist(),prob.tolist())
            print(case['name']+' '+name+' verified'+(' with independent in-memory refit' if refit else ''),flush=True)
        cases_checked.append({'case':case['name'],'models':len(entry['models']),'train_event_cutoff':cutoff,
                              'train_rows':entry['n_train'],'validation_rows':entry['n_validation']})
    # Run the unchanged audit with its sole output redirected into memory.
    # Neither the saved audit nor any numerical or fitted source is modified.
    captured=[];writer=na.save_json
    try:
        na.save_json=lambda path,value:captured.append((str(path),value))
        rebuilt=na.audit(112,'model')
    finally:
        na.save_json=writer
    assert len(captured)==1 and captured[0][1]==rebuilt
    same(rebuilt,read(audit_path))
    for s in rebuilt['screens']:
        metric=s['metric'];temporal='temporal' in s['comparison'] or 'half-step' in s['comparison']
        if metric=='target_relative_RMS': key='half_step_target_relative_RMS_max' if temporal else 'target_relative_RMS_max'
        elif metric=='event_Jaccard': key='event_Jaccard_min'
        elif metric.endswith('_new_probability_relative_RMS'):key='new_model_probability_relative_RMS_max'
        elif metric=='spatial_AP_contrast_drift':key='new_model_spatial_AP_contrast_'+('macro' if s['seed']=='macro' else 'per_run')+'_drift_max'
        elif metric=='RMSE_contrast_drift_relative_validation_A':key='new_model_RMSE_contrast_relative_validation_A_max'
        elif metric.startswith(('A_full_','B_history_')):key='frozen_old_prediction_relative_RMS_max'
        else:key=metric+'_max'
        bound=p['numeric_screens'][key];direction='min' if metric=='event_Jaccard' else 'max'
        passed=s['value']>=bound if direction=='min' else s['value']<=bound
        assert s['bound']==bound and s['direction']==direction and s['pass']==passed and passed
    assert rebuilt['all_screens_pass']
    if pretest:fm.assert_test_absent(p)
    assert sha(mp)==lock['manifest_sha256'] and sha(audit_path)==lock['model_numeric_audit_sha256']
    result={'verified_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'all_checks_pass':True,'grid':112,'train_archives':32,'validation_archives':16,
            'cases':cases_checked,'model_pairs':39,'fitted_estimators':78,
            'in_memory_refit_verified':refit,'test_absence_checked':pretest,
            'reconstruction_max_absolute_difference':max(differences,default=0.),
            'reconstruction_tolerance':0. if exact else 1e-12,
            'numerical_screens_passed':sum(x['pass'] for x in rebuilt['screens']),
            'numerical_screens_total':len(rebuilt['screens']),
            'protocol_sha256':sha(F/'protocol.json'),'manifest_sha256':sha(mp),
            'model_audit_sha256':sha(audit_path),'test_lock_sha256':sha(F/'test-lock.json'),
            'verification_code_sha256':sha(Path(__file__)),'archives':archives,
            'sources':m['sources'],'model_files':{v['file']:v['sha256'] for c in m['cases'].values() for v in c['models'].values()},
            'scope':'Numerical accuracy and model-lock integrity only. No untouched-test estimate or confidence interval.'}
    print(json.dumps({k:v for k,v in result.items() if k not in ['archives','sources','model_files','cases']}),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--refit',action='store_true');ap.add_argument('--pretest',action='store_true');ap.add_argument('--exact',action='store_true');ap.add_argument('--receipt',type=Path)
    args=ap.parse_args()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='No positive class found in y_true')
        result=verify(args.refit,args.pretest,args.exact)
    if args.receipt:
        with args.receipt.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2);stream.write('\n')
