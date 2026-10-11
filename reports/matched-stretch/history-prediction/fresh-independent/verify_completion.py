"""Independent archive, split, frozen-model and saved-prediction audit."""
import json
import numpy as np
from fit_models import HERE,load_models,paths_for,sha,save_json
from fresh_features import make_dataset,make_matrices
from analyze_test import evaluate_case


def verify():
    p=json.loads((HERE/'protocol.json').read_text());lock=json.loads((HERE/'test-lock.json').read_text())
    manifest=HERE/lock['manifest'];assert sha(manifest)==lock['manifest_sha256']
    m,models=load_models(manifest);r=json.loads((HERE/'results.json').read_text())
    sets=[set(p[k]) for k in ['train_seeds','validation_seeds','test_seeds']]+[set(p['pilot']['seeds'])]
    assert all(not sets[i]&sets[j] for i in range(4) for j in range(i+1,4))
    assert m['test_absent_at_freeze'] and lock['test_absent_at_lock']
    paths=paths_for(p['test_seeds'],m['grid'],p);errors={}
    for path in paths:
        meta=json.loads(path.with_suffix('.json').read_text());assert meta['started_UTC']>lock['locked_at_UTC']
    for name,case in m['cases'].items():
        d=make_dataset(paths,p,case['case']['options']);xs=make_matrices(d)
        with np.load(HERE/(name+'-predictions.npz'),allow_pickle=False) as z:
            for key in ['target','group','anchor','particle','position']:np.testing.assert_array_equal(z[key],d[key])
            for model,pair in models[name].items():
                key=model[5:] if model.startswith('tree_') else model
                for kind,pred in [('prediction',pair[0].predict(xs[key])),('probability',pair[1].predict_proba(xs[key])[:,1])]:
                    err=float(abs(pred-z[model+'_'+kind]).max());assert err<1e-12;errors[name+'_'+model+'_'+kind]=err
        rerun,_=evaluate_case(paths,p,case,models[name])
        # Also independently reproduce metrics and paired uncertainty, not only scores.
        assert json.dumps(rerun,sort_keys=True)==json.dumps(r['cases'][name],sort_keys=True)
    inventory=json.loads((HERE/'execution.json').read_text())
    names=[x['name'] for x in inventory['records']];assert len(names)==len(set(names))<=p['max_new_simulations']
    for row in inventory['records']:
        path=HERE/'recorded-data'/(row['name']+'.npz')
        from record_history import validate
        meta=validate(path,row['config']);assert meta['sha256']==row['sha256']
    assert r['test_files']=={x.name:sha(x) for x in paths}
    save_json(HERE/'verification.json',{'verified':True,'independent_test_runs':len(paths),
        'protocol_sha256':sha(HERE/'protocol.json'),'manifest_sha256':sha(manifest),
        'results_sha256':sha(HERE/'results.json'),'checked_archives':len(names),
        'max_prediction_reproduction_error':max(errors.values()),'prediction_errors':errors,
        'all_metrics_and_paired_intervals_reproduced':True})


if __name__=='__main__':verify()
