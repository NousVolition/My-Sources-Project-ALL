"""Archive original-platform validation probabilities without fitting or simulation."""
from pathlib import Path
import sys,json,hashlib,datetime,warnings,os
os.environ.setdefault("OMP_NUM_THREADS","1")
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
import numpy as np
HERE=Path(__file__).resolve().parent;F=HERE.parent
sys.path.insert(0,str(F))
import fit_models as fm
from evaluate import summarize
from fresh_statistics import spatial_scores

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=json.loads((F/'protocol.json').read_text());fm.assert_test_absent(p)
    mp=F/'models-112/manifest.json';m,models=fm.load_models(mp)
    paths=fm.paths_for(p['train_seeds']+p['validation_seeds'],112,p)
    out=HERE/'validation-reference';out.mkdir(exist_ok=True)
    assert not (HERE/'validation-reference.json').exists()
    primary=None;rows=[]
    for case in fm.cases(p):
        if not case['options'] and primary is not None:d=primary
        else:d=fm.make_dataset(paths,p,case['options'])
        if not case['options']:primary=d
        va=np.isin(d['group'],p['validation_seeds']);entry=m['cases'][case['name']]
        event=(d['target'][va]>=entry['train_event_cutoff']).astype(int);xs=fm.make_matrices(d)
        for name,info in entry['models'].items():
            key=name[5:] if name.startswith('tree_') else name;pair=models[case['name']][name]
            pred=pair[0].predict(xs[key][va]);prob=pair[1].predict_proba(xs[key][va])[:,1]
            value=summarize(d['target'][va],event,pred,prob,d['group'][va],info['decision_threshold'])
            value['spatial']=spatial_scores(event,prob,d['group'][va],d['anchor'][va])
            assert value==info['validation'],(case['name'],name)
            path=out/(case['name']+'__'+name+'.npz');assert not path.exists()
            np.savez_compressed(path,probability=prob)
            rows.append({'case':case['name'],'model':name,'path':path.relative_to(HERE).as_posix(),
                         'sha256':sha(path),'rows':len(prob),'model_sha256':info['sha256']})
        print(case['name']+' reference probabilities exactly reproduce saved validation metrics',flush=True)
    fm.assert_test_absent(p)
    result={'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest_sha256':sha(mp),
            'export_code_sha256':sha(Path(__file__)),'all_original_validation_metrics_exact':True,
            'purpose':'Original Windows probability bytes for independently bounding cross-platform rounding at locked decision thresholds; no fitting, simulation, criterion change or threshold change.',
            'files':rows}
    with (HERE/'validation-reference.json').open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')

if __name__=='__main__':
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='No positive class found in y_true')
        main()
