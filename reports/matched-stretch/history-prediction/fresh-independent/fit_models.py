"""Fit only new train/validation runs, then seal every model before test generation."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,argparse,datetime
import numpy as np
import joblib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from evaluate import fit_pair, summarize
from features import shuffle_blocks
from simulate import sha,save_json
from fresh_features import make_dataset,make_matrices,cases
from fresh_statistics import spatial_scores
from record_history import name_for,validate,sources,source_hash_matches


def assert_test_absent(p):
    for seed in p['test_seeds']:
        for folder in [HERE/'recorded-data',HERE/'claims',HERE/'recovery']:
            if folder.exists() and list(folder.glob(f's{seed}_*')):
                raise ValueError('Test simulation already generated/claimed: fitting forbidden')


def paths_for(seeds,n,p,dt=None,nu=None):
    config={**p['base'],'n':n,'dt':dt or p['base']['dt'],'nu':nu or p['base']['nu']}
    paths=[HERE/'recorded-data'/(name_for(seed,config)+'.npz') for seed in seeds]
    for path in paths:validate(path,config)
    return paths


def model_names(primary):
    if not primary:return ['A_full','B_history']
    return ['A_full','B_history','A_matched_current','A_weak','B_weak',
        'B_past_only','B_no_past','B_no_exchange','B_reversed','B_irrelevant',
        'B_shuffled_91','B_shuffled_92','B_shuffled_93','B_permuted_labels',
        'B_mirrored_history','B_mirrored_reversed','tree_A_full','tree_B_history','tree_A_matched_current']


def fit_case(d,p,case,folder):
    tr=np.isin(d['group'],p['train_seeds']);va=np.isin(d['group'],p['validation_seeds'])
    if not tr.any() or not va.any() or not np.all(tr|va) or np.any(tr&va):raise ValueError('Invalid whole-run fitting split')
    cutoff=float(np.quantile(d['target'][tr],case['quantile']))
    event=(d['target']>=cutoff).astype(int);Xs=make_matrices(d)
    result={'train_event_cutoff':cutoff,'case':case,'n_train':int(tr.sum()),'n_validation':int(va.sum()),'models':{}}
    for name in model_names(case['name']=='primary'):
        tree=name.startswith('tree_');key=name[5:] if tree else name
        perm=None
        if name=='B_permuted_labels':
            perm={False:shuffle_blocks(d['target'],d['group'],d['anchor'],61),
                  True:shuffle_blocks(event,d['group'],d['anchor'],61)}
        # The third mask is validation again, never an independent test mask.
        # fit_pair's returned evaluation arrays are validation diagnostics only.
        info,models=fit_pair(Xs[key],d['target'],event,[tr,va,va],d['group'],
                            family='tree' if tree else 'linear',permuted=perm)
        pred=info.pop('prediction');prob=info.pop('probability');info.pop('validation_prediction')
        info['validation']=summarize(d['target'][va],event[va],pred,prob,d['group'][va],info['decision_threshold'])
        info['validation']['spatial']=spatial_scores(event[va],prob,d['group'][va],d['anchor'][va])
        info['input_dimensions']=Xs[key].shape[1]
        path=folder/(case['name']+'__'+name+'.joblib');joblib.dump(models,path,compress=3)
        info['file']=path.name;info['sha256']=sha(path)
        result['models'][name]=info
        print(case['name']+' '+name+' fitted on train/validation only',flush=True)
    return result


def fit(n):
    p=json.loads((HERE/'protocol.json').read_text());assert_test_absent(p)
    frozen=sources();folder=HERE/f'models-{n}';folder.mkdir(exist_ok=True)
    if (folder/'manifest.json').exists():raise ValueError('Existing fitted manifest: do not refit silently')
    paths=paths_for(p['train_seeds']+p['validation_seeds'],n,p)
    manifest={'grid':n,'protocol_sha256':sha(HERE/'protocol.json'),'sources':frozen,
        'fit_code_sha256':sha(HERE/'fit_models.py'),'files':{x.name:sha(x) for x in paths},
        'training_only_scalers':True,'no_train_validation_refit':True,'cases':{}}
    for case in cases(p):
        assert_test_absent(p)
        d=make_dataset(paths,p,case['options']);manifest['cases'][case['name']]=fit_case(d,p,case,folder)
        save_json(folder/'fitting-progress.json',manifest)
    if sources()!=frozen:raise ValueError('Numerical/features source changed during fitting')
    assert_test_absent(p)
    manifest['models_frozen_at_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    manifest['test_absent_at_freeze']=True
    save_json(folder/'manifest.json',manifest)
    return folder/'manifest.json'


def load_models(manifest_path):
    manifest_path=Path(manifest_path);m=json.loads(manifest_path.read_text())
    if m['protocol_sha256']!=sha(HERE/'protocol.json') or set(m['sources'])!=set(sources()):raise ValueError('Frozen model/protocol source changed')
    for rel,digest in m['sources'].items():
        if not source_hash_matches((HERE.parent.parent/rel).read_bytes(),digest):raise ValueError('Frozen source changed: '+rel)
    if not source_hash_matches((HERE/'fit_models.py').read_bytes(),m['fit_code_sha256']):raise ValueError('Fitting source changed')
    for name,digest in m['files'].items():
        if sha(HERE/'recorded-data'/name)!=digest:raise ValueError('Fitting archive changed: '+name)
    out={}
    for case,entry in m['cases'].items():
        out[case]={}
        for name,model in entry['models'].items():
            path=manifest_path.parent/model['file']
            if sha(path)!=model['sha256']:raise ValueError('Frozen model checksum changed')
            # Only this study's locally generated, checksum-verified model files.
            out[case][name]=joblib.load(path)
    return m,out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--grid',type=int,required=True)
    print(fit(ap.parse_args().grid))
