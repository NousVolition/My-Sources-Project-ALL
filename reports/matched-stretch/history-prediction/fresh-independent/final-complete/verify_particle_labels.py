"""Check consistent relabeling on the first preregistered test run; never simulate or refit."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,tempfile,hashlib,argparse
import numpy as np
HERE=Path(__file__).resolve().parent;F=HERE.parent;sys.path.insert(0,str(F))
import fit_models as fm
def verify():
    p=json.loads((F/'protocol.json').read_text());lock=json.loads((F/'test-lock.json').read_text())
    m,models=fm.load_models(F/lock['manifest']);path=fm.paths_for([p['test_seeds'][0]],112,p)[0]
    d=fm.make_dataset([path],p);perm=np.random.default_rng(70).permutation(p['base']['markers'])
    order=np.concatenate([i*len(perm)+perm for i in range(len(p['anchors']))])
    with tempfile.TemporaryDirectory(prefix='fresh-id-audit-') as folder:
        target=Path(folder)/path.name
        with np.load(path,allow_pickle=False) as z:
            data={k:z[k] if k=='time' else z[k][:,perm] for k in z.files}
        np.savez_compressed(target,**data);target.with_suffix('.json').write_text(json.dumps({'seed':p['test_seeds'][0]}))
        changed=fm.make_dataset([target],p)
    errors={}
    for key in ['present','history','reverse','irrelevant','mirrored','mirrored_reverse','target']:
        np.testing.assert_allclose(d[key][order],changed[key],atol=2e-8,rtol=2e-8)
        errors[key]=float(np.max(abs(d[key][order]-changed[key])))
    a=fm.make_matrices(d);b=fm.make_matrices(changed)
    for name in ['A_full','B_history']:
        for kind,model in zip(['continuous','event'],models['primary'][name]):
            x=model.predict(a[name]) if kind=='continuous' else model.predict_proba(a[name])[:,1]
            y=model.predict(b[name]) if kind=='continuous' else model.predict_proba(b[name])[:,1]
            np.testing.assert_allclose(x[order],y,atol=2e-8,rtol=2e-8)
            errors[name+'_'+kind]=float(np.max(abs(x[order]-y)))
    return {'all_checks_pass':True,'seed':p['test_seeds'][0],'permutation_seed':70,'markers':512,'anchors':5,
            'tolerance':'atol=rtol=2e-8, unchanged from the frozen relabeling fixture test',
            'no_particle_ID_supplied_as_predictor':True,'maximum_absolute_errors':errors,
            'source_archive_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();r=verify()
    if args.output:
        with args.output.open('x',encoding='utf-8') as s:json.dump(r,s,indent=2);s.write('\n')
    print(json.dumps(r),flush=True)
