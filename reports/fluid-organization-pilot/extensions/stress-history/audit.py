"""Audit complete delivered bytes, freeze, numerical checks and every held-out forecast."""
import csv,json
from pathlib import Path
import numpy as np
from common import ROOT,digest,read_protocol
from forecast import verify_freeze,load,families,null_errors,bootstrap,bootstrap_null
from features import inputs,target
from models import predict
def main():
    manifest=json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'))
    for f in manifest['files']:
        q=ROOT/f['path'];assert q.stat().st_size==f['size'];assert digest(q)==f['sha256'],f['path']
    p=read_protocol();m=verify_freeze()
    with (ROOT/'results/predictions.csv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    summary=json.loads((ROOT/'results/prediction_summary.json').read_text(encoding='utf-8'))
    checks=json.loads((ROOT/'results/numerical_verification.json').read_text(encoding='utf-8'))
    assert checks['all_numerical_passed'] and checks['all_refinements_passed']
    def equal(actual,expected):
        if isinstance(actual,dict):
            assert actual.keys()==expected.keys()
            for key in actual:equal(actual[key],expected[key])
        elif isinstance(actual,(list,float,int)) and not isinstance(actual,bool):np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-12)
        else:assert actual==expected
    with np.load(ROOT/'data/forecast/shuffled_errors.npz') as saved:
      for regime_index,regime in enumerate(p['regimes']):
        rr=[r for r in rows if r['regime']==regime]
        c,a,h,y,hashes=load(regime,'test');x=families(c,a,h,regime,'test')
        predicted={name:predict(model,x[name]) for name,model in m['models'][regime].items()};predicted['prior_by_strength']=np.tile(m['prior_by_strength'][regime],len(y)//2)
        assert len(rr)==len(y)
        for j,r in enumerate(rr):
            assert int(r['seed'])==p['splits']['test'][j//2];assert float(r['relative_strength'])==p['kick_rms_over_baseline_rms'][j%2]
            assert abs(float(r['observed'])-y[j])<1e-12
            for name,v in predicted.items():assert abs(float(r[name])-v[j])<1e-10
        errors={name:np.sqrt(((v-y)**2).reshape(-1,2).mean(1)) for name,v in predicted.items()}
        matrix=null_errors(m['models'][regime],a,h,y);np.testing.assert_allclose(matrix,saved[regime],rtol=1e-10,atol=1e-12)
        contrasts={'history_vs_present':bootstrap(errors['history'],errors['present'],p,1200+regime_index),'history_vs_shuffled_mean':bootstrap_null(errors['history'],matrix.mean(0),p,1300+regime_index),'present_vs_conventional':bootstrap(errors['present'],errors['conventional'],p,1400+regime_index),'history_vs_training_prior':bootstrap(errors['history'],errors['prior_by_strength'],p,1500+regime_index)}
        equal(contrasts,summary[regime]['contrasts']);equal(hashes,summary[regime]['test_data_sha256']);equal({k:float(v.mean()) for k,v in errors.items()},summary[regime]['mean_field_rmse'])
    print(json.dumps({'manifest_files':len(manifest['files']),'heldout_prediction_rows':len(rows),'frozen_models_verified':True,'numerical_recording_files':len(list((ROOT/'data').rglob('*.npz')))}))
if __name__=='__main__':main()
