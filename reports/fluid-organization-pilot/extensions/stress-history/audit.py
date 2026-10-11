"""Audit complete delivered bytes, freeze, numerical checks and every held-out forecast."""
import csv,json
from pathlib import Path
import numpy as np
from common import ROOT,digest,read_protocol
from forecast import verify_freeze
from features import inputs,target
from models import predict
def main():
    manifest=json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'))
    for f in manifest['files']:
        q=ROOT/f['path'];assert q.stat().st_size==f['size'];assert digest(q)==f['sha256'],f['path']
    p=read_protocol();m=verify_freeze()
    with (ROOT/'results/predictions.csv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    for regime in p['regimes']:
        rr=[r for r in rows if r['regime']==regime]
        for i,seed in enumerate(p['splits']['test']):
            with np.load(ROOT/'data/runs'/f'{regime}_base_{seed}.npz') as z:
                _,a,h=inputs(z,p['history_lags']);y=target(z)
            predicted=[predict(m['models'][regime]['present'],a),predict(m['models'][regime]['history'],np.c_[a,h])]
            for j in range(2):
                r=rr[2*i+j];assert int(r['seed'])==seed;assert abs(float(r['observed'])-y[j])<1e-12
                for name,v in zip(['present','history'],predicted):assert abs(float(r[name])-v[j])<1e-10
    print(json.dumps({'manifest_files':len(manifest['files']),'heldout_prediction_rows':len(rows),'frozen_models_verified':True,'numerical_recording_files':len(list((ROOT/'data').rglob('*.npz')))}))
if __name__=='__main__':main()
