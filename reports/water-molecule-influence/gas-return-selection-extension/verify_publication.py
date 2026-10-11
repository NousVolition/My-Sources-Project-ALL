"""Verify the published snapshot; --numeric independently rebuilds forecasts."""
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=json.loads((ROOT/'manifest_sha256.json').read_text())
for name,digest in manifest.items():assert sha(ROOT/name)==digest,name
r=json.loads((ROOT/'data/results.json').read_text())
source=ROOT.parent/'gas-damping-extension/data'
for name,digest in r['input_sha256'].items():assert sha(source/name)==digest,name
lock=json.loads((ROOT/'data/locked_training.json').read_text())
assert r['validation']['passed'] and lock==r['training']
assert lock['script_sha256']==sha(ROOT/'analyze_selection.py')
assert lock['protocol_sha256']==sha(ROOT/'data/protocol.json')
for method in ['direction','velocity']:
    best=min([x for x in lock['grid'] if x['method']==method],key=lambda x:(x['mean_mse'],x['weight']))
    assert best['weight']==lock['weights'][method]
checks=0
if '--numeric' in sys.argv:
    import numpy as np
    methods=['position','direction','velocity','hold','linear'];H=np.array([.02,.05,.1])
    for rep in range(3):
        for g in [.2,1,10]:
            for fine in ([False,True] if rep==1 else [False]):
                suffix='_fine' if fine else '';dt=.001 if fine else .01
                with np.load(source/f'vapor_r{rep}_gamma{g:g}{suffix}.npz') as z:p=z['state'][:,:2]
                with np.load(ROOT/f'data/predictions_r{rep}_gamma{g:g}{suffix}.npz') as d:
                    a=d['position__time'];ids=np.rint(a/dt).astype(int);current=p[ids]
                    truth=p[ids[:,None]+np.rint(H/dt).astype(int)]
                    for method in methods:
                        assert np.array_equal(a,d[method+'__time'])
                        if method in methods[:3]:
                            b=d[method+'__old_time'];assert np.all(b+.1<=a+1e-10)
                            assert np.all(a-b<=1+1e-10)
                            times=np.arange(len(p))*dt
                            old=np.column_stack([np.interp(b,times,p[:,j]) for j in range(2)])
                            continuation=np.stack([np.column_stack([np.interp(b+h,times,p[:,j]) for j in range(2)]) for h in H],axis=1)
                            prediction=current[:,None]+continuation-old[:,None]
                            assert np.all(d[method+'__initial_velocity_defined'])
                            assert np.max(d[method+'__initial_distance'])<=.02+1e-10
                        elif method=='hold':prediction=np.repeat(current[:,None,:],3,axis=1)
                        else:prediction=current[:,None]+H[None,:,None]*(p[ids]-p[ids-1])[:,None]/dt
                        error=np.sum((prediction-truth)**2,axis=-1)
                        assert np.allclose(prediction,d[method+'__forecast'],atol=1e-10)
                        assert np.allclose(error,d[method+'__error2'],atol=1e-10)
                        if fine:summary=next(x for x in r['fine_recording_check'] if x['friction_per_ps']==g)['native'][method]
                        else:summary=next(x for x in r['rows'] if x['start']==rep+1 and x['friction_per_ps']==g)['methods'][method]
                        assert summary['n']==len(a)
                        assert np.allclose(summary['mse'],error.mean(axis=0),atol=1e-10)
                    checks+=1
print(json.dumps(dict(passed=True,published_hashes=len(manifest),input_hashes=len(r['input_sha256']),independent_reconstructions=checks)))
