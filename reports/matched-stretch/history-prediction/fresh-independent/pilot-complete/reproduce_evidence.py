"""Verify the completed pilot from saved archives; no simulation or fitting."""
from pathlib import Path
import sys, json, hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
FRESH=HERE.parent
sys.path.insert(0,str(FRESH))
import record_history as rh
import numeric_audit as na

def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def same(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:same(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):same(x,y)
    elif isinstance(a,float):np.testing.assert_allclose(a,b,rtol=0,atol=1e-12)
    else:assert a==b

def verify():
    v=read(HERE/'verification.json');p=read(FRESH/'protocol.json')
    saved=read(HERE/'pilot-audit-112.json')
    assert sha(FRESH/'protocol.json')==v['protocol_sha256']
    assert sha(HERE/'pilot-audit-112.json')==v['pilot_audit_sha256']
    assert rh.sources()==v['source_hashes']
    inventory=read(HERE/'source-inventory.json')
    for rel,digest in inventory['fresh_implementation'].items():assert sha(FRESH/rel)==digest,rel
    for rel,digest in inventory['predecessor'].items():assert sha(FRESH.parent/rel)==digest,rel
    for e in v['pilot_archives']:
        path=FRESH/'recorded-data'/(e['name']+'.npz')
        m=rh.validate(path,e['config'])
        assert sha(path)==e['sha256'] and sha(path.with_suffix('.json'))==e['metadata_sha256']
        assert m['study_sources']==v['source_hashes']
    oldpath=FRESH.parent/'results/primary_models.json'
    assert sha(oldpath)==saved['old_models_sha256'];old=read(oldpath)
    seeds=p['pilot']['seeds'];half=p['pilot']['half_step_seeds']
    data={n:na.datasets(seeds,n,p) for n in [80,112,160]}
    rows=na.pair_rows(data[112],data[160],p,p['pilot']['event_cutoff'],old)
    temporal=na.pair_rows(na.datasets(half,112,p),na.datasets(half,112,p,p['pilot']['half_step_dt']),p,p['pilot']['event_cutoff'],old)
    before=na.pair_rows(data[80],data[112],p,p['pilot']['event_cutoff'],old)
    screens=na.screens_for(rows,p,'112\u2192160')+na.volumes(seeds,112,p)+na.volumes(seeds,160,p)
    screens+=na.screens_for(temporal,p,'112 half-step',True,rows)+na.volumes(half,112,p,p['pilot']['half_step_dt'])
    screens+=na.screens_for(before,p,'80\u2192112')
    for a,b in zip(rows,before):
        x=a['target_RMS']/max(b['target_RMS'],1e-30)
        screens.append({'comparison':'adjacent contraction','seed':a['seed'],'metric':'target_RMS_ratio','value':x,'bound':1.,'direction':'max','pass':x<1.})
    # Absolute roundoff allowance for reproduction on another BLAS/platform only.
    # The accuracy bounds and pass decisions below remain exactly preregistered.
    same(rows,saved['rows']);same(temporal,saved['temporal']);same(screens,saved['screens'])
    for s in screens:
        metric=s['metric'];direction='min' if metric=='event_Jaccard' else 'max'
        if metric=='target_RMS_ratio':bound=1.;passed=s['value']<bound
        else:
            key=('half_step_target_relative_RMS_max' if 'half-step' in s['comparison'] else 'target_relative_RMS_max') if metric=='target_relative_RMS' else ('event_Jaccard_min' if direction=='min' else ('frozen_old_prediction_relative_RMS_max' if metric.startswith(('A_full_','B_history_')) else metric+'_max'))
            bound=p['numeric_screens'][key]
            passed=s['value']>=bound if direction=='min' else s['value']<=bound
        assert s['bound']==bound and s['direction']==direction and s['pass']==passed and passed
    assert len(screens)==116
    # Independent target and feature RMS aggregation from the frozen descriptors.
    for a,b,rr in [(data[112],data[160],rows),(data[80],data[112],before)]:
        for row in rr:
            use=a['group']==row['seed'];rms=lambda x:float(np.sqrt(np.mean(x*x)))
            assert rms(a['target'][use]-b['target'][use])==row['target_RMS']
            for label,scale in [('present',np.array(old['models']['A_full']['continuous']['scale'])),('history',np.array(old['models']['B_history']['continuous']['scale'])[31:])]:
                assert rms((a[label][use]-b[label][use])/scale)==row[label+'_original_train_standardized_RMS']
    print('22 pilot archives verified; 116/116 screens reproduced within1e-12; frozen thresholds pass. No simulations or fits run.')

if __name__=='__main__':verify()
