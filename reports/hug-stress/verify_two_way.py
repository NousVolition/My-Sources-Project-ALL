"""Audit completed two-way results against saved trajectories and source hashes."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.integrate import trapezoid
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'two-way-data'
def read(n):return json.loads((DATA/n).read_text())
def main():
    protocol=read('protocol.json');controls=read('controls.json');rows=read('paired.json');stress=read('stress.json');opening=read('opening-tests.json')
    for name,digest in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((ROOT/'test_two_way_opening.py').read_bytes()).hexdigest()==opening['source_sha256']
    assert len(controls)==7 and all(c['passed'] for c in controls)
    assert len(rows)==12 and len(stress)==12 and len(opening['cases'])==12 and opening['all_passed']
    for row in rows:
        with np.load(DATA/(row['name']+'.npz')) as a:
            assert a['time'][0]==0 and a['time'][-1]==80
            mask=a['time']>=20;t=a['time'][mask]
            for key,column in [('state','metrics'),('refined','refined')]:
                y=a[key][mask];rms=np.sqrt(trapezoid(np.sum(y[:,:2]**2,axis=1),t)/(t[-1]-t[0]))
                assert abs(rms-row[column]['activity_rms'])<1e-12
                assert row[column]['memory_bounds'] and row[column]['passive_feedback']
                assert row[column]['scaled_activity_budget_error']<1e-7
    for p in DATA.glob('*.npz'):
        with np.load(p) as a:assert all(np.isfinite(a[k]).all() for k in a.files)
    summary=read('summary.json')
    assert summary['all_bounds_and_passivity'] and summary['constant_pressure_hug_independence_error']<1e-7
    assert len(summary['effects'])==3
    for s in summary['effects']:
        for metric,effect in s['effects'].items():
            changes=[]
            for seed in range(3):
                one=next(r for r in rows if r['seed']==seed and r['strength']==0)
                two=next(r for r in rows if r['seed']==seed and r['strength']==s['strength'])
                changes.append(100*(two['metrics'][metric]/one['metrics'][metric]-1))
            assert np.allclose(changes,effect['percent_changes'],rtol=0,atol=1e-12)
    result=dict(audit='passed',paired_configurations=12,stress_configurations=12,short_controls=7,opening_configurations=12,
        integrations=84,physical_measurements=0,baseline_model_sha256=hashlib.sha256((ROOT/'baseline/hug_model.py').read_bytes()).hexdigest())
    (DATA/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
