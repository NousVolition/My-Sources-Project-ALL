"""Audit saved evidence and provenance; numerical failures must remain visible."""
from pathlib import Path
import json,hashlib,sys
import numpy as np

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def read(path):return json.loads(path.read_text())
def main():
    rows=read(DATA/'screen.json');summary=read(DATA/'summary.json');protocol=read(DATA/'protocol.json')
    assert len(rows)==len(protocol['cases'])==224
    assert len({r['name'] for r in rows})==224
    assert sum(r['group']=='matrix' for r in rows)==144
    assert sum(r['false_convergence'] for r in rows)==summary['false_convergence']==52
    assert sum(r['extent']['outside_cube'] for r in rows)==summary['outside_cube']==29
    for path,digest in protocol['source_hashes'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((ROOT/'run_stress.py').read_bytes()+(ROOT/'stress_solver.py').read_bytes()).hexdigest()==protocol['runner_fingerprint']
    for r,c in zip(rows,protocol['cases']):
        assert all(r[k]==v for k,v in c.items())
        assert read(DATA/(r['name']+'.json'))==r
        with np.load(DATA/(r['name']+'.npz')) as a:
            for tag in ('coarse','fine','reference'):
                t=a[tag+'_time'];y=a[tag+'_state']
                assert y.shape==(len(t),5) and np.isfinite(y).all() and np.all(np.diff(t)>0)
                assert t[0]==0 and np.allclose(y[0],[r['q0'],r['v0'],0,0,0],rtol=0,atol=1e-12)
                if tag=='reference':assert t[-1]==r['end']
                else:
                    f=r['fixed'][0 if tag=='coarse' else 1]
                    assert (t[-1]<r['end'])==bool(f['failure'])
                    assert np.max(abs(y))<=1e12
    independent=read(DATA/'independent-solvers.json');tight=read(DATA/'long-refinement.json');controls=read(DATA/'controls.json')
    assert len(independent)==43 and sum(x['passed'] for x in independent)==41
    assert {x['name'] for x in independent if not x['passed']}=={x['name'] for x in tight}
    assert len(tight)==2 and all(x['passed'] and x['scaled_discrepancy']<2e-6 for x in tight)
    assert all(x['rk4_errors']['fine']<=1e-3<x['rk4_errors']['coarse'] for x in tight)
    assert len(controls)==11 and all(x['passed'] for x in controls)
    assert max(tight[0]['exact_errors'].values())<2e-6 and tight[0]['exact_initial_error']<1e-12
    ld=ROOT/'lorenz-data';ls=read(ld/'summary.json');lp=read(ld/'protocol.json')
    assert ls['all_checks_passed'] and ls['configurations']==84
    assert len(read(ld/'screen.json'))==84
    for name,digest in lp['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    assert all(x['passed'] for x in read(ld/'sweep-independent.json')) and len(read(ld/'sweep-independent.json'))==7
    assert read(ld/'event-check.json')['passed']
    matched=read(ld/'matched-control.json');assert matched['passed']
    assert hashlib.sha256((ROOT/'align_lorenz_control.py').read_bytes()).hexdigest()==matched['source_sha256']
    with np.load(ld/'chaotic.npz') as a,np.load(ld/'no-feedback-matched.npz') as b:
        assert np.array_equal(a['pressure'],b['pressure'])
        assert np.array_equal(a['state'][:,:3],b['state'][:,:3]) and np.array_equal(a['state'][:,5],b['state'][:,5])
    assert len(read(ld/'lyapunov.json'))==4 and ls['lyapunov_min']>0
    for path in ld.glob('*.npz'):
        with np.load(path) as arrays:
            assert all(np.isfinite(arrays[key]).all() for key in arrays.files)
    manifest=ROOT/'checksums.json';hashed=0
    if manifest.exists():
        for name,digest in read(manifest).items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest;hashed+=1
    result=dict(evidence_audit='passed',cases=224,original_reference_checks_passed=41,original_reference_checks_total=43,
                refined_reference_checks_passed=2,controls_passed=11,original_failures_preserved=True,lorenz_cases=84,
                lorenz_independent_presets=5,lorenz_independent_sweep_cases=7,lorenz_opening_event_checked=True,manifest_files_verified=hashed)
    if '--record' in sys.argv:
        result.pop('manifest_files_verified');(DATA/'package-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
