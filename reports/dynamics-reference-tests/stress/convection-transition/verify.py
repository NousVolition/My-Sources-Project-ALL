"""Audit delivered current arrays and provenance; rejected history is separate."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'


def main():
    summary=json.loads((DATA/'summary.json').read_text())
    for name,sha in summary['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    with (DATA/'cases.csv').open(newline='') as f:
        cases=list(csv.DictReader(f))
    assert len(cases)==summary['pde_runs']==len({r['label'] for r in cases})
    for row in cases:
        with np.load(DATA/(row['label']+'.npz')) as a:
            assert all(np.isfinite(a[k]).all() for k in a.files),row['label']
            r=a['records'];fields=a['fields']
            assert abs(r[-1,0]-float(row['duration']))<1e-12
            assert fields.shape==(3,3,int(row['n']),int(row['n']))
            for f in fields:
                reflected=f[:,:,(-np.arange(f.shape[-1]))%f.shape[-1]]
                scale=max(float(np.max(abs(f))),1e-300)
                assert np.max(abs(f[0]-reflected[0]))/scale<1e-12
                assert np.max(abs(f[1:]+reflected[1:]))/scale<1e-12
            for energy,work,diss,residual in [(3,5,6,9),(4,7,8,10)]:
                computed=r[:,energy]-r[0,energy]-r[:,work]+r[:,diss]
                np.testing.assert_allclose(computed,r[:,residual],atol=1e-25,rtol=1e-12)
    modal=json.loads((DATA/'inertia.json').read_text())
    for row in modal:
        with np.load(DATA/(row['label']+'.npz')) as a:
            assert all(np.isfinite(a[k]).all() for k in a.files)
    expected={r['label']+'.npz' for r in cases+modal}
    assert {p.name for p in DATA.glob('*.npz')}==expected
    archive=list((DATA/'rejected-boundary-run').glob('*.npz'))
    output=dict(passed=True,pde_recordings=len(cases),modal_recordings=len(modal),
                current_arrays_finite=True,source_fingerprints_match=True,
                saved_field_parity_verified=True,budgets_recomputed=True,
                rejected_recordings_excluded=len(archive))
    (DATA/'delivery-audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output))

if __name__=='__main__':main()
