"""Check integrity and reproducibility of the delivered research record."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent


def included():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not any(part in ['__pycache__','.pytest_cache','reproduction'] for part in p.relative_to(ROOT).parts)
                  and p.name not in ['manifest.json','package-verification.json'])


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--make-manifest',action='store_true'); ap.add_argument('--compare',type=Path)
    args=ap.parse_args()
    if args.make_manifest:
        payload={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in included()}
        (ROOT/'manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    for p in (ROOT/'data').glob('*.npz'):
        with np.load(p) as z:
            assert all(np.isfinite(z[k]).all() for k in z.files),p.name
    with open(ROOT/'data/column_budgets.csv',newline='') as f:
        for r in csv.DictReader(f):
            assert sum(int(r[k]) for k in ['suspended','deposited','attached'])==8192
    with open(ROOT/'data/particle_timeseries.csv',newline='') as f:
        for r in csv.DictReader(f): assert float(r['count'])==8192
    compared=0
    if args.compare:
        for other in sorted((args.compare/'data').iterdir()):
            original=ROOT/'data'/other.name
            assert original.exists(),other.name
            if other.suffix=='.npz':
                with np.load(original) as a,np.load(other) as b:
                    assert a.files==b.files
                    for k in a.files: np.testing.assert_allclose(a[k],b[k],rtol=1e-12,atol=1e-13,err_msg=f'{other.name}:{k}')
            else:
                assert original.read_bytes()==other.read_bytes(),other.name
            compared+=1
    result={'sha256_files_verified':len(manifest),'budgets':'passed','finite_arrays':'passed','independent_reproduction_files_compared':compared}
    (ROOT/'package-verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
