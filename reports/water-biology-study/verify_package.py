"""Check integrity and reproducibility of the delivered research record."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent
# These remain hash-verified package inputs/derived outputs, but run_study.py
# does not generate them: one literature transcription, one separate transient
# audit, and three tables made by analyze.py.
NON_RUN_DATA={'published_budke_koop_table_A2.csv','scalar_transient_convergence.csv',
              'freezing_conditional_probabilities.csv','freezing_simulated_survival.csv',
              'solute_particle_overlap_seed100.csv'}


def included():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not any(part in ['__pycache__','.pytest_cache','reproduction'] for part in p.relative_to(ROOT).parts)
                  and p.relative_to(ROOT).as_posix() not in ['manifest.json','package-verification.json'])


def compare_reproduction(actual_root, manifest):
    expected={name.removeprefix('data/') for name in manifest if name.startswith('data/')} - NON_RUN_DATA
    if len(expected)!=48:
        raise ValueError(f'Expected 48 documented raw-data files, manifest lists {len(expected)}')
    data=actual_root/'data'
    if not data.is_dir():
        raise ValueError('Reproduction data directory is missing')
    actual={p.relative_to(data).as_posix() for p in data.rglob('*') if p.is_file()}
    if actual!=expected:
        raise ValueError(f'Reproduction file set differs: missing={sorted(expected-actual)}, unexpected={sorted(actual-expected)}')
    for name in sorted(expected):
        original=ROOT/'data'/name;other=data/name
        if other.suffix=='.npz':
            with np.load(original,allow_pickle=False) as a,np.load(other,allow_pickle=False) as b:
                if set(a.files)!=set(b.files):
                    raise ValueError(f'{name}: array names differ')
                for key in a.files:
                    if a[key].shape!=b[key].shape:
                        raise ValueError(f'{name}:{key}: array shapes differ')
                    np.testing.assert_allclose(a[key],b[key],rtol=1e-12,atol=1e-13,equal_nan=False,err_msg=f'{name}:{key}')
        elif original.read_bytes()!=other.read_bytes():
            raise ValueError(f'{name}: file contents differ')
    return len(expected)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--make-manifest',action='store_true'); ap.add_argument('--compare',type=Path)
    args=ap.parse_args()
    if args.make_manifest:
        payload={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in included()}
        (ROOT/'manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    manifest=json.loads((ROOT/'manifest.json').read_text())
    actual_included={p.relative_to(ROOT).as_posix() for p in included()}
    if set(manifest)!=actual_included:
        raise ValueError(f'Manifest coverage differs: unlisted={sorted(actual_included-set(manifest))}, missing={sorted(set(manifest)-actual_included)}')
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
        compared=compare_reproduction(args.compare,manifest)
    result={'sha256_files_verified':len(manifest),'budgets':'passed','finite_arrays':'passed','independent_reproduction_files_compared':compared}
    (ROOT/'package-verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
