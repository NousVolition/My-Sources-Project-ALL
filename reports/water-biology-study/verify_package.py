"""Check integrity and reproducibility of the delivered research record."""
from pathlib import Path
import argparse,csv,hashlib,json,math,re
import numpy as np

ROOT=Path(__file__).resolve().parent
# These remain hash-verified package inputs/derived outputs, but run_study.py
# does not generate them: one literature transcription, one separate transient
# audit, and three tables made by analyze.py.
NON_RUN_DATA={'published_budke_koop_table_A2.csv','scalar_transient_convergence.csv',
              'freezing_conditional_probabilities.csv','freezing_simulated_survival.csv',
              'solute_particle_overlap_seed100.csv'}

RTOL=1e-12
ATOL=1e-13


def compare_number(a,b,where,stats):
    if not math.isfinite(a) or not math.isfinite(b):
        raise ValueError(f'{where}: nonfinite reproduced number')
    difference=abs(b-a); tolerance=ATOL+RTOL*abs(a)
    stats['numeric_values_checked']+=1
    stats['numeric_values_different']+=int(a!=b)
    if difference>stats['max_absolute_difference']:
        stats['max_absolute_difference']=difference;stats['largest_absolute_difference_at']=where
    stats['max_fraction_of_tolerance']=max(stats['max_fraction_of_tolerance'],difference/tolerance)
    if difference>tolerance:
        raise ValueError(f'{where}: file contents differ: expected {a!r}, got {b!r}, tolerance {tolerance!r}')


def compare_json(a,b,where,stats):
    if isinstance(a,dict) and isinstance(b,dict):
        if set(a)!=set(b):raise ValueError(f'{where}: JSON keys differ')
        for key in a:compare_json(a[key],b[key],f'{where}/{key}',stats)
    elif isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b):raise ValueError(f'{where}: JSON lengths differ')
        for i,(x,y) in enumerate(zip(a,b)):compare_json(x,y,f'{where}/{i}',stats)
    elif type(a) is float and type(b) in (float,int):compare_number(a,b,where,stats)
    elif type(a)!=type(b) or a!=b:
        raise ValueError(f'{where}: file contents differ: {a!r} versus {b!r}')


def compare_text(original,other,stats):
    if original.suffix=='.json':
        compare_json(json.loads(original.read_text(encoding='utf-8')),
                     json.loads(other.read_text(encoding='utf-8')),original.name,stats)
    elif original.suffix=='.csv':
        with original.open(newline='',encoding='utf-8') as a,other.open(newline='',encoding='utf-8') as b:
            left=list(csv.reader(a));right=list(csv.reader(b))
        if not left or not right or len(left)!=len(right) or left[0]!=right[0]:raise ValueError(f'{original.name}: CSV rows or header differ')
        for row,(xs,ys) in enumerate(zip(left[1:],right[1:]),2):
            if len(xs)!=len(ys):raise ValueError(f'{original.name}:{row}: CSV columns differ')
            for col,(x,y) in enumerate(zip(xs,ys)):
                where=f'{original.name}:{row}:{left[0][col]}'
                if re.fullmatch(r'[+-]?\d+',x):
                    if not re.fullmatch(r'[+-]?\d+',y) or int(x)!=int(y):raise ValueError(f'{where}: integer differs')
                    continue
                try: number=float(x)
                except ValueError:
                    if x!=y:raise ValueError(f'{where}: label differs')
                else:
                    try:actual=float(y)
                    except ValueError:raise ValueError(f'{where}: expected a number') from None
                    compare_number(number,actual,where,stats)
    elif original.read_bytes().replace(b'\r\n',b'\n')!=other.read_bytes().replace(b'\r\n',b'\n'):
        raise ValueError(f'{original.name}: file contents differ')


def included():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not any(part in ['__pycache__','.pytest_cache','reproduction'] for part in p.relative_to(ROOT).parts)
                  and p.relative_to(ROOT).as_posix() not in ['manifest.json','package-verification.json'])


def compare_reproduction(actual_root, manifest):
    stats={'numeric_values_checked':0,'numeric_values_different':0,'max_absolute_difference':0.,'max_fraction_of_tolerance':0.}
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
                    np.testing.assert_allclose(b[key],a[key],rtol=RTOL,atol=ATOL,equal_nan=False,err_msg=f'{name}:{key}')
        else:compare_text(original,other,stats)
    print(json.dumps({'text_reproduction_comparison':stats},indent=2))
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
