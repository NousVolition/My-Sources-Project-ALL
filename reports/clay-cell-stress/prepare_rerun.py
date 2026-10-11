"""Create a fresh experiment directory without reusing computed recordings."""
from pathlib import Path
import argparse,shutil
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('destination',type=Path);args=ap.parse_args()
    target=args.destination.resolve()
    if target.exists():raise SystemExit('Choose a new destination; existing files are never replaced.')
    target.mkdir(parents=True);(target/'data').mkdir()
    for p in ROOT.iterdir():
        if p.is_file() and (p.suffix in ['.py','.cjs','.txt','.md'] or p.name=='molding-template.html'):
            shutil.copy2(p,target/p.name)
    for directory in ['solver-versions','initial-attempts']:
        shutil.copytree(ROOT/directory,target/directory,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for direction in ['one','two']:
        for suffix in ['.npz','.json']:
            name='prior-extreme-'+direction+suffix;shutil.copy2(ROOT/'data'/name,target/'data'/name)
    print(target)
if __name__=='__main__':main()
