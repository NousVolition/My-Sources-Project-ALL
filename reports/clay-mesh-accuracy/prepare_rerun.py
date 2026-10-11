"""Copy exact sources into a clean directory without old computed recordings."""
from pathlib import Path
import argparse,shutil
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('destination',type=Path);a=ap.parse_args();dest=a.destination.resolve()
    if dest.exists():raise ValueError('Destination must be new; never overwrite existing research')
    dest.mkdir(parents=True)
    for p in ROOT.rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts:continue
        rel=p.relative_to(ROOT)
        if p.suffix in ['.py','.cjs','.md','.txt'] or p.name in ['molding-template.html','executed-protocol.json'] or 'historical' in rel.parts:
            out=dest/rel;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
    for part in ['','model','refined','production','resolved','solver']:(dest/part/'data').mkdir(parents=True,exist_ok=True)
    print(dest)
if __name__=='__main__':main()
