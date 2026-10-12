from pathlib import Path
import argparse,shutil
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('destination',type=Path);a=ap.parse_args();dest=a.destination.resolve()
    if dest.exists():raise ValueError('Destination must be new')
    dest.mkdir(parents=True)
    for p in ROOT.rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts:continue
        rel=p.relative_to(ROOT)
        if p.suffix in ['.py','.cjs','.md','.txt','.css'] or p.name in ['executed-protocol.json','molding-template.html'] or rel.parts[0]=='parents':
            out=dest/rel;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
    (dest/'solver/data').mkdir(parents=True,exist_ok=True);(dest/'data').mkdir(exist_ok=True)
    print(dest)
