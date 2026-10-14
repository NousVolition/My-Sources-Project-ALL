"""Make a result-free copy with frozen sources, inputs and exact run protocol."""
from pathlib import Path
import sys,shutil
ROOT=Path(__file__).resolve().parent
target=Path(sys.argv[1]).resolve();assert not target.exists();target.mkdir()
for d in ['solver','legacy_source','checkpoint_source','parents','data','figures']:(target/d).mkdir()
for d in ['solver','legacy_source','checkpoint_source']:
 for p in (ROOT/d).glob('*.py'):shutil.copy2(p,target/d/p.name)
for folder in (ROOT/'solver-versions').glob('*'):
 if not folder.is_dir():continue
 dest=target/'solver-versions'/folder.name;dest.mkdir(parents=True);(dest/'data').mkdir()
 for p in folder.glob('*.py'):shutil.copy2(p,dest/p.name)
for p in (ROOT/'parents').glob('*'):
 if p.is_file() and p.suffix in ['.json','.npz']:shutil.copy2(p,target/'parents'/p.name)
for p in ROOT.glob('*'):
 if p.is_file() and p.suffix in ['.py','.md','.txt','.html','.css'] and p.name not in ['report.html','molding-fragment.html']:shutil.copy2(p,target/p.name)
shutil.copy2(ROOT/'data/executed-protocol.json',target/'data/executed-protocol.json')
(target/'solver/data').mkdir();print(target)
