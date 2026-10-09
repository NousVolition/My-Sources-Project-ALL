"""Repeat four 65-grid runs in a new folder; retain published reference data."""
from pathlib import Path
import shutil,subprocess,sys
root=Path(__file__).resolve().parent;dest=root/'rerun65'
if dest.exists():raise SystemExit('rerun65 already exists; inspect it before resuming.')
dest.mkdir()
for name in ['study.py','refine65.py','analyze.py','analyze65.py','training-freeze.json','fitted-training-models.json','fitted-record-models.json','protocol.json']:
    shutil.copy2(root/name,dest/name)
shutil.copytree(root/'source',dest/'source',ignore=shutil.ignore_patterns('__pycache__'))
for p in (root/'runs').glob('*/result.json'):
    if p.parent.name.endswith('-n65'):continue
    out=dest/'runs'/p.parent.name/p.name;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
for script in ['refine65.py','analyze65.py']:subprocess.run([sys.executable,script],cwd=dest,check=True)
