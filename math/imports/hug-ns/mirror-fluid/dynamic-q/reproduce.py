"""Rerun the declared 26 jobs in a new directory, preserving supplied results."""
from pathlib import Path
import shutil, subprocess, sys
ROOT=Path(__file__).resolve().parent
dest=ROOT/'rerun'
if dest.exists():raise SystemExit('rerun already exists; inspect its saved state before resuming.')
dest.mkdir()
for name in ['study.py','analyze.py']:
    shutil.copy2(ROOT/name,dest/name)
shutil.copytree(ROOT/'source',dest/'source',ignore=shutil.ignore_patterns('__pycache__'))
for arguments in [['study.py','--audit'],['study.py','--workers','2'],['analyze.py']]:
    subprocess.run([sys.executable,*arguments],cwd=dest,check=True)
