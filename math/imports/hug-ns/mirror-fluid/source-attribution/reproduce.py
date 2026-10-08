"""Reproduce source controls in rerun/, keeping published measurements intact."""
from pathlib import Path
import shutil
import subprocess
import sys
here=Path(__file__).resolve().parent
target=here/'rerun'
target.mkdir(exist_ok=True)
for name in ('measure_sources.py','analyze_sources.py','reference-fluid-sources.json'):
    dest=target/name
    if dest.exists() and dest.read_bytes()!=(here/name).read_bytes():
        raise RuntimeError(f'{dest} differs; use a fresh directory to keep earlier runs intact.')
    shutil.copyfile(here/name,dest)
(target/'source').mkdir(exist_ok=True)
for p in (here/'source').glob('*.py'):
    dest=target/'source'/p.name
    if dest.exists() and dest.read_bytes()!=p.read_bytes():
        raise RuntimeError(f'{dest} differs; use a fresh directory.')
    shutil.copyfile(p,dest)
subprocess.run([sys.executable,str(target/'measure_sources.py'),'--workers','2'],check=True)
subprocess.run([sys.executable,str(target/'analyze_sources.py')],check=True)
