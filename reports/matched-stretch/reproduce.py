"""Generate starting arrays and run the recorded matrix in a fresh directory."""
from pathlib import Path
import subprocess
import sys
import numpy as np
from numerics import Flow
p=Path(__file__).resolve().parent
for n in (64,128,256):
    f=p/f'initial-n{n}.npy'
    if not f.exists(): np.save(f,Flow(n).initial())
subprocess.run([sys.executable,str(p/'run_suite.py'),'--workers','2'],check=True)
