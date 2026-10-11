"""Reproduce numerical studies, figures, report and audit from shipped checkpoints."""
import subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
    for args in [('long_cycles.py',),('ramps.py',),('ramps.py','--refine'),('noise.py',),('capacitance.py',),('prediction.py',),('run_fluid.py',),('prediction.py','--fluid'),('fluid_confirmation.py',),('report.py',),('verify.py',)]:
        subprocess.run([sys.executable,str(ROOT/args[0]),*args[1:]],check=True,cwd=ROOT)
