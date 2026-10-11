"""Sequential GPU follow-ups and final verification/report workflow."""
import os,subprocess,sys
from pathlib import Path

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1'
    for script in ['matched_sham.py','delayed_probe.py','followup_analysis.py','template_analysis.py','impulse_energy.py','verify.py','report.py']:
        print('Running '+script,flush=True)
        subprocess.run([sys.executable,'-B',script],cwd=root,env=env,check=True)
