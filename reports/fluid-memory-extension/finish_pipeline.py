"""Finish local deliverables after both running numerical pipelines complete."""
import time,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for name in ['gpu-pipeline-complete.json','cpu-pipeline-complete.json']:
    deadline=time.monotonic()+86400
    while not (ROOT/'results'/name).exists():
        if (ROOT/'results'/'gpu-pipeline-failure.json').exists():raise RuntimeError('Simulation failure recorded; no completed-study report generated')
        if time.monotonic()>deadline:raise RuntimeError('Completion wait timed out')
        time.sleep(10)
for script in ['followup_analysis.py','report_extension.py','verify_extension.py']:
    subprocess.run([sys.executable,script],cwd=ROOT,check=True)
print('Local results, figures and verification complete. External publication remains a separate step.',flush=True)
