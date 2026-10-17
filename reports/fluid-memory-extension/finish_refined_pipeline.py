"""Wait for all original and additional tests before interpreting memory."""
import time,subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for name in ['gpu-pipeline-complete.json','cpu-pipeline-complete.json','refinement-complete.json','refined-observation-validation.json']:
    deadline=time.monotonic()+86400
    def ready():
        p=ROOT/'results'/name
        if not p.exists():return False
        if name=='cpu-pipeline-complete.json':
            plan=json.loads((ROOT/'fresh-grid-refinement.json').read_text())
            return len(json.loads(p.read_text())['completed'])==12+plan.get('additional_independent_comparisons',0)
        return True
    while not ready():
        for failure in ['gpu-pipeline-failure.json','cpu-pipeline-failure.json','refinement-failure.json','refined-cpu-failure.json','refined-observation-failure.json']:
            if (ROOT/'results'/failure).exists():raise RuntimeError('Pipeline failure: '+failure)
        if time.monotonic()>deadline:raise RuntimeError('Completion wait timed out')
        time.sleep(10)
subprocess.run([sys.executable,'final_analysis.py'],cwd=ROOT,check=True)
print('All original and refined results verified. External publication remains.',flush=True)
