"""Explicitly dated follow-up; do not run concurrently with run_gpu.py."""
import json
from run import ROOT,dump,digest
import run_gpu
from refinement_diagnostics import RefinementSolver,validate_diagnostics
if __name__=='__main__':
    checks=validate_diagnostics();dump(ROOT/'refinement-diagnostics-validation.json',dict(checks=checks,passed=all(c['passed'] for c in checks)))
    assert all(c['passed'] for c in checks),'GPU diagnostics must match CPU definitions'
    run_gpu.GPUSolver=RefinementSolver
    for task in json.loads((ROOT/'refinement-plan.json').read_text())['cases']:
        print(run_gpu.run_gpu(task),flush=True)
        path=ROOT/'data-gpu'/run_gpu.name(task)/'result.json'
        result=json.loads(path.read_text());result['diagnostics_backend']='GPU, checked against CPU'
        for name in ['refine.py','refinement_diagnostics.py','refinement-plan.json']:
            result['source_sha256'][name]=digest(ROOT/name)
        dump(path,result)
