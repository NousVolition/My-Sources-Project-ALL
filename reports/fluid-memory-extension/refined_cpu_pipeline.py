"""Run two refined checks while the original CPU controller remains active.

Unique files prevent duplicate writers. This controller joins both completion
records only after the original twelve and the additional two finish. For fresh
reruns, parallel_cpu_pipeline.py performs the same fourteen checks directly.
"""
import json,time,traceback
from threadpoolctl import threadpool_limits
from final_experiment import ROOT,dump,digest
from refined_independent_check import check

def main():
    completed=[]
    for tight in [False,True]:
        with threadpool_limits(1):r=check(82000,'forced_nu002',tight)
        completed.append(dict(seed=82000,regime='forced_nu002',tight=tight,kind='refined_full'))
        dump(ROOT/'results'/'refined-cpu-progress.json',dict(completed=completed,planned=2))
    dump(ROOT/'results'/'refined-cpu-complete.json',dict(status='complete',completed=completed,controller_sha256=digest(ROOT/'refined_cpu_pipeline.py')))
    deadline=time.monotonic()+86400
    original=ROOT/'results'/'cpu-pipeline-complete.json'
    while not original.exists():
        if (ROOT/'results'/'cpu-pipeline-failure.json').exists():raise RuntimeError('Original independent CPU pipeline failed')
        if time.monotonic()>deadline:raise RuntimeError('Original independent CPU wait timed out')
        time.sleep(10)
    info=json.loads(original.read_text());assert len(info['completed'])==12
    dump(ROOT/'results'/'cpu-original-complete.json',info)
    info.update(completed=info['completed']+completed,additional_refined_comparisons=2,
        refinement_controller_sha256=digest(ROOT/'refined_cpu_pipeline.py'))
    dump(original,info)
    print('All fourteen independent comparisons complete.',flush=True)

if __name__=='__main__':
    try:main()
    except Exception:
        dump(ROOT/'results'/'refined-cpu-failure.json',dict(status='failed',traceback=traceback.format_exc()));raise
