"""Complete the pre-effect fresh-flow refinement after the original matrix."""
import json,time,traceback
from final_experiment import ROOT,run_one,dump,digest
from common_probe import run_one as common
from unforced_readout import run_one as off
from common_unforced import run_one as common_off
from safe_common_probe import run_one as bounded_common
from safe_unforced_readout import run_one as bounded_off
from safe_common_unforced import run_one as bounded_common_off
from bounded_final_experiment_v2 import run_one as bounded_main

def main():
    deadline=time.monotonic()+86400
    while not (ROOT/'results'/'gpu-pipeline-complete.json').exists():
        if (ROOT/'results'/'gpu-pipeline-failure.json').exists():raise RuntimeError('Original GPU pipeline failed')
        if time.monotonic()>deadline:raise RuntimeError('Original matrix wait timed out')
        time.sleep(10)
    from validate_bounded_streaming import validate
    validate()
    file=ROOT/'fresh-grid-refinement.json';plan=json.loads(file.read_text());completed=[]
    methods={'main':run_one,'common':common,'unforced':off,'common_unforced':common_off}
    bounded={'main':bounded_main,'common':bounded_common,'unforced':bounded_off,'common_unforced':bounded_common_off}
    for task in plan['tasks']:
        args=dict(task);kind=args.pop('kind');result=(bounded if task['n']>=128 else methods)[kind](**args)
        completed.append(dict(task=task,config=result['config'],raw_sha256=result['raw_sha256']))
        dump(ROOT/'results'/'refinement-progress.json',dict(planned=len(plan['tasks']),completed=completed,plan_sha256=digest(file)))
    dump(ROOT/'results'/'refinement-complete.json',dict(status='complete',completed=completed,planned=len(plan['tasks']),
        plan_sha256=digest(file),controller_sha256=digest(ROOT/'refinement_gpu_pipeline.py')))

if __name__=='__main__':
    try:main()
    except Exception:
        dump(ROOT/'results'/'refinement-failure.json',dict(status='failed',traceback=traceback.format_exc()));raise
