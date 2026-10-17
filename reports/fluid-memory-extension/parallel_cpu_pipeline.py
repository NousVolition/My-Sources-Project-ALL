"""Original twelve plus two refined independent comparisons, two jobs at once.

Use only after the sequential controller has stopped. Unique normal/tight
output files keep jobs independent; existing results undergo cache verification.
No equations, settings, sampling or statistical criteria change.
"""
import concurrent.futures as futures,json,time
from final_experiment import ROOT,P,dump,digest
from analysis import effective_regimes,filename

def matrix():
    tasks=[]
    for regime,seeds in [('smooth',P['smooth_seeds'][:2]),('forced_nu004',P['forced_seeds'][:1]),('forced_nu002',P['forced_seeds'][:1])]:
        for seed in seeds:
            for tight in [False,True]:tasks.append(dict(seed=seed,regime=regime,tight=tight,kind='full'))
    for regime,seed in [('smooth',P['smooth_seeds'][0]),('forced_nu002',P['forced_seeds'][0])]:
        for tight in [False,True]:tasks.append(dict(seed=seed,regime=regime,tight=tight,kind='common'))
    if (ROOT/'fresh-grid-refinement.json').exists():
        plan=json.loads((ROOT/'fresh-grid-refinement.json').read_text())
        if plan.get('additional_independent_comparisons')==2:
            for tight in [False,True]:tasks.append(dict(seed=82000,regime='forced_nu002',tight=tight,kind='refined_full'))
    return tasks

def input_ready(task):
    n=32 if task['regime']=='smooth' else max(64,min(effective_regimes()[task['regime']]['grids']))
    if task['kind']=='refined_full':n=96
    if task['kind'] in ['full','refined_full']:p=filename(task['seed'],task['regime'],'coincident',n)
    else:
        dt=effective_regimes()[task['regime']]['dt']
        p=ROOT/'common-probe-data'/f"{task['regime']}-coincident-s{task['seed']}-n{n}-dt{dt:g}.npz"
    return p.with_suffix('.json').exists()

def run_task(task):
    from threadpoolctl import threadpool_limits
    from independent_check_full import check
    from independent_common import check as common
    from refined_independent_check import check as refined
    with threadpool_limits(1):
        result={'full':check,'common':common,'refined_full':refined}[task['kind']](task['seed'],task['regime'],task['tight'])
    return dict(task=task,passed=result['passed'])

def main():
    pending=matrix();planned=len(pending);completed=[];running={};deadline=time.monotonic()+86400
    with futures.ProcessPoolExecutor(max_workers=2) as pool:
        while pending or running:
            if (ROOT/'results'/'gpu-pipeline-failure.json').exists():raise RuntimeError('GPU pipeline failed')
            if time.monotonic()>deadline:raise RuntimeError('Independent comparison wait timed out')
            for task in pending.copy():
                if len(running)>=2:break
                if input_ready(task):running[pool.submit(run_task,task)]=task;pending.remove(task)
            if not running:time.sleep(10);continue
            done,_=futures.wait(running,timeout=10,return_when=futures.FIRST_COMPLETED)
            for future in done:
                result=future.result();completed.append(running.pop(future))
                dump(ROOT/'results'/'cpu-pipeline-progress.json',dict(completed=completed,planned=planned,schedule='At most two independent CPU jobs'))
    dump(ROOT/'results'/'cpu-pipeline-complete.json',dict(status='complete',completed=completed,controller_sha256=digest(ROOT/'parallel_cpu_pipeline.py'),
        additional_refined_comparisons=planned-12,original_controller_sha256=digest(ROOT/'cpu_pipeline.py')))

if __name__=='__main__':
    try:main()
    except Exception:
        import traceback
        dump(ROOT/'results'/'cpu-pipeline-failure.json',dict(status='failed',traceback=traceback.format_exc()))
        raise
