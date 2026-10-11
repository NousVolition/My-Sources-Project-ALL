import time,json
from threadpoolctl import threadpool_limits
from experiment import ROOT,P,dump,digest
from analysis import filename,effective_regimes
from independent_check_full import check
from independent_common import check as common

def wait_file(path):
    deadline=time.time()+21600
    while not path.with_suffix('.json').exists():
        if (ROOT/'results'/'gpu-pipeline-failure.json').exists():raise RuntimeError('GPU pipeline failed; inspect its diagnostic')
        if time.time()>deadline:raise RuntimeError('Input checkpoint wait timeout')
        time.sleep(10)
    return path

completed=[]
with threadpool_limits(1):
    while not (ROOT/'resolution-amendment.json').exists():time.sleep(10)
    for regime,seeds in [('smooth',P['smooth_seeds'][:2]),('forced_nu004',P['forced_seeds'][:1]),('forced_nu002',P['forced_seeds'][:1])]:
        n=32 if regime=='smooth' else max(64,min(effective_regimes()[regime]['grids']))
        # Amendment is produced before any production checkpoint exists.
        while not (ROOT/'resolution-amendment.json').exists():time.sleep(10)
        for seed in seeds:
            wait_file(filename(seed,regime,'coincident',n))
            for tight in [False,True]:
                r=check(seed,regime,tight);completed.append(dict(seed=seed,regime=regime,tight=tight,kind='full'))
                dump(ROOT/'results'/'cpu-pipeline-progress.json',dict(completed=completed,planned=12))
    for regime,seed in [('smooth',P['smooth_seeds'][0]),('forced_nu002',P['forced_seeds'][0])]:
        n=32 if regime=='smooth' else max(64,min(effective_regimes()[regime]['grids']))
        dt=effective_regimes()[regime]['dt'];wait_file(ROOT/'common-probe-data'/f'{regime}-coincident-s{seed}-n{n}-dt{dt:g}.npz')
        for tight in [False,True]:
            r=common(seed,regime,tight);completed.append(dict(seed=seed,regime=regime,tight=tight,kind='common'))
            dump(ROOT/'results'/'cpu-pipeline-progress.json',dict(completed=completed,planned=12))
dump(ROOT/'results'/'cpu-pipeline-complete.json',dict(status='complete',completed=completed,controller_sha256=digest(ROOT/'cpu_pipeline.py')))
