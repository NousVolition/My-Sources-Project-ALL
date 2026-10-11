import json
from experiment import P,ROOT,run_one,dump
from analysis import pilot_audit

completed=[]
for regime,cfg in P['regimes'].items():
    for seed in P['pilot_seeds']:
        for n in cfg['grids']:
            completed.append(run_one(seed,regime,'coincident',n)['config'])
            dump(ROOT/'results'/'pilot-execution.json',dict(planned=12,completed=completed))
        audit=pilot_audit();print('pilot screens: '+json.dumps(audit['pairs'][-1]),flush=True)
audit=pilot_audit();print('All pilot gates: '+str(audit['all_passed']),flush=True)
