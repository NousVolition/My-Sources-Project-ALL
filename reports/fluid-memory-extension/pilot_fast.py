import json
from accelerated_experiment import P,ROOT,run_one,dump
from analysis import pilot_audit,read,response
import numpy as np
completed=[]
for regime,cfg in P['regimes'].items():
    for seed in P['pilot_seeds']:
        for n in cfg['grids']:
            completed.append(run_one(seed,regime,'coincident',n)['config'])
            dump(ROOT/'results'/'pilot-fast-execution.json',dict(planned=12,completed=completed))
        audit=pilot_audit();print('pilot screens: '+json.dumps(audit['pairs'][-1]),flush=True)
audit=pilot_audit();print('All pilot gates: '+str(audit['all_passed']),flush=True)
# Numerical equivalence with every original completed pilot configuration.
checks=[]
for original in (ROOT/'data').glob('*.npz'):
    if not original.with_suffix('.json').exists():continue
    fast=ROOT/'fast-data'/original.name
    if not fast.exists():continue
    a,ma=read(original,['time','tangent','positions','end_fields']);b,mb=read(fast,['time','tangent','positions','end_fields'])
    ra,rb=response(a,ma),response(b,mb)
    error=max(float(abs(ra[m][arm]-rb[m][arm]).max()) for m in P['probe_methods'] for arm in P['arms'])
    checks.append(dict(file=original.name,response_error=error,position_error=float(abs(a['positions']-b['positions']).max()),field_error=float(abs(a['end_fields']-b['end_fields']).max()),passed=error<1e-10))
dump(ROOT/'results'/'accelerated-pilot-equivalence.json',dict(all_passed=all(v['passed'] for v in checks),checks=checks))
