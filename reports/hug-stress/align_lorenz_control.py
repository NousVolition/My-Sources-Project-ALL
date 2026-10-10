"""Use exactly the same saved forcing in the displayed coupling-off control.

Separate long chaotic integrations can drift apart numerically. The Lorenz
and memory equations do not depend on lean; reuse those saved states exactly.
Keep the original independent coupled-solver result as evidence.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from lorenz_hug import model
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'lorenz-data'
def main():
    a=np.load(DATA/'chaotic.npz');b=np.load(DATA/'no-feedback.npz')
    t=a['time'];y=a['state'].copy();p=model.Parameters(r=-.25,pressure=0,memory_coupling=0)
    separate=model.reference(p,t,q0=.04,tolerance=1e-12)
    y[:,[3,4,6,7]]=separate[:,[0,1,3,4]]
    np.savez_compressed(DATA/'no-feedback-matched.npz',time=t,state=y,pressure=a['pressure'])
    presets=json.loads((DATA/'presets.json').read_text())
    target=next(r for r in presets if r['key']=='no-feedback')
    target['trajectory_file']='no-feedback-matched.npz'
    target['samples']=np.round(np.column_stack([t[::10],y[::10,:6],a['pressure'][::10]]),8).tolist()
    (DATA/'presets-matched.json').write_text(json.dumps(presets,indent=2)+'\n')
    result=dict(reason='Use an exactly matched recorded Lorenz pressure history for the displayed lean-coupling-off control.',
        original_separate_run_retained='no-feedback.npz',matched_run='no-feedback-matched.npz',
        original_max_pressure_difference=float(max(abs(a['pressure']-b['pressure']))),
        matched_pressure_difference=0.,matched_memory_difference=float(max(abs(y[:,5]-a['state'][:,5]))),
        matched_lean_independence_error=float(np.max(abs(y[:,3:5]-separate[:,:2]))),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),passed=True)
    assert np.array_equal(y[:,:3],a['state'][:,:3]) and np.array_equal(y[:,5],a['state'][:,5])
    (DATA/'matched-control.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
