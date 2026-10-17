"""Compare direct CPU and CUDA observations on actual refined fluid fields.

This validates the new 128^3 observation size, using one saved start and endpoint
field. It does not measure a history contrast or select physical effects.
"""
import json,time,traceback
import numpy as np
from final_experiment import ROOT,P,dump,digest
from corrected_solver import Solver
from corrected_observation import FastSolver
from analysis import read,filename

def main():
    source=filename(82000,'forced_nu002','coincident',128,.002)
    deadline=time.monotonic()+172800
    while not (source.with_suffix('.json').exists() and (ROOT/'results'/'refinement-complete.json').exists()):
        for name in ['gpu-pipeline-failure.json','refinement-failure.json']:
            if (ROOT/'results'/name).exists():raise RuntimeError('Source simulation failed: '+name)
        if time.monotonic()>deadline:raise RuntimeError('Refined observation input wait timed out')
        time.sleep(10)
    a,m=read(source,['start_fields','end_fields','modes'])
    reference=Solver(128,.02,'cpu');fast=FastSolver(128,.02,'gpu');xp=fast.xp
    np.testing.assert_array_equal(a['modes'],reference.modes_cpu)
    positions=np.random.default_rng(7).uniform(0,2*np.pi,(1,10,3));rows=[]
    for key in ['start_fields','end_fields']:
        coeff=a[key][:1];v,j=reference.sample(reference.unpack(coeff),positions,10)
        vv,jj=fast.sample(fast.unpack(coeff),xp.asarray(positions),10)
        errors=dict(velocity=float(np.max(abs(v-fast.cpu(vv)))),gradient=float(np.max(abs(j-fast.cpu(jj)))))
        rows.append(dict(field=key,**errors,passed=max(errors.values())<1e-11))
    result=dict(all_passed=all(r['passed'] for r in rows),n=128,rows=rows,tolerance=1e-11,
        source_raw_sha256=m['raw_sha256'],scope='Observation validation only; one prepared and endpoint field, no history effect estimate',
        sources={f:digest(ROOT/f) for f in ['validate_refined_observation.py','corrected_solver.py','corrected_observation.py']})
    dump(ROOT/'results'/'refined-observation-validation.json',result)
    print(json.dumps(result),flush=True)
    xp.get_default_memory_pool().free_all_blocks()

if __name__=='__main__':
    try:main()
    except Exception:
        dump(ROOT/'results'/'refined-observation-failure.json',dict(status='failed',traceback=traceback.format_exc()));raise
