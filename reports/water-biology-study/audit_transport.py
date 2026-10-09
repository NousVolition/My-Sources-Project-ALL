"""Check transient scalar convergence, not only the near-uniform final field."""
from pathlib import Path
import json
import numpy as np
from core import scalar_run
from run_study import write_csv

ROOT=Path(__file__).resolve().parent


def main():
    runs={key:scalar_run(n=key[0],dt=key[1]) for key in [(32,.01),(64,.01),(128,.01),(128,.02),(128,.005)]}
    rows=[]; errors={}
    for a,b,label in [((32,.01),(64,.01),'grid32_64'),((64,.01),(128,.01),'grid64_128'),
                       ((128,.02),(128,.01),'dt02_01'),((128,.01),(128,.005),'dt01_005')]:
        ra,ma=runs[a]; rb,mb=runs[b]
        ratio=b[0]//a[0]
        e=np.sqrt(np.mean((ma-mb[:,::ratio,::ratio])**2,axis=(1,2)))
        errors[label]=float(e.max())
        rows.extend([[label,float(t),float(v)] for t,v in zip(ra[:,0],e)])
    times={f'n{n}_dt{dt}':float(r[0][np.flatnonzero(r[0][:,3]>=.95)[0],0]) for (n,dt),r in runs.items()}
    # Predefined numeric screens for representative seed 100.
    assert errors['grid64_128']<1e-4
    assert errors['dt01_005']<1e-3
    assert abs(times['n64_dt0.01']-times['n128_dt0.005'])<=.10000001
    result={'seed':100,'maximum_over_saved_times_RMS':errors,'t95_s':times,
            'dt_order_from_max_error':float(np.log2(errors['dt02_01']/errors['dt01_005'])),
            'status':'passed; sampled transient field comparison, not all initial conditions'}
    write_csv(ROOT/'data/scalar_transient_convergence.csv',['comparison','t','RMS'],rows)
    (ROOT/'transport-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
