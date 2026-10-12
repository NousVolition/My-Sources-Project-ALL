"""Compare optimized contact to exhaustive pairs on actual hard-forcing states."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
from contact_hug import Mesh as Exhaustive

def main():
    checks=[];timings=[]
    for name in ['mesh24-one','mesh72-one','hard72-one','stress128-one']:
        r=json.loads((ROOT/'parents'/(name+'.json')).read_text());raw=np.load(ROOT/'parents'/(name+'.npz'));m=Mesh(Parameters(**r['parameters']))
        for i in np.unique(np.linspace(0,len(raw['time'])-1,4).astype(int)):
            y=raw['state'][i];x=m.unpack(y)[0]
            exact=Exhaustive.self_gap(m,x)[0];bound=m.self_gap(x,False)
            near=exact<m.p.self_width
            assert np.max(bound-exact)<1e-12 and np.max(abs(bound[near]-exact[near]),initial=0)<1e-12
            started=time.perf_counter();fast=m.rhs(0,y);fast_seconds=time.perf_counter()-started
            optimized=m.self_gap;m.self_gap=lambda x,gradient=True:Exhaustive.self_gap(m,x,gradient)
            started=time.perf_counter();slow=m.rhs(0,y);slow_seconds=time.perf_counter()-started
            m.self_gap=optimized
            error=float(np.max(abs(fast-slow)/np.maximum(1,abs(slow))))
            assert error<1e-12
            checks.append(dict(run=name,time=float(raw['time'][i]),rhs_scaled_error=error,lower_bound=True,passed=True))
            timings.append(dict(fast=fast_seconds,exhaustive=slow_seconds))
    result=dict(all_passed=True,checks=checks,timings=timings,aggregate_speed_ratio=sum(t['exhaustive'] for t in timings)/sum(t['fast'] for t in timings),test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'data/contact-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
