"""Independent adaptive integrator and time/tolerance refinement checks."""
from pathlib import Path
from dataclasses import replace
import json,time,sys,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from contact_hug import *
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def main():
    checks=[]
    for key,p in [('strong',Parameters(preload=.35,alpha=8)),('highrho',Parameters(preload=.35,alpha=4,rho=160))]:
        m=Mesh(p);t=np.linspace(0,4,101);start=time.perf_counter()
        path=DATA/('independent-'+key+'.npz');meta=DATA/('independent-'+key+'.json')
        digest=hashlib.sha256((ROOT/'contact_hug.py').read_bytes()).hexdigest()
        if path.exists() and meta.exists():
            assert json.loads(meta.read_text())['model_sha256']==digest
            y=np.load(path)['state']
        else:
            sol=solve_ivp(m.rhs,(0,4),m.pack((-8,8,27)),method='DOP853',rtol=2e-10,atol=2e-12,max_step=.0005,t_eval=t)
            assert sol.success
            y=sol.y.T
            np.savez_compressed(path,time=t,state=y)
            meta.write_text(json.dumps(dict(model_sha256=digest,method='DOP853',rtol=2e-10,atol=2e-12,max_step=.0005,parameters=p.__dict__,initial=[-8,8,27],end=4,success=bool(sol.success),seconds=time.perf_counter()-start),indent=2)+'\n')
            print(key,'independent integration recorded',flush=True)
        if '--record-only' in sys.argv:continue
        fine=np.load(DATA/('refine-'+key+'-fine.npz'))['state']
        base=np.load(DATA/('refine-'+key+'-base.npz'))['state']
        scale=np.maximum(1,np.max(abs(y[:,:-4]),axis=0))
        for name,a,b in [('step-and-tolerance',base,fine),('independent-DOP853',fine,y)]:
            error=float(np.max(abs(a[:,:-4]-b[:,:-4])/scale))
            checks.append(dict(name=key+'-'+name,error=error,target=1e-4,passed=bool(error<1e-4)))
        checks.append(dict(name=key+'-independent-domain',minimum_clearance=min(m.diagnostics(v)['certified_clearance'] for v in y),passed=bool(all(m.diagnostics(v)['certified_clearance']>0 for v in y)),seconds=time.perf_counter()-start))
    if '--record-only' in sys.argv:return
    a=np.load(DATA/'zero-port-one.npz')['state'];b=np.load(DATA/'zero-port-two.npz')['state']
    error=float(np.max(abs(a-b)));checks.append(dict(name='zero-port-one-two-identical',error=error,target=1e-12,passed=bool(error<1e-12)))
    result=dict(checks=checks,all_passed=all(c['passed'] for c in checks))
    (DATA/'solver-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    if not result['all_passed']:raise SystemExit(1)
if __name__=='__main__':main()
