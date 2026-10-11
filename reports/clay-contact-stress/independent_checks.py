"""Independent adaptive integrator and time/tolerance refinement checks."""
from pathlib import Path
from dataclasses import replace
import json,time,sys,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from contact_hug import *
from scipy.integrate import trapezoid
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def main():
    checks=[]
    for key,p,end in [('strong',Parameters(preload=.35,alpha=8),2),('highrho',Parameters(preload=.35,alpha=4,rho=160),2),('self-repair',Parameters(preload=.35,alpha=4,feedback=False),6)]:
        m=Mesh(p);t=np.linspace(0,end,round(end/.01)+1);start=time.perf_counter()
        path=DATA/('independent-'+key+'.npz');meta=DATA/('independent-'+key+'.json')
        digest=hashlib.sha256((ROOT/'contact_hug.py').read_bytes()).hexdigest()
        if path.exists() and meta.exists():
            assert json.loads(meta.read_text())['model_sha256']==digest
            y=np.load(path)['state']
        else:
            sol=solve_ivp(m.rhs,(0,end),m.pack((-8,8,27)),method='DOP853',rtol=2e-10,atol=2e-12,max_step=.0005,t_eval=t)
            assert sol.success
            y=sol.y.T
            np.savez_compressed(path,time=t,state=y)
            meta.write_text(json.dumps(dict(model_sha256=digest,method='DOP853',rtol=2e-10,atol=2e-12,max_step=.0005,parameters=p.__dict__,initial=[-8,8,27],end=end,record_step=.01,success=bool(sol.success),seconds=time.perf_counter()-start),indent=2)+'\n')
            print(key,'independent integration recorded',flush=True)
        if '--record-only' in sys.argv:continue
        reference=y[::4]
        if key=='self-repair':
            fine=np.load(DATA/'amplitude4-one.npz')['state'];pairs=[('independent-DOP853',fine,reference)]
        else:
            fine=np.load(DATA/('refine-'+key+'-fine.npz'))['state'];base=np.load(DATA/('refine-'+key+'-base.npz'))['state']
            pairs=[('step-and-tolerance',base,fine),('independent-DOP853',fine,reference)]
        scale=np.maximum(1,np.max(abs(reference[:,:-4]),axis=0))
        for name,a,b in pairs:
            error=float(np.max(abs(a[:,:-4]-b[:,:-4])/scale))
            checks.append(dict(name=key+'-'+name,error=error,target=1e-4,passed=bool(error<1e-4)))
        checks.append(dict(name=key+'-independent-domain',minimum_clearance=min(m.diagnostics(v)['certified_clearance'] for v in y),passed=bool(all(m.diagnostics(v)['certified_clearance']>0 for v in y)),seconds=time.perf_counter()-start))
        def stats(tt,states):
            use=tt>=end/2;tt=tt[use];ly=states[use,-7:-4];mean=lambda a:trapezoid(a,tt)/(tt[-1]-tt[0]);zm=mean(ly[:,2])
            return np.array([np.sqrt(mean(np.sum(ly[:,:2]**2,axis=1))),np.sqrt(mean((ly[:,2]-zm)**2))])
        a=stats(t,y);b=stats(t[::4],reference);err=float(np.max(abs(a-b)/np.maximum(abs(a),1e-10)))
        checks.append(dict(name=key+'-output-sampling',error=err,target=.01,passed=bool(err<.01)))
    if '--record-only' in sys.argv:return
    a=np.load(DATA/'zero-port-one.npz')['state'];b=np.load(DATA/'zero-port-two.npz')['state']
    error=float(np.max(abs(a-b)));checks.append(dict(name='zero-port-one-two-identical',error=error,target=1e-12,passed=bool(error<1e-12)))
    result=dict(checks=checks,all_passed=all(c['passed'] for c in checks))
    (DATA/'solver-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    if not result['all_passed']:raise SystemExit(1)
if __name__=='__main__':main()
