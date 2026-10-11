"""Compare CPU and GPU fluid RHS, direct observations and joint trajectories."""
import json
from pathlib import Path
import numpy as np
from solver import Solver, random_field, curl_gaussian
from run import digest, dump, marker_cloud, P

def main():
    checks=[]
    for seed in P['pilot_seeds']:
        c,g=Solver(16,.12,'cpu'),Solver(16,.12,'gpu')
        h=np.stack([random_field(c,seed)]*2)
        p=np.stack([marker_cloud(seed)]*2)
        f=np.tile(np.eye(3),(2,P['centers'],1,1))
        hc,pc,fc=h.copy(),p.copy(),f.copy()
        hg,pg,fg=g.xp.asarray(h),g.xp.asarray(p),g.xp.asarray(f)
        uv,jv=c.sample(h,p,P['centers']);ug,jg=g.sample(hg,pg,P['centers'])
        field_sample_error=float(abs(uv-g.cpu(ug)).max());gradient_error=float(abs(jv-g.cpu(jg)).max())
        rc,_,_=c.rhs(h,np.zeros_like(h));rg,_,_=g.rhs(hg,g.xp.zeros_like(hg))
        rhs_error=float(np.sqrt(c.inner(rc-g.cpu(rg),rc-g.cpu(rg)).max()/c.inner(rc,rc).max()))
        for _ in range(20):
            hc,pc,fc,_,_=c.step(hc,pc,fc,.01)
            hg,pg,fg,_,_=g.step(hg,pg,fg,.01)
        field_error=float(np.sqrt(c.inner(hc-g.cpu(hg),hc-g.cpu(hg)).max()/c.inner(hc,hc).max()))
        position_error=float(abs(pc-g.cpu(pg)).max());tangent_error=float(abs(fc-g.cpu(fg)).max())
        checks.append(dict(seed=seed,sample_error=field_sample_error,gradient_error=gradient_error,rhs_relative_error=rhs_error,
                           field_relative_error=field_error,position_error=position_error,tangent_error=tangent_error,
                           passed=bool(max(field_sample_error,gradient_error,rhs_error,field_error,position_error,tangent_error)<1e-11)))
    result=dict(checks=checks,all_passed=all(c['passed'] for c in checks),sources={p:digest(Path(__file__).parent/p) for p in ['solver.py','validate_backend.py']})
    dump(Path(__file__).parent/'backend-validation.json',result)
    print(json.dumps(result,indent=2),flush=True)
    if not result['all_passed']:raise SystemExit(1)

if __name__=='__main__':main()
