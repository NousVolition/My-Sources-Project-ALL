"""Exact direct and CUDA observations with integer Fourier wave vectors."""
import json,numpy as np
from corrected_solver import Solver,random_field
from corrected_observation import FastSolver
from final_experiment import ROOT,dump,digest
def validate():
    rows=[]
    for n in [32,64,80,96]:
        a=Solver(n,.04,'gpu');b=FastSolver(n,.04,'gpu');xp=a.xp
        h=xp.stack([random_field(a,80000),random_field(a,80001)])
        x=xp.asarray(np.random.default_rng(7).uniform(0,2*np.pi,(2,10,3)));f=xp.asarray(np.tile(np.eye(3),(2,10,1,1)))
        v,j=a.sample(h,x,10);vv,jj=b.sample(h,x,10)
        aa=a.step(h,x,f,.002);bb=b.step(h,x,f,.002)
        errors=[float(a.cpu(xp.max(abs(y-z)))) for y,z in zip(aa,bb)]
        row=dict(n=n,velocity=float(a.cpu(xp.max(abs(v-vv)))),gradient=float(a.cpu(xp.max(abs(j-jj)))),coupled_errors=errors)
        row['passed']=row['velocity']<1e-12 and row['gradient']<1e-12 and max(errors)<1e-12
        rows.append(row);print('corrected observation '+json.dumps(row),flush=True)
        del a,b,h,x,f,v,j,vv,jj,aa,bb;xp.get_default_memory_pool().free_all_blocks()
    result=dict(all_passed=all(x['passed'] for x in rows),rows=rows,sources={f:digest(ROOT/f) for f in ['validate_corrected.py','corrected_solver.py','corrected_observation.py']})
    dump(ROOT/'results'/'corrected-validation.json',result)
    if not result['all_passed']:raise RuntimeError('Corrected observation validation failed')
    return result
if __name__=='__main__':validate()
