import json,time,numpy as np
from reference_solver import Solver,random_field
from fast_observation import FastSolver
from experiment import ROOT,dump,digest

results=[]
for n in [24,32,64,96]:
    old=Solver(n,.04,'gpu');new=FastSolver(n,.04,'gpu');xp=old.xp
    h=xp.stack([random_field(old,81000),random_field(old,81001)])
    points=xp.asarray(np.random.default_rng(3).uniform(0,2*np.pi,(2,10,3)));f=xp.asarray(np.tile(np.eye(3),(2,10,1,1)))
    v,j=old.sample(h,points,10);vv,jj=new.sample(h,points,10)
    row=dict(n=n,velocity_error=float(old.cpu(xp.max(abs(v-vv)))),gradient_error=float(old.cpu(xp.max(abs(j-jj)))))
    start=time.time()
    for _ in range(10):new.sample(h,points,10)
    xp.cuda.get_current_stream().synchronize();row['fast_seconds_per_sample']=(time.time()-start)/10
    start=time.time()
    for _ in range(10):old.sample(h,points,10)
    xp.cuda.get_current_stream().synchronize();row['reference_seconds_per_sample']=(time.time()-start)/10
    hh,xx,ff,bb,cc=old.step(h,points,f,.002);hhh,xxx,fff,bbb,ccc=new.step(h,points,f,.002)
    row['field_step_error']=float(old.cpu(xp.max(abs(hh-hhh))))
    row['position_step_error']=float(old.cpu(xp.max(abs(xx-xxx))))
    row['tangent_step_error']=float(old.cpu(xp.max(abs(ff-fff))))
    row['pass']=row['velocity_error']<1e-12 and row['gradient_error']<1e-12 and row['field_step_error']==0 and row['position_step_error']<1e-13 and row['tangent_step_error']<1e-13
    results.append(row);print(json.dumps(row),flush=True);del h,points,f;xp.get_default_memory_pool().free_all_blocks()
dump(ROOT/'results'/'fast-observation-validation.json',dict(all_passed=all(r['pass'] for r in results),checks=results,sources={f:digest(ROOT/f) for f in ['reference_solver.py','fast_observation.py','validate_fast.py']}))
