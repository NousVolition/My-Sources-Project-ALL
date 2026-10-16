"""Equivalent one-way system with an independently solved common Lorenz input.

Only material and work variables are integrated in this solver. Saved records
combine their actual accepted states with the independently computed Lorenz
state. No material position, rest length or work value is projected or clipped.
"""
import time,numpy as np
from scipy.integrate import BDF,Radau
from mesh_hug import Mesh,swept_shape
from jacobian import coupled_jacobian
from reference_drive import reference
from work_blocks import install
def integrate(p,end,max_step,initial,state,method,rtol,atol,record_step,mesh_type,max_rhs,max_steps,max_seconds,start_time):
 if p.feedback or not 0<=start_time<end:raise ValueError('Require a one-way case and increasing clock')
 if state is None and start_time!=0:raise ValueError('A nonzero clock requires a checkpoint')
 started=time.perf_counter();m=mesh_type(p);horizon=max(.5,end);drive=reference(horizon,initial,p.rho);keep=np.r_[np.arange(m.size-7),np.arange(m.size-4,m.size)]
 def full(t,z):return np.r_[z[:-4],drive.sol(t),z[-4:]]
 first=m.pack(initial) if state is None else np.array(state,copy=True)
 if first.shape!=(m.size,) or not np.isfinite(first).all():raise ValueError('Invalid checkpoint state')
 if not np.array_equal(first[-7:-4],drive.sol(start_time)):raise ValueError('Checkpoint does not exactly match the common Lorenz input')
 z=first[keep];t=float(start_time);times=[t];states=[first.copy()];calls=accepted=rejects=0;failure=None;solver=None
 cc=float(np.min(m.clearance(m.unpack(states[0])[0],False)));sc=float(np.min(m.self_gap(m.unpack(states[0])[0],False)));minimum=dict(clearance=cc,self_clearance=sc,area_ratio=1.,quality_ratio=1.,altitude_ratio=1.,step=max_step)
 area,q,h=swept_shape(m,first,first)
 minimum.update(area_ratio=area,quality_ratio=q,altitude_ratio=h)
 forward=returned=0.;first_step=min(1e-5,max_step,end-t)
 def rhs(t,z):
  nonlocal calls
  calls+=1
  if calls>max_rhs:raise RuntimeError('Declared RHS evaluation budget exceeded')
  if time.perf_counter()-started>max_seconds:raise RuntimeError('Declared wall-time budget exceeded')
  return m.rhs(t,full(t,z))[keep]
 def jac(t,z):
  J=coupled_jacobian(m,m.rhs,t,full(t,z));return J[keep][:,keep].tocsc()
 previous_power=float(rhs(t,z)[-1]);cls=Radau if method=='Radau' else BDF
 first_record=(np.floor(start_time/record_step+1e-9)+1)*record_step
 targets=np.r_[np.arange(first_record,end-1e-12,record_step),end]
 for target in targets:
  while t<target-1e-12:
   if accepted>=max_steps:failure=dict(time=t,message='Declared accepted-step budget exceeded');break
   old_t=t;old_z=z.copy();attempt=min(max_step,target-t)
   if first_step<1e-11:failure=dict(time=t,message='Restart step below minimum',step=first_step);break
   try:
    if solver is None:
     solver=cls(rhs,t,z,end,max_step=attempt,first_step=min(first_step,attempt,end-t),rtol=rtol,atol=atol,jac=jac);install(solver,m)
    solver.max_step=attempt;message=solver.step()
    if solver.status=='failed':raise FloatingPointError(message or 'Implicit solver failed')
    candidate=solver.y.copy();new_t=float(solver.t)
    if not np.isfinite(candidate).all():raise FloatingPointError('Nonfinite state')
    power=float(rhs(new_t,candidate)[-1]);old_y=full(old_t,old_z);new_y=full(new_t,candidate)
    clear,_,self_clear=m.swept_clearance(old_y,new_y);area,q,h=swept_shape(m,old_y,new_y)
    if min(clear,self_clear)<=0 or area<=p.area_floor or q<=p.shape_floor:raise FloatingPointError('Swept material/contact certificate failed')
   except FloatingPointError:
    rejects+=1;attempted_h=float(solver.h_abs) if solver is not None else first_step;first_step=min(first_step,attempt,attempted_h)*.5;solver=None;continue
   except RuntimeError as error:failure=dict(time=t,message=str(error));break
   t=new_t;z=candidate;step=t-old_t;accepted+=1
   if accepted%500==0:print('integration progress',round(t,6),'steps',accepted,'calls',calls,'restarts',rejects,flush=True)
   forward+=step/2*(max(previous_power,0)+max(power,0));returned+=step/2*(max(-previous_power,0)+max(-power,0));previous_power=power
   for key,value in [('clearance',clear),('self_clearance',self_clear),('area_ratio',area),('quality_ratio',q),('altitude_ratio',h),('step',step)]:minimum[key]=min(minimum[key],value)
   first_step=min(max_step,max(1e-11,step))
  if failure:
   if t>times[-1]+1e-12:times.append(t);states.append(full(t,z))
   break
  assert abs(t-target)<1e-10
  times.append(t);states.append(full(t,z))
 audit=dict(method=method,accepted_steps=accepted,restarted_trials=rejects,rhs_calls=calls,minimum_swept_clearance=minimum['clearance'],minimum_swept_self_clearance=minimum['self_clearance'],minimum_swept_area_ratio=minimum['area_ratio'],minimum_swept_quality_ratio=minimum['quality_ratio'],minimum_swept_altitude_ratio=minimum['altitude_ratio'],minimum_step=minimum['step'],forward_work_step_trapezoid=float(forward),return_work_step_trapezoid=float(returned),rtol=rtol,atol=atol,max_step=max_step,max_rhs=max_rhs,max_steps=max_steps,max_seconds=max_seconds,seconds=time.perf_counter()-started,split_one_way=True,lorenz_method='DOP853',lorenz_rtol=1e-12,lorenz_atol=1e-14,lorenz_max_step=.002,lorenz_reference_end=horizon,lorenz_rhs_calls=drive.nfev,work_block_factorization=True,checkpoint_clock_preserved=True)
 return m,np.array(times),np.array(states),failure,audit
