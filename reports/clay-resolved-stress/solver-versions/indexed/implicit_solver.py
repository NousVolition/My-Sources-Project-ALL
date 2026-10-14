"""Guarded implicit integration of the unchanged work-balanced ODE.

Only positive-domain stages and certified straight accepted paths are retained.
Failed stages restart from the preceding state with a smaller trial step.
This wrapper does not certify an implicit dense-output polynomial.
"""
import time
import numpy as np
from scipy.integrate import Radau,BDF
from mesh_hug import Mesh,swept_shape
from jacobian import coupled_jacobian

def integrate(p,end=4,max_step=.005,initial=(-8.,8.,27.),state=None,method='Radau',
              rtol=2e-7,atol=2e-9,record_step=.04,mesh_type=Mesh,max_rhs=4000000,max_steps=200000,max_seconds=900,start_time=0.):
    if method not in ['Radau','BDF']:raise ValueError('Unknown implicit method')
    if end<=start_time:raise ValueError('End must exceed the starting time')
    m=mesh_type(p);y=m.pack(initial) if state is None else np.array(state,copy=True)
    m.energy_force(*m.unpack(y)[:2]);started=time.perf_counter()
    t=float(start_time);times=[t];states=[y.copy()];calls=accepted=rejects=0;failure=None;solver=None
    floor=p.area_floor if isinstance(m,Mesh) else .01;qfloor=p.shape_floor if isinstance(m,Mesh) else 0.
    cc=float(np.min(m.clearance(m.unpack(y)[0],False)));sc=float(np.min(m.self_gap(m.unpack(y)[0],False)))
    minimum=dict(clearance=cc,self_clearance=sc,area_ratio=1.,quality_ratio=1.,altitude_ratio=1.,step=max_step)
    forward=returned=0.;first_step=min(1e-5,max_step,end)
    def rhs(t,y):
        nonlocal calls
        calls+=1
        if calls>max_rhs:raise RuntimeError('Declared RHS evaluation budget exceeded')
        if time.perf_counter()-started>max_seconds:raise RuntimeError('Declared wall-time budget exceeded')
        return m.rhs(t,y)
    previous_power=float(rhs(t,y)[-1]);cls=Radau if method=='Radau' else BDF
    targets=np.r_[np.arange((np.floor(start_time/record_step+1e-9)+1)*record_step,end-1e-12,record_step),end]
    for target in targets:
        while t<target-1e-12:
            if accepted>=max_steps:
                failure=dict(time=t,message='Declared accepted-step budget exceeded');break
            old_t=t;old_y=y.copy();attempt=min(max_step,target-t)
            if first_step<1e-11:
                failure=dict(time=t,message='Restart step below minimum',step=first_step);break
            try:
                if solver is None:
                    solver=cls(rhs,t,y,end,max_step=attempt,first_step=min(first_step,attempt,end-t),rtol=rtol,atol=atol,jac=lambda a,b:coupled_jacobian(m,rhs,a,b))
                solver.max_step=attempt
                message=solver.step()
                if solver.status=='failed':raise FloatingPointError(message or 'Implicit solver failed')
                candidate=solver.y.copy();new_t=float(solver.t)
                if not np.isfinite(candidate).all():raise FloatingPointError('Nonfinite state')
                # RHS validates the same energy domain and supplies power; avoid
                # recomputing its complete force once more on accepted states.
                power=float(rhs(new_t,candidate)[-1])
                clear,_,self_clear=m.swept_clearance(old_y,candidate)
                area,q,h=swept_shape(m,old_y,candidate)
                if min(clear,self_clear)<=0 or area<=floor or q<=qfloor:
                    raise FloatingPointError('Swept material/contact certificate failed')
            except FloatingPointError:
                rejects+=1
                attempted_h=float(solver.h_abs) if solver is not None else first_step
                first_step=min(first_step,attempt,attempted_h)*.5;solver=None;continue
            except RuntimeError as error:
                failure=dict(time=t,message=str(error));break
            t=new_t;y=candidate;step=t-old_t;accepted+=1
            if accepted%500==0:print('integration progress',round(t,6),'steps',accepted,'calls',calls,'restarts',rejects,flush=True)
            forward+=step/2*(max(previous_power,0)+max(power,0))
            returned+=step/2*(max(-previous_power,0)+max(-power,0));previous_power=power
            for key,value in [('clearance',clear),('self_clearance',self_clear),('area_ratio',area),('quality_ratio',q),('altitude_ratio',h),('step',step)]:
                minimum[key]=min(minimum[key],value)
            first_step=min(max_step,max(1e-11,step))
        if failure:
            if t>times[-1]+1e-12:times.append(t);states.append(y.copy())
            break
        assert abs(t-target)<1e-10,(t,target)
        times.append(t);states.append(y.copy())
    audit=dict(method=method,accepted_steps=accepted,restarted_trials=rejects,rhs_calls=calls,
        minimum_swept_clearance=minimum['clearance'],minimum_swept_self_clearance=minimum['self_clearance'],
        minimum_swept_area_ratio=minimum['area_ratio'],minimum_swept_quality_ratio=minimum['quality_ratio'],
        minimum_swept_altitude_ratio=minimum['altitude_ratio'],minimum_step=minimum['step'],
        forward_work_step_trapezoid=float(forward),return_work_step_trapezoid=float(returned),
        rtol=rtol,atol=atol,max_step=max_step,max_rhs=max_rhs,max_steps=max_steps,max_seconds=max_seconds,seconds=time.perf_counter()-started)
    return m,np.array(times),np.array(states),failure,audit

