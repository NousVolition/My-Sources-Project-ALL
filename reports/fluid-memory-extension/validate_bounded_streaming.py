"""Check broad-spectrum input packing and a coupled step across branch batches."""
import json,numpy as np
from final_experiment import ROOT,DrivenSolver,dump,digest
from bounded_final_experiment_v2 import DrivenSolver as BoundedSolver

def validate():
    s=DrivenSolver(128,.02,.2,'gpu');bounded=BoundedSolver(128,.02,.2,'gpu');xp=s.xp
    real=np.random.default_rng(97001).normal(size=(8,3,128,128,128))
    h=s.project(xp.fft.rfftn(xp.asarray(real),axes=(-3,-2,-1)))
    h*=.8/xp.sqrt(s.inner(h,h))[:,None,None,None,None]
    packed=s.pack(h);roundtrip=bounded.unpack(packed)
    input_error=float(s.cpu(xp.max(abs(h-roundtrip))))
    del real,roundtrip
    xp.get_default_memory_pool().free_all_blocks()
    x=xp.asarray(np.random.default_rng(97002).uniform(0,2*np.pi,(8,10,3)))
    f=xp.asarray(np.tile(np.eye(3),(8,10,1,1)))
    reference=s.step(h,x,f,.0005)
    xp.get_default_memory_pool().free_all_blocks()
    pieces=[bounded.step(bounded.unpack(packed[i:i+6]),x[i:i+6],f[i:i+6],.0005) for i in range(0,8,6)]
    errors=[float(s.cpu(xp.max(abs(a-xp.concatenate([r[j] for r in pieces],axis=0))))) for j,a in enumerate(reference)]
    result=dict(passed=input_error==0 and max(errors)<1e-11,n=128,branches=8,comparison_batches=[6,2],production_batch=8,input_seed=97001,position_seed=97002,
        gpu_pool_bytes=xp.get_default_memory_pool().total_bytes(),
        input_roundtrip_max_error=input_error,coupled_step_max_errors=errors,tolerance=1e-11,
        scope='Broad-spectrum real velocity input; exact complex128 packing and one coupled fluid/position/deformation/energy step. Resource-scheduling check, not a fluid memory test.',
        sources={p:digest(ROOT/p) for p in ['validate_bounded_streaming.py','final_experiment.py','bounded_final_experiment_v2.py','corrected_solver.py','corrected_observation.py']})
    dump(ROOT/'results'/'bounded-streaming-validation.json',result);print(json.dumps(result),flush=True)
    if not result['passed']:raise RuntimeError('Packed streaming/batch equivalence failed')
    xp.get_default_memory_pool().free_all_blocks();return result

if __name__=='__main__':validate()
