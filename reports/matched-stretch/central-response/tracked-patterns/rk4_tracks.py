"""RK4 diagnostic tracks in stored original velocity fields, not a new fluid run."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy import fft
from analyze import trilinear,save,sha,wrap
ROOT=Path(__file__).resolve().parent
FOLDER=ROOT.parents[1]/'runs/baseline-n64-base'
def rk4(x,dt,velocity,t):
    k1=velocity(x,t);k2=velocity(x+dt*k1/2,t+dt/2);k3=velocity(x+dt*k2/2,t+dt/2);k4=velocity(x+dt*k3,t+dt)
    return x+dt*(k1+2*k2+2*k3+k4)/6
def tracks():
    r=json.loads((FOLDER/'result.json').read_text());assert r['status']=='complete'
    seed=np.array(np.meshgrid([-.15,0,.15],[-.15,0,.15],[-.3,0,.3],indexing='ij')).reshape(3,-1).T
    result={'initial_positions':seed.tolist(),'labels':27,'source':'baseline-n64-base','domain_side':6,'external_force':0,'fluid_evolution_changed':False,'method':'RK4 passive tracks with periodic trilinear spatial and linear temporal interpolation of stored velocity','runs':{},'verified_fields':[]}
    cache={}
    def field(index):
        if index not in cache:
            p=FOLDER/f'field-{index:03d}.npy';h=np.load(p,allow_pickle=False);assert np.isfinite(h).all()
            cache[index]=fft.irfftn(h,s=(64,)*3,axes=(-3,-2,-1),workers=1)
            if index not in checked:result['verified_fields'].append({'file':p.name,'sha256':sha(p)});checked.add(index)
        return cache[index]
    checked=set()
    for name,stride,substeps in [('base',1,4),('half_step',1,8),('quarter_step',1,16),('eighth_step',1,32),('coarser_saved_times',2,64)]:
        x=seed.copy();frames=[{'t':0.,'positions':x.tolist()}];cache.clear()
        for start in range(0,40,stride):
            u0=field(start);u1=field(start+stride);duration=.01*stride;dt=duration/substeps
            def velocity(points,time):
                fraction=np.clip(time/duration,0,1)
                return (1-fraction)*trilinear(u0,points)+fraction*trilinear(u1,points)
            for step in range(substeps):x=rk4(x,dt,velocity,step*dt)
            frames.append({'t':round(.01*(start+stride),8),'positions':x.tolist()})
            cache.pop(start,None)
        result['runs'][name]={'dt':dt,'saved_velocity_interval':duration,'frames':frames}
        print('RK4 '+name+' complete',flush=True)
    fine=result['runs']['eighth_step']['frames'];base=result['runs']['quarter_step']['frames'];coarse=result['runs']['coarser_saved_times']['frames']
    err_step=np.linalg.norm(wrap(np.array([f['positions'] for f in base])-np.array([f['positions'] for f in fine])),axis=2)
    err_samples=np.linalg.norm(wrap(np.array([f['positions'] for f in coarse])-np.array([f['positions'] for f in fine[::2]])),axis=2)
    # Independent known-solution control for the RK4 time integration.
    initial=np.array([[1.,0.,0.]]);errors=[]
    for dt in (.1,.05):
        x=initial.copy()
        for k in range(round(1/dt)):x=rk4(x,dt,lambda p,t:np.stack([-p[:,1],p[:,0],p[:,2]*0],axis=1),k*dt)
        errors.append(float(np.linalg.norm(x-[np.cos(1),np.sin(1),0])))
    assert errors[0]/errors[1]>15 and np.isfinite(err_step).all()
    ladder=[]
    for a,b in zip(['base','half_step','quarter_step'],['half_step','quarter_step','eighth_step']):
        aa=result['runs'][a];bb=result['runs'][b];dif=np.linalg.norm(wrap(np.array([f['positions'] for f in aa['frames']])-np.array([f['positions'] for f in bb['frames']])),axis=2)
        ladder.append({'dt':aa['dt'],'half_dt':bb['dt'],'max_position_difference':float(dif.max())})
    result['verification']={'step_ladder':ladder,'step_halving_max_position_difference':float(err_step.max()),'saved_time_coarsening_max_position_difference':float(err_samples.max()),'through_0_1_step_difference':float(err_step[:11].max()),'through_0_1_saved_time_difference':float(err_samples[:6].max()),'spatial_cell_size':6/64,'rotation_exact_solution_errors':errors,'rotation_error_reduction':errors[0]/errors[1]}
    result['limits']=['Traces sample an already computed velocity field; RK4 does not repair its unresolved spin peak.','Only the initial central neighborhood is seeded; these 27 labels do not represent every molecule.','Physical x-y plots project 3D time-dependent trajectories; projected crossings need not be collisions.','Small RK4 error does not bound errors from spatial or saved-time interpolation.']
    save(ROOT/'rk4-tracks.json',result);print(json.dumps(result['verification']))
if __name__=='__main__':tracks()
