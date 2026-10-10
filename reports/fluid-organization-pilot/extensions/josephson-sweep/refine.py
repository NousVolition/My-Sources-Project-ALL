"""Follow-up to the failed full-sweep step gate; no original gates are relaxed.

Replay whole histories with adaptive DOP853, then cross-check the sensitive
alpha=.5 downward S=200 history with finer fixed steps and tighter tolerances.
Declared follow-up choices are saved before its simulations start.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from model import to_state,from_state,stable_rhs,stable_velocity,dissipation,integrate_state,voltage,power_residual
ROOT=Path(__file__).resolve().parent


def adaptive(s0,sign,bias,alpha,duration,rtol=1e-11,atol=1e-12,max_step=.5,stride=5):
    shape=s0.shape;size=s0.size;batchshape=shape[:-1]
    def f(t,y):
        state=y[:size].reshape(shape)
        return np.concatenate([stable_rhs(state,bias,alpha).ravel(),
                               dissipation(stable_velocity(state,sign,bias,alpha),alpha).ravel()])
    sol=solve_ivp(f,(0,duration),np.concatenate([s0.ravel(),np.zeros(s0[...,0].size)]),
                  method='DOP853',rtol=rtol,atol=atol,max_step=max_step,t_eval=np.arange(0,duration+stride/2,stride))
    if not sol.success:raise RuntimeError(sol.message)
    return sol.t,sol.y[:size].T.reshape((-1,)+shape),sol.y[size:].T.reshape((-1,)+batchshape)


def main():
    c=json.loads((ROOT/'protocol.json').read_text());out=ROOT/'data'
    followup={'reason':'Original full-sweep fixed-step gate failed at alpha=0.5 in the S=200 downward sweep, despite independent fresh-start checks passing.',
              'adaptive_reference':{'method':'DOP853','rtol':1e-11,'atol':1e-12,'max_step':.5},
              'targeted_fixed_steps':[.0125,.00625],
              'targeted_tighter_reference':{'alpha':.5,'mode':'down','settle':200,'rtol':1e-12,'atol':1e-13,'max_step':.25},
              'voltage_acceptance_limit':2e-5,'reference_tolerance_change_limit':2e-6,
              'window_step_check':'Replay every W=800 reference from its same settled state at RK4 dt=.025; compare total voltage and individual voltages to adaptive reference.'}
    (ROOT/'refinement_protocol.json').write_text(json.dumps(followup,indent=2)+'\n')
    initial=np.load(out/'initial_conditions.npz')['phases'];n=len(initial);alphas=np.array(c['alphas']);biases=np.array(c['biases'])
    init=np.broadcast_to(initial,(3,len(alphas),n,2));init_state,sign=to_state(init)
    aa=np.broadcast_to(alphas[None,:,None],init_state.shape[:-1])
    for settle in [100,200]:
        path=out/f'reference_S{settle}.npz'
        if path.exists():continue
        previous=init_state.copy();starts=[];settled=[];ends=[];means=[];res=[];states_all=[];ds=[];orders=[]
        for j in range(len(biases)):
            previous[0]=init_state[0]
            ib=np.broadcast_to(np.array([biases[j],biases[j],biases[-j-1]])[:,None,None],aa.shape)
            t,states,d=adaptive(previous,sign,ib,aa,settle+200)
            p=from_state(states,sign);idx=settle//5
            starts.append(previous.copy());settled.append(states[idx]);ends.append(states[-1]);means.append(voltage(p[idx],p[-1],200))
            res.append(power_residual(p,d,ib));states_all.append(states);ds.append(d);orders.append(ib[:,0,0])
            previous=states[-1].copy()
            if j%4==0:print(f'Adaptive reference S={settle}, index {j+1}/17',flush=True)
        np.savez_compressed(path,initial=np.array(starts),settled=np.array(settled),endpoint=np.array(ends),
                            voltage=np.array(means),residual=np.array(res),states=np.array(states_all),sign=sign,
                            dissipated=np.array(ds),time=t,bias_order=np.array(orders),alphas=alphas,settle=settle,window=200)
    with np.load(out/'reference_S200.npz') as z:p0=z['settled'];order=z['bias_order'];sign=z['sign']
    ib=np.broadcast_to(order[:,:,None,None],p0.shape[:-1]);aa=np.broadcast_to(alphas[None,None,:,None],p0.shape[:-1])
    for name,method in [('reference_windows','adaptive'),('reference_windows_RK4','fixed')]:
        path=out/f'{name}.npz'
        if path.exists():continue
        if method=='adaptive':t,states,d=adaptive(p0,sign,ib,aa,800)
        else:t,states,d=integrate_state(p0,sign,ib,aa,800,.025,stride=5)
        p=from_state(states,sign)
        np.savez_compressed(path,time=t,states=states,sign=sign,dissipated=d,bias_order=order,alphas=alphas,residual=power_residual(p,d,ib))
        print(name+' done',flush=True)
    s0,sign=to_state(initial)
    for label,dt in [('dt0125',.0125),('dt00625',.00625),('tighter',None)]:
        path=out/f'targeted_down_{label}.npz'
        if path.exists():continue
        previous=s0.copy();means=[];states_all=[];res=[]
        for j,b in enumerate(biases[::-1]):
            if dt is None:t,states,d=adaptive(previous,sign,b,.5,400,rtol=1e-12,atol=1e-13,max_step=.25)
            else:t,states,d=integrate_state(previous,sign,b,.5,400,dt,stride=5)
            p=from_state(states,sign);means.append(voltage(p[40],p[-1],200));states_all.append(states);res.append(power_residual(p,d,b))
            previous=states[-1].copy()
            if j%4==0:print(f'Targeted {label}, index {j+1}/17',flush=True)
        np.savez_compressed(path,voltage=np.array(means),states=np.array(states_all),sign=sign,bias_order=biases[::-1],residual=np.array(res),time=t)
    with np.load(out/'reference_S200.npz') as z:reference=z['voltage'][:,2,1]
    comparisons={};records={}
    for label in ['dt0125','dt00625','tighter']:
        with np.load(out/f'targeted_down_{label}.npz') as z:records[label]=z['voltage']
        delta=records[label]-reference
        comparisons[label]={'max_individual_voltage_difference':float(abs(delta).max()),'max_total_voltage_difference':float(abs(delta.sum(-1)).max())}
    comparisons['finest_step_pair']={'max_total_voltage_difference':float(abs((records['dt00625']-records['dt0125']).sum(-1)).max())}
    with np.load(out/'reference_windows.npz') as zr,np.load(out/'reference_windows_RK4.npz') as zf:
        pr=from_state(zr['states'],zr['sign']);pf=from_state(zf['states'],zf['sign'])
        errors=[]
        for w in [200,400,800]:
            diff=voltage(pr[0],pr[w//5],w)-voltage(pf[0],pf[w//5],w)
            errors.extend([float(abs(diff).max()),float(abs(diff.sum(-1)).max())])
        comparisons['window_RK4_reference_max_error']=max(errors)
    gates={
        'finest_fixed_reference':max(comparisons['dt00625'].values())<2e-5,
        'finer_step_pair':comparisons['finest_step_pair']['max_total_voltage_difference']<2e-5,
        'reference_tolerance_change':max(comparisons['tighter'].values())<2e-6,
        'all_window_cross_solver':comparisons['window_RK4_reference_max_error']<2e-5}
    result={'comparisons':comparisons,'gates':gates,'all_passed':all(gates.values()),
            'scope':'Additional high-accuracy full sweep histories at S=100 and 200. Targeted further fixed-step/tolerance checks cover alpha=.5 downward S=200 only. All W=800 replays cross-check RK4 against adaptive integration.'}
    (ROOT/'results/refinement.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['all_passed']:raise RuntimeError('Additional refinement has not met its declared limits')


if __name__=='__main__':main()
