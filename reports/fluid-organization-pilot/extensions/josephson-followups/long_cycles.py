"""Two whole turns in mean-phase coordinates, then verified periodic continuation."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from support import ROOT,P,save,table,interval,evolve
from circuit import stable_velocity,dissipation

def cycles(state,bias,alpha,max_step=.02):
    state=np.asarray(state);p0=state[:,0]%(2*np.pi);z0=state[:,1];n=len(state)
    bias=np.broadcast_to(bias,(n,));k=1+2*np.broadcast_to(alpha,(n,))
    def fun(u,y):
        z=y[:n];den=bias+np.sin(p0+u)*np.tanh(z)
        v=stable_velocity(np.stack([p0+u,z],axis=-1),1,bias,(k-1)/2)
        return np.r_[-k*np.cos(p0+u)/den,k/den,dissipation(v,(k-1)/2)*k/den]
    sol=solve_ivp(fun,(0,4*np.pi),np.r_[z0,np.zeros(2*n)],method='DOP853',rtol=1e-11,atol=1e-12,max_step=max_step,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    one=sol.sol(2*np.pi);two=sol.sol(4*np.pi)
    return sol,one[n:2*n],np.maximum(abs(one[:n]-z0),abs(two[:n]-z0)),abs(two[n:2*n]-2*one[n:2*n])

def continued_voltage(sol,period,windows):
    n=len(period);out=[]
    for w in windows:
        count=np.floor(w/period).astype(int);rem=w-count*period
        phase=np.array([brentq(lambda u:sol.sol(u)[n+j]-rem[j],0,2*np.pi,xtol=1e-12) if rem[j]>1e-12 else 0 for j in range(n)])
        out.append(2*(2*np.pi*count+phase)/w)
    return np.array(out)

def main():
    c=P['long_cycles'];z=np.load(ROOT/'data/checkpoints.npz');shape=z['state'].shape[:-1]
    state=z['state'].reshape(-1,2);bias=np.broadcast_to(z['bias'][:,None,None,None],shape).ravel();alpha=np.broadcast_to(z['alpha'][None,None,:,None],shape).ravel()
    sol,period,ret,perr=cycles(state,bias,alpha,c['phase_max_step'])
    values=continued_voltage(sol,period,c['windows']);whole=4*np.pi/period
    # Six fixed, pre-outcome selections cover every load and current.
    selected=[0,36+5,72+9,108+24+3,216+12+7,216+60+11]
    times=np.array(c['windows']);direct=evolve(state[selected],bias[selected],alpha[selected],float(times[-1]),np.r_[0,times],max_step=.5,rtol=2e-12,atol=2e-13)
    dv=2*(direct[1:,:,0]-direct[0,:,0])/times[:,None];derr=abs(dv-values[:,selected])
    u=np.linspace(0,2*np.pi,1001);sample=sol.sol(u)
    n=len(state);du=2*np.cos(state[:,0])*(np.tanh(sol.sol(2*np.pi)[:n])-np.tanh(state[:,1]))
    balance=(4*np.pi*bias-du-sol.sol(2*np.pi)[2*n:])/(1+4*np.pi*bias)
    np.savez_compressed(ROOT/'data/long_cycles.npz',state=state,bias=bias,alpha=alpha,period=period,return_error=ret,period_error=perr,windows=times,voltage=values,whole_voltage=whole,phase=u,cycle_z=sample[:n],cycle_time=sample[n:2*n],cycle_dissipation=sample[2*n:],cycle_balance=balance,direct_indices=selected,direct_states=direct,direct_voltage=dv)
    rows=[];vv=values.reshape((len(times),)+shape);ww=whole.reshape(shape)
    for a,av in enumerate(z['alpha']):
        for b,bv in enumerate(z['bias']):
            for j,w in enumerate(times):
                ci=interval(abs(vv[j,b,1,a]-vv[j,b,2,a]));rows.append({'alpha':float(av),'bias':float(bv),'window':int(w),'mean_abs_up_down_gap':ci['mean'],'ci_low':ci['ci95'][0],'ci_high':ci['ci95'][1]})
            ci=interval(abs(ww[b,1,a]-ww[b,2,a]));rows.append({'alpha':float(av),'bias':float(bv),'window':'whole_cycle','mean_abs_up_down_gap':ci['mean'],'ci_low':ci['ci95'][0],'ci_high':ci['ci95'][1]})
    table('long_cycles.csv',rows)
    checks={'max_log_separation_return_error':float(ret.max()),'max_second_minus_first_period':float(perr.max()),'max_direct_time_voltage_error':float(derr.max()),'max_alpha_zero_exact_voltage_error':float(abs(whole[alpha==0]-2*np.sqrt(bias[alpha==0]**2-1)).max()),'max_normalized_cycle_energy_error':float(abs(balance).max())}
    save('results/long_cycles.json',{'cases':len(state),'direct_time_cases':len(selected),'direct_duration':int(times[-1]),'continuation':'Whole-cycle continuation after checking phase/log-separation return; six independent direct-time comparisons','checks':checks,'gates':{'cycle_return':bool(ret.max()<c['cycle_return_gate']),'direct_voltage':bool(derr.max()<c['voltage_gate']),'alpha_zero_exact':bool(checks['max_alpha_zero_exact_voltage_error']<c['voltage_gate']),'cycle_energy_1e_6':bool(abs(balance).max()<1e-6)},'regimes':rows})
    print(checks,flush=True)
if __name__=='__main__':main()
