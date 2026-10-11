"""RCSJ pair restricted explicitly to its exactly synchronized invariant branch."""
import numpy as np
from scipy.integrate import solve_ivp
from support import ROOT,P,save,table,interval

def energy(state,beta):
    return np.asarray(beta)*state[...,1]**2-2*np.cos(state[...,0])
def segment(state,bias,beta,alpha,settle,window,max_step=.2):
    shape=state.shape;count=int(np.prod(shape[:-1]));beta=np.broadcast_to(beta,shape[:-1]);k=1+2*alpha
    start=np.concatenate([state,np.zeros(shape)],axis=-1)
    def fun(t,y):
        q=y.reshape(shape[:-1]+(4,));p,v=q[...,0],q[...,1]
        return np.stack([v,(bias-np.sin(p)-k*v)/beta,2*bias*v,2*k*v*v],axis=-1).ravel()
    sol=solve_ivp(fun,(0,settle+window),start.ravel(),method='DOP853',rtol=1e-10,atol=1e-12,max_step=max_step,t_eval=[0,settle,settle+window])
    if not sol.success:raise RuntimeError(sol.message)
    q=sol.y.T.reshape((3,)+shape[:-1]+(4,));work=q[-1,...,2];diss=q[-1,...,3]
    residual=(energy(q[-1,...,:2],beta)-energy(state,beta)-work+diss)/(1+abs(work)+diss)
    voltage=2*(q[-1,...,0]-q[1,...,0])/window
    return q,voltage,residual

def sweep(initial,betas,alpha,biases,settle,window,max_step):
    state=initial.copy();data=[];volts=[];checks=[]
    for current in list(biases)+list(biases[::-1]):
        q,v,r=segment(state,current,betas,alpha,settle,window,max_step);data.append(q);volts.append(v);checks.append(r);state=q[-1,...,:2]
    return np.array(data),np.array(volts),np.array(checks)

def main():
    c=P['capacitance'];rng=np.random.default_rng(c['seed']);beta=np.array(c['betas'])[:,None]
    initial=np.zeros((len(beta),c['starts'],2));initial[:,:,0]=rng.uniform(-.2,.2,c['starts'])
    arrays={};results={}
    for name,settle,step in [('base',c['settle'],c['max_step']),('half_step',c['settle'],c['max_step']/2),('double_settle',2*c['settle'],c['max_step'])]:
        q,v,r=sweep(initial,beta,c['alpha'],np.array(c['biases']),settle,c['window'],step)
        arrays.update({name+'_states':q,name+'_voltage':v,name+'_balance':r});results[name]={'settle':settle,'max_step':step,'max_normalized_balance_error':float(abs(r).max())}
    rows=[];n=len(c['biases'])
    for a,b in enumerate(c['betas']):
        for j,i in enumerate(c['biases']):
            up=arrays['base_voltage'][j,a];down=arrays['base_voltage'][2*n-1-j,a];ci=interval(abs(up-down))
            rows.append({'beta':b,'bias':i,'mean_up_voltage':float(up.mean()),'mean_down_voltage':float(down.mean()),'mean_abs_gap':ci['mean'],'ci_low':ci['ci95'][0],'ci_high':ci['ci95'][1],'step_change_max':float(abs(arrays['base_voltage'][:,a]-arrays['half_step_voltage'][:,a]).max()),'settling_change_at_bias':float(abs(arrays['base_voltage'][[j,2*n-1-j],a]-arrays['double_settle_voltage'][[j,2*n-1-j],a]).max())})
    step_error=float(abs(arrays['base_voltage']-arrays['half_step_voltage']).max())
    np.savez_compressed(ROOT/'data/capacitance.npz',initial=initial,betas=c['betas'],biases=c['biases'],**arrays)
    table('capacitance.csv',rows);save('results/capacitance.json',{'scope':c['scope'],'variants':results,'regimes':rows,'max_step_voltage_change':step_error,'max_settling_voltage_change':float(abs(arrays['base_voltage']-arrays['double_settle_voltage']).max()),'gates':{'power':max(r['max_normalized_balance_error'] for r in results.values())<c['power_gate'],'time_step_2e_5':step_error<2e-5},'normalization':'beta=2 e Ic r^2 C/hbar; synchronized equation beta*p_ddot+(1+2*alpha)*p_dot+sin(p)=i; total voltage=2*p_dot'})
    print('capacitance',results,'step',step_error,flush=True)
if __name__=='__main__':main()
