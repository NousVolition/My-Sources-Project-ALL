"""Dimensionless thermal-noise/mismatch sensitivity, not physically calibrated noise."""
import numpy as np
from support import ROOT,P,save,table,interval,evolve
from circuit import to_state

def mobility(alpha):
    return np.eye(2)-alpha/(1+2*alpha)*np.ones((2,2))
def drift(phi,bias,alpha,critical):
    b=np.asarray(bias)[...,None]-critical*np.sin(phi)
    return b@mobility(alpha)
def heun(phi,bias,alpha,critical,dt,increment):
    f=drift(phi,bias,alpha,critical)
    pred=phi+dt*f+increment
    return phi+.5*dt*(f+drift(pred,bias,alpha,critical))+increment

def main():
    c=P['noise'];rng=np.random.default_rng(c['seed']);m=mobility(c['alpha']);chol=np.linalg.cholesky(m)
    conditions=np.array(c['conditions']);gam=conditions[:,0];delta=conditions[:,1]
    init=rng.uniform(-np.pi,np.pi,(c['starts'],2));shape=(len(c['biases']),len(gam),c['starts'],2)
    fine=np.broadcast_to(init,shape).copy();coarse=fine.copy();b=np.array(c['biases'])[:,None,None]
    critical=np.stack([1-delta,1+delta],axis=-1)[None,:,None,:]
    scale=np.sqrt(2*gam)[None,:,None,None];dt=c['fine_dt'];steps=round(c['duration']/dt);every=round(c['record_stride']/dt)
    fs=[fine.copy()];cs=[coarse.copy()];sumdw=np.zeros_like(fine);thermal_work=np.zeros(shape[:-1]);drift_work=np.zeros(shape[:-1]);coarse_work=np.zeros(shape[:-1])
    for j in range(steps):
        dw=rng.normal(size=(len(c['biases']),1,c['starts'],2))*np.sqrt(dt)
        inc=scale*(dw@chol.T);old=fine.copy();f0=drift(old,b,c['alpha'],critical)
        predictor=old+dt*f0+inc;f1=drift(predictor,b,c['alpha'],critical)
        fine=old+.5*dt*(f0+f1)+inc;sumdw+=inc
        # Stratonovich potential budget, using midpoint restoring force.
        grad=critical*np.sin((old+fine)/2)
        drift_work+=(grad*(.5*dt*(f0+f1))).sum(-1);thermal_work+=(grad*inc).sum(-1)
        if (j+1)%2==0:
            oldc=coarse.copy();coarse=heun(coarse,b,c['alpha'],critical,2*dt,sumdw)
            coarse_work+=(critical*np.sin((oldc+coarse)/2)*(coarse-oldc)).sum(-1);sumdw.fill(0)
        if (j+1)%every==0:fs.append(fine.copy());cs.append(coarse.copy())
    fs=np.array(fs);cs=np.array(cs);t=np.arange(len(fs))*c['record_stride'];j=int(np.argmin(abs(t-c['settle'])))
    vf=(fs[-1]-fs[j]).sum(-1)/(c['duration']-c['settle']);vc=(cs[-1]-cs[j]).sum(-1)/(c['duration']-c['settle'])
    err=np.sqrt(np.mean((vf-vc)**2,axis=-1));rows=[]
    for bi,bias in enumerate(c['biases']):
        for k,(gamma,mismatch) in enumerate(conditions):
            ci=interval(vf[bi,k]);rows.append({'bias':bias,'Gamma':gamma,'mismatch':mismatch,'mean_total_voltage':ci['mean'],'ci_low':ci['ci95'][0],'ci_high':ci['ci95'][1],'paired_resolution_rms':float(err[bi,k])})
    initial_state,sign=to_state(init);bs=np.array(c['biases'])[:,None]
    reference=evolve(np.broadcast_to(initial_state,(len(bs),)+initial_state.shape).copy(),bs,c['alpha'],c['duration'],[c['settle'],c['duration']])
    noiseless=2*(reference[-1,...,0]-reference[0,...,0])/(c['duration']-c['settle'])
    du=(-critical*np.cos(fs[-1])).sum(-1)-(-critical*np.cos(fs[0])).sum(-1)
    balance=(du-drift_work-thermal_work)/(1+abs(du)+abs(drift_work+thermal_work))
    coarse_du=(-critical*np.cos(cs[-1])).sum(-1)-(-critical*np.cos(cs[0])).sum(-1)
    coarse_balance=(coarse_du-coarse_work)/(1+abs(coarse_du)+abs(coarse_work))
    for row in rows:
        bi=c['biases'].index(row['bias']);row['matched_identical_noiseless_mean']=float(noiseless[bi].mean())
    np.savez_compressed(ROOT/'data/noise.npz',time=t,fine_phase=fs,coarse_phase=cs,initial=init,bias=c['biases'],conditions=conditions,fine_voltage=vf,coarse_voltage=vc,noiseless_voltage=noiseless,drift_work=drift_work,thermal_work=thermal_work,potential_budget_residual=balance,coarse_potential_budget_residual=coarse_balance)
    # Independent free-diffusion covariance control: no Josephson restoring force.
    control=np.random.default_rng(43112).normal(size=(200000,2))@chol.T*np.sqrt(2*.001)
    cov=np.cov(control,rowvar=False);expected=2*.001*m;cov_error=np.max(abs(cov-expected))/np.max(abs(expected))
    table('noise.csv',rows);save('results/noise.json',{'regimes':rows,'independent_runs_per_regime':c['starts'],'number_of_conditions':len(rows),'averaging_window':[c['settle'],c['duration']],'common_random_numbers':'Same Brownian drivers across sensitivity conditions; fine/coarse increments coupled exactly','max_resolution_rms':float(err.max()),'resolution_gate':bool(err.max()<c['voltage_step_rms_gate']),'free_diffusion_covariance':cov.tolist(),'expected_covariance':expected.tolist(),'covariance_relative_error':float(cov_error),'covariance_gate_0_02':bool(cov_error<.02),'max_normalized_potential_budget_residual':float(abs(balance).max()),'budget_description':'Finite-step Stratonovich potential budget includes both drift work and thermal work; diagnostic, not an ordinary deterministic power balance.','interpretation':'Exploratory finite-window dimensionless sensitivity. White-noise voltage is measured through phase increments, not instantaneous derivatives. No deterministic pathwise energy-conservation claim.'})
    print('noise rms',err.max(),'covariance error',cov_error,flush=True)
if __name__=='__main__':main()
