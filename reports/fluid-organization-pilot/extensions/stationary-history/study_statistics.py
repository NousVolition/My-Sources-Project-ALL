"""Run-cluster uncertainty. Snapshots, strengths and history donors are dependent."""
import numpy as np

def contrast(a,b,reps=10000,seed=501):
 a=np.asarray(a);b=np.asarray(b);rng=np.random.default_rng(seed);idx=rng.integers(0,len(a),(reps,len(a)))
 diff=(a-b)[idx].mean(1);relative=100*(1-a[idx].mean(1)/b[idx].mean(1))
 return dict(mean_error_difference=float(np.mean(a-b)),difference_interval=np.quantile(diff,[.025,.975]).tolist(),error_reduction_percent=float(100*(1-a.mean()/b.mean())),reduction_interval=np.quantile(relative,[.025,.975]).tolist(),independent_runs=len(a))

def donor_contrast(a,matrix,reps=10000,seed=502):
 n=len(a);mask=~np.eye(n,dtype=bool);base=float(matrix[mask].mean());rng=np.random.default_rng(seed);diff=[];relative=[]
 for start in range(0,reps,250):
  idx=rng.integers(0,n,(min(250,reps-start),n));valid=idx[:,:,None]!=idx[:,None,:];denom=valid.sum((1,2));assert np.all(denom>0)
  b=(matrix[idx[:,:,None],idx[:,None,:]]*valid).sum((1,2))/denom;v=a[idx].mean(1);diff.extend(v-b);relative.extend(100*(1-v/b))
 return dict(mean_error_difference=float(np.mean(a)-base),difference_interval=np.quantile(diff,[.025,.975]).tolist(),error_reduction_percent=float(100*(1-np.mean(a)/base)),reduction_interval=np.quantile(relative,[.025,.975]).tolist(),independent_runs=n,null_mean_run_rmse=base)

def stationarity(records,p):
 windows=p['stationarity_windows'];keys=['energy','dissipation','integral_length'];means=[];late=[]
 for record in records:
  rows=record['baseline'];pair=[]
  for a,b in windows:
   selected=[r for r in rows if a<=r['t']<b];assert selected
   pair.append([np.mean([r[k] for r in selected]) for k in keys])
  means.append(pair);late.append([r for r in rows if windows[0][0]<=r['t']<windows[-1][-1]])
 means=np.array(means);n=len(means);rng=np.random.default_rng(702);idx=rng.integers(0,n,(10000,n));out={}
 for j,k in enumerate(keys):
  x=means[:,0,j];y=means[:,1,j];relative=2*(y[idx].mean(1)-x[idx].mean(1))/(y[idx].mean(1)+x[idx].mean(1));ci=np.quantile(relative,[.025,.975]);out[k]=dict(early_mean=float(x.mean()),late_mean=float(y.mean()),relative_drift=float(2*(y.mean()-x.mean())/(y.mean()+x.mean())),drift_interval=ci.tolist(),equivalence_pass=bool(np.max(abs(ci))<=p['gates']['stationarity_relative_drift']))
 diss=np.array([np.mean([r['dissipation'] for r in rows]) for rows in late]);ci=np.quantile(diss[idx].mean(1)/p['power']-1,[.025,.975]);out['power_balance']=dict(relative_mismatch=float(diss.mean()/p['power']-1),interval=ci.tolist(),pass_check=bool(np.max(abs(ci))<=p['gates']['power_relative_mismatch']))
 relam=np.array([np.mean([r['re_lambda'] for r in rows]) for rows in late]);flux=np.array([np.mean([r['flux_k4'] for r in rows]) for rows in late]);comp=np.array([np.mean([r['component_energies'] for r in rows],axis=0) for rows in late]);turn=np.array([np.trapezoid([1/r['turnover'] for r in rec['baseline'] if r['t']<=p['events'][0]],[r['t'] for r in rec['baseline'] if r['t']<=p['events'][0]]) for rec in records])
 out.update(independent_runs=n,re_lambda_mean=float(relam.mean()),re_lambda_interval=np.quantile(relam[idx].mean(1),[.025,.975]).tolist(),flux_mean=float(flux.mean()),flux_interval=np.quantile(flux[idx].mean(1),[.025,.975]).tolist(),component_ratio_to_isotropic=(comp.mean(0)/comp.mean()).tolist(),minimum_spin_turnovers=float(turn.min()))
 spectra=np.array([r['window_spectra'] for r in records]);distance=[]
 for offset in range(0,len(idx),250):
  means_sp=spectra[idx[offset:offset+250]].mean(1);distance.extend(.5*np.abs(means_sp[:,1]-means_sp[:,0]).sum(1))
 bound=float(np.quantile(distance,.95));out['spectrum']=dict(total_variation=float(.5*np.abs(spectra[:,1].mean(0)-spectra[:,0].mean(0)).sum()),upper_95_percent_bound=bound,pass_check=bound<=p['gates']['spectral_total_variation_max'])
 out['passed']=bool(all(out[k]['equivalence_pass'] for k in keys) and out['power_balance']['pass_check'] and out['spectrum']['pass_check'] and relam.mean()>p['gates']['mean_re_lambda_min'] and out['flux_interval'][0]>p['gates']['flux_fraction_min']*p['power'] and turn.min()>p['gates']['spin_turnovers_min'])
 return out
