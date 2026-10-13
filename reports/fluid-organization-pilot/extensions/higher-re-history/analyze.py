"""Descriptive physics, matched refinement forecasts and scientific figures."""
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from common import ROOT,read,save,protocol
from features import inputs,target
from models import predict,run_errors
from forecast import verify_freeze
from study_statistics import contrast

def table(name,rows):
 with (ROOT/'results'/name).open('w',newline='',encoding='utf-8') as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def band(v,reps=3000,seed=220):
 rng=np.random.default_rng(seed);idx=rng.integers(0,len(v),(reps,len(v)));means=v[idx].mean(1)
 return v.mean(0),np.quantile(means,[.025,.975],axis=0)
def analyze():
 p=protocol();m=verify_freeze();r=read('results/prediction_summary.json');checks=read('results/numerical_checks.json');folder=ROOT/'results/figures';folder.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':135,'savefig.dpi':170})
 metadata=[];data=[];runrows=[]
 for seed in p['splits']['test']:
  q=read(f'data/runs/base_{seed}.json');metadata.append(q)
  with np.load(ROOT/'data/runs'/f'base_{seed}.npz') as z:data.append({k:z[k] for k in z.files})
  late=[row for row in q['baseline'] if 20<=row['t']<44]
  runrows.append(dict(seed=seed,**{key:float(np.mean([x[key] for x in late])) for key in ['energy','dissipation','re_lambda','integral_re','integral_length','turnover','flux_k4','tail_energy']},minimum_baseline_keta=min(x['keta'] for x in q['baseline']),maximum_baseline_cfl=max(x['cfl'] for x in q['baseline'])))
 table('flow_regimes.csv',runrows)
 physics=[]
 curves=np.array([z['energy_ratio'].mean(0) for z in data]);budgets=np.array([z['difference_budget'].mean(0) for z in data]);tt=data[0]['response_time'][0]
 for k,alpha in enumerate(p['strengths']):
  final=curves[:,-1,k];area=np.array([target(z).reshape(len(p['events']),-1)[:,k].mean() for z in data]);pdq=budgets[:,-1,k,:]
  _,final_ci=band(final);_,area_ci=band(area)
  physics.append(dict(relative_impulse=alpha,independent_runs=len(data),events_per_run=len(p['events']),mean_final_energy_ratio=float(final.mean()),final_ratio_ci_low=float(final_ci[0]),final_ratio_ci_high=float(final_ci[1]),mean_integrated_ratio=float(area.mean()),integrated_ratio_ci_low=float(area_ci[0]),integrated_ratio_ci_high=float(area_ci[1]),mean_production_integral=float(pdq[:,0].mean()),mean_dissipation_integral=float(pdq[:,1].mean()),mean_forcing_difference_integral=float(pdq[:,2].mean()),fraction_of_events_with_gain_above_1_01=float(np.mean([np.max(z['energy_ratio'][:,:,k],axis=1)>1.01 for z in data]))))
 table('response_regimes.csv',physics)
 secondary=[]
 with np.load(ROOT/'data/predictions.npz') as z:
  y=z['truth'].reshape(len(data),len(p['events']),len(p['strengths']));a=z['history'].reshape(y.shape);b=z['present'].reshape(y.shape)
  for k,alpha in enumerate(p['strengths']):
   ea=np.sqrt(((a[:,:,k]-y[:,:,k])**2).mean(1));eb=np.sqrt(((b[:,:,k]-y[:,:,k])**2).mean(1));q=contrast(ea,eb,p['bootstrap_replicates'],1900+k)
   secondary.append(dict(strength=alpha,history_mean_run_rmse=float(ea.mean()),present_mean_run_rmse=float(eb.mean()),reduction_percent=q['error_reduction_percent'],ci_low=q['reduction_interval'][0],ci_high=q['reduction_interval'][1],interpretation='Descriptive subgroup comparison; not the joint confirmatory support criterion'))
 table('prediction_by_strength.csv',secondary)
 # Refinement uses frozen models and the same paired events, never refits.
 refinement=[]
 for seed in p['refinement_seeds']:
  for label in ['base','half','fine']:
   with np.load(ROOT/'data/runs'/f'{label}_{seed}.npz') as z:
    selected={k:z[k] for k in z.files}
    if label=='base':
     ids=[int(np.flatnonzero(z['events']==t)[0]) for t in p['refinement_events']]
     for k in ['history','events','initial_rates','response_time','energy_ratio']:selected[k]=selected[k][ids]
    c,a,h=inputs(selected);y=target(selected)
   row=dict(seed=seed,label=label)
   for name,x in [('global_present',c),('present',a),('history',np.c_[a,h])]:row[name+'_rmse']=float(np.sqrt(np.mean((predict(m['models'][name],x)-y)**2)))
   row['history_minus_present']=row['history_rmse']-row['present_rmse'];refinement.append(row)
 table('refinement_predictions.csv',refinement)
 drift=[]
 for seed in p['refinement_seeds']:
  a={x['label']:x for x in refinement if x['seed']==seed}
  drift.append(dict(seed=seed,step_contrast_change=a['half']['history_minus_present']-a['base']['history_minus_present'],grid_contrast_change=a['fine']['history_minus_present']-a['half']['history_minus_present']))
 save('results/refinement_prediction_summary.json',dict(runs=drift,max_absolute_step_contrast_change=max(abs(x['step_contrast_change']) for x in drift),max_absolute_grid_contrast_change=max(abs(x['grid_contrast_change']) for x in drift),scope='Four dependent matched restart controls; not population-wide continuum extrapolation'))
 save('results/physical_summary.json',dict(independent_runs=len(data),mean_statistics={key:float(np.mean([r[key] for r in runrows])) for key in ['energy','dissipation','re_lambda','integral_re','integral_length','turnover','flux_k4','tail_energy']},minimum_keta=min([r['minimum_baseline_keta'] for r in runrows]+[float(z['response_checks'][...,0].min()) for z in data]),maximum_response_cfl=max(float(z['response_cfl'].max()) for z in data),maximum_difference_budget_residual=max(float(z['response_checks'][...,3].max()) for z in data),responses=physics))
 # Statistics are averaged within a run before bootstrap; frames are not independent.
 times=data[0]['baseline_time'];e=np.array([[q['energy'] for q in m['baseline']] for m in metadata]);d=np.array([[q['dissipation'] for q in m['baseline']] for m in metadata])
 fig,ax=plt.subplots(2,2,figsize=(11,7.4),layout='constrained')
 for values,a,title in [(e,ax[0,0],'Kinetic energy'),(d,ax[0,1],'Dissipation')]:
  avg,ci=band(values);a.plot(times,avg,color='#176b87');a.fill_between(times,*ci,color='#176b87',alpha=.18);a.axvspan(20,32,color='#d8e3eb',alpha=.25);a.axvspan(32,44,color='#f3ddbf',alpha=.25);a.set(xlabel='Baseline time',ylabel=title,title=title+' across independent runs')
 ax[0,1].axhline(p['power'],color='#c46b22',ls='--',label='Input power');ax[0,1].legend()
 relam=np.array([[q['re_lambda'] for q in m['baseline']] for m in metadata]);avg,ci=band(relam);ax[1,0].plot(times,avg,color='#176b87');ax[1,0].fill_between(times,*ci,color='#176b87',alpha=.18);ax[1,0].axhline(p['upstream']['previous_measured_re_lambda'],color='#777',ls=':',label='Earlier stationary mean');ax[1,0].set(xlabel='Baseline time',ylabel='Approximate Taylor Reynolds number',title='Higher Reynolds number with sustained forcing');ax[1,0].legend()
 spectra=np.array([z['spectra'][(times>=20)&(times<44)].mean(0) for z in data]);avg,ci=band(spectra);k=np.arange(len(avg))*2*np.pi/6;mask=(k>0)&(avg>1e-18);ax[1,1].loglog(k[mask],avg[mask],color='#176b87');ax[1,1].fill_between(k[mask],ci[0,mask],ci[1,mask],alpha=.18,color='#176b87');ax[1,1].set(xlabel='Physical wavenumber',ylabel='Energy per shell',title='Spectrum during the assessment windows')
 fig.suptitle('Stationarity diagnostics · shading is a 95% bootstrap interval across runs',fontsize=12);fig.savefig(folder/'stationarity.png');plt.close(fig)
 labels=['global_present','present','history','training_prior','event_prior'];names=['Global\ncurrent','Richer\ncurrent','Matched\nhistory','Strength\nprior','Time and\nstrength prior'];values=[r['mean_run_rmse'][k] for k in labels];values.append(r['contrasts']['wrong_history']['null_mean_run_rmse']);names.append('Wrong\nhistory')
 fig,ax=plt.subplots(1,2,figsize=(11.3,4.8),layout='constrained');ax[0].bar(names,values,color=['#8196a3','#477689','#197a66','#b6b6b6','#999999','#cf964f']);ax[0].set(ylabel='Mean run RMSE',title='Independent test prediction error');ax[0].tick_params(axis='x',labelsize=8)
 order=['present','global_present','wrong_history','training_prior','event_prior'];textnames=['Richer current','Global current','Wrong history','Strength prior','Time and strength prior']
 for j,key in enumerate(order):
  q=r['contrasts'][key];v=q['error_reduction_percent'];lo,hi=q['reduction_interval'];ax[1].errorbar(v,j,xerr=[[v-lo],[hi-v]],fmt='o',color='#197a66',capsize=4)
 ax[1].axvline(0,color='#777',ls='--');ax[1].set(yticks=range(len(order)),yticklabels=textnames,xlabel='Error reduction from matched history (%)',title='Paired 95% intervals across independent runs');ax[1].invert_yaxis()
 fig.suptitle('History benefit '+('meets the statistical comparison criterion' if r['history_supported'] else 'does not meet the statistical comparison criterion'),fontsize=12);fig.savefig(folder/'prediction.png');plt.close(fig)
 fig,ax=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
 for j,alpha in enumerate(p['strengths']):
  avg,ci=band(curves[:,:,j]);ax[0].plot(tt,avg,label=f'Impulse {alpha:g} × RMS speed');ax[0].fill_between(tt,*ci,alpha=.15)
  for col,name,style in [(0,'Nonlinear production','-'),(1,'Viscous loss','--'),(2,'Forcing difference',':')]:ax[1].plot(tt,budgets[:,:,j,col].mean(0),ls=style,label=f'{alpha:g}: {name}')
 ax[0].axhline(1,color='#888',ls=':');ax[0].set(xlabel='Time after impulse',ylabel='Difference energy / initial difference energy',title='Disturbance growth or decay');ax[0].legend(fontsize=8)
 ax[1].set(xlabel='Time after impulse',ylabel='Integrated budget term / initial difference energy',title='Production, dissipation and forcing');ax[1].legend(fontsize=7)
 fig.savefig(folder/'response_budgets.png');plt.close(fig)
 fig,ax=plt.subplots(1,2,figsize=(10.5,4.6),layout='constrained');pairs=checks['refinements']
 for label,color in [('half','#176b87'),('fine','#c46b22')]:
  vals=[q['max_curve_relative_change'] for q in pairs if q['label']==label];ax[0].plot(range(len(vals)),vals,'o-',color=color,label='Half step vs base' if label=='half' else 'Finer grid vs half step')
 ax[0].set(yscale='log',xticks=range(len(p['refinement_seeds'])),xticklabels=[str(x) for x in p['refinement_seeds']],ylabel='Maximum normalized curve change',xlabel='Matched restart run',title='Time and grid checks');ax[0].legend(fontsize=8)
 for label,color in [('base','#777'),('half','#176b87'),('fine','#c46b22')]:ax[1].plot(range(len(p['refinement_seeds'])),[q['history_minus_present'] for q in refinement if q['label']==label],'o-',label=label,color=color)
 ax[1].axhline(0,color='#aaa',ls=':');ax[1].set(xticks=range(len(p['refinement_seeds'])),xticklabels=[str(x) for x in p['refinement_seeds']],xlabel='Matched restart run',ylabel='History minus richer-current RMSE',title='Frozen predictions under refinement');ax[1].legend(fontsize=8)
 fig.savefig(folder/'resolution.png');plt.close(fig)
 print('Tables and four figures generated.')
if __name__=='__main__':
 with threadpool_limits(limits=1):analyze()
