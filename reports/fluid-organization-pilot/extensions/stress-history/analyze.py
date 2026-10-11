"""Summarize completed independent test fields; descriptive intervals resample fields."""
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,save,table,read_protocol
from features import inputs,target
from models import predict
from forecast import verify_freeze
def interval(v,seed=17):
    a=np.asarray(v);rng=np.random.default_rng(seed);means=a[rng.integers(0,len(a),(4000,len(a)))].mean(1)
    return {'mean':float(a.mean()),'ci95':np.quantile(means,[.025,.975]).tolist(),'min':float(a.min()),'max':float(a.max())}
def main():
    p=read_protocol();m=verify_freeze();summary=json.loads((ROOT/'results/prediction_summary.json').read_text(encoding='utf-8'));figs=ROOT/'results/figures';figs.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})
    arrays={};regimes=[];responses=[];controls=[];descriptive={};strength_rows=[]
    with (ROOT/'results/predictions.csv').open(encoding='utf-8') as f:prediction_rows=list(csv.DictReader(f))
    for i,regime in enumerate(p['regimes']):
        for j,strength in enumerate(p['kick_rms_over_baseline_rms']):
            rows=[r for r in prediction_rows if r['regime']==regime and float(r['relative_strength'])==strength]
            error={name:np.array([float(r[name])-float(r['observed']) for r in rows]) for name in ['present','history']};rms={name:float(np.sqrt(np.mean(v*v))) for name,v in error.items()}
            rng=np.random.default_rng(3500+10*i+j);idx=rng.integers(0,len(rows),(4000,len(rows)));a=np.sqrt(np.mean(error['history'][idx]**2,axis=1));b=np.sqrt(np.mean(error['present'][idx]**2,axis=1));lo,hi=np.quantile(100*(1-a/b),[.025,.975])
            strength_rows.append({'regime':regime,'relative_strength':strength,'independent_fields':len(rows),'strong_present_rmse':rms['present'],'history_rmse':rms['history'],'error_reduction_percent':100*(1-rms['history']/rms['present']),'descriptive_reduction_ci95_low':float(lo),'descriptive_reduction_ci95_high':float(hi)})
    table('prediction_by_strength.csv',strength_rows)
    save('results/strength_prediction_scope.json',{'scope':'Descriptive strength breakdown; pooled per-strength RMSE across fields, distinct from the primary mean per-field RMSE across both strengths. Unadjusted paired 95% intervals, conditional on frozen models. Four comparisons are not additional confirmatory tests.'})
    for regime,c in p['regimes'].items():
        ee=[];rr=[];ss=[];tr=[];mm=[];diag=[];spectra=[];grad=[]
        for seed in p['splits']['test']:
            path=ROOT/'data/runs'/f'{regime}_base_{seed}.npz';meta=json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
            with np.load(path) as z:
                time=z['time'];e=z['perturbation_energy']/z['initial_perturbation_energy'];ee.append(e);rr.append(z['radius']);ss.append(z['participation']);tr.append(z['transfer']/z['initial_perturbation_energy'][None,:,None]);
                edges=z['positions'][:,:,:,1:]-z['positions'][:,:,:,0:1];logs=np.log(np.maximum(np.linalg.norm(edges,axis=-1)/p['marker_edge'],1e-12));mm.append(logs.mean((2,3)))
            dd=meta['diagnostics'];diag.append(dd[0][0]);spectra.append([q[0]['spectrum'] for q in dd]);grad.append([q[0]['vortex_stretching'] for q in dd])
        a={'energy':np.array(ee),'radius':np.array(rr),'participation':np.array(ss),'transfer':np.array(tr),'marker':np.array(mm),'spectra':np.array(spectra),'stretch':np.array(grad)};arrays[regime]=a
        quality={k:interval([d[k] for d in diag]) for k in ['Taylor_Re','integral_Re','fraction_k_ge_4','flux_k4','nonlinear_over_viscous_rms','vortex_stretching']};descriptive[regime]={'intervention_time_baseline':quality,'response':{}}
        regimes.append({'regime':regime,**c,'train_fields':len(p['splits']['train']),'validation_fields':len(p['splits']['validation']),'heldout_fields':len(p['splits']['test']),**{k:quality[k]['mean'] for k in quality}})
        for j,strength in enumerate(p['kick_rms_over_baseline_rms']):
            vals={'energy_area':np.trapezoid(a['energy'][:,:,j],time,axis=1),'final_energy_ratio':a['energy'][:,-1,j],'peak_energy_ratio':a['energy'][:,:,j].max(1),'radius_change':a['radius'][:,-1,j]-a['radius'][:,0,j],'participation_final':a['participation'][:,-1,j],'integrated_production_over_initial_energy':a['transfer'][:,-1,j,0],'integrated_dissipation_over_initial_energy':a['transfer'][:,-1,j,1],'marker_log_length_excess':a['marker'][:,-1,j+1]-a['marker'][:,-1,0]}
            quantities={k:interval(v) for k,v in vals.items()};descriptive[regime]['response'][str(strength)]=quantities
            responses.append({'regime':regime,'relative_strength':strength,'fields_with_sampled_energy_growth_above_1_01':int((vals['peak_energy_ratio']>1.01).sum()),**{k:v['mean'] for k,v in quantities.items()}})
        for seed in p['refinement_seeds']:
            for label in ['base','half','fine']:
                with np.load(ROOT/'data/runs'/f'{regime}_{label}_{seed}.npz') as z:
                    _,a0,h0=inputs(z,p['history_lags']);y=target(z);pr=predict(m['models'][regime]['present'],a0);hr=predict(m['models'][regime]['history'],np.c_[a0,h0]);errors=[float(np.sqrt(np.mean((v-y)**2))) for v in [pr,hr]]
                controls.append({'regime':regime,'seed':seed,'settings':label,'present_rmse':errors[0],'history_rmse':errors[1],'history_minus_present':errors[1]-errors[0],'target_weak':float(y[0]),'target_strong':float(y[1])})
    table('parameter_regimes.csv',regimes);table('response_regimes.csv',responses);table('prediction_refinements.csv',controls);save('results/flow_summary.json',descriptive)
    refinements={}
    for regime in p['regimes']:
        group=[q for q in controls if q['regime']==regime];means={label:float(np.mean([q['history_minus_present'] for q in group if q['settings']==label])) for label in ['base','half','fine']};base={q['seed']:q for q in group if q['settings']=='base'}
        drift=max(abs(q['history_minus_present']-base[q['seed']]['history_minus_present']) for q in group if q['settings']=='fine')
        refinements[regime]={'mean_history_minus_present_on_four_fields':means,'maximum_grid_contrast_drift':float(drift),'primary_mean_contrast':summary[regime]['contrasts']['history_vs_present']['mean_rmse_difference'],'drift_exceeds_primary_contrast_magnitude':drift>abs(summary[regime]['contrasts']['history_vs_present']['mean_rmse_difference']),'coarse_fine_mean_sign_changed':means['base']*means['fine']<0,'scope':'Descriptive dependent four-field checks; not an independent significance test or ensemble continuum proof'}
    save('results/refinement_forecast_summary.json',refinements)
    fig,axs=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    colors=['#287eaf','#d35d24']
    for row,(regime,a) in enumerate(arrays.items()):
        for j,strength in enumerate(p['kick_rms_over_baseline_rms']):
            for col,key in enumerate(['energy','radius','participation']):
                v=a[key][:,:,j];mean=v.mean(0);se=v.std(0,ddof=1)/np.sqrt(len(v));axs[row,col].plot(time,mean,color=colors[j],label=f'Impulse RMS / baseline RMS = {strength}');axs[row,col].fill_between(time,mean-1.96*se,mean+1.96*se,color=colors[j],alpha=.14)
        axs[row,0].axhline(1,color='.5',lw=1,ls=':');axs[row,0].set_ylabel(regime+' flow');axs[row,0].set_title('Perturbation energy / initial value');axs[row,1].set_title('Energy-weighted periodic radius');axs[row,2].set_title('Spatial participation fraction')
        for ax in axs[row]:ax.set_xlabel('Time after impulse')
    axs[0,0].legend(fontsize=8);fig.suptitle('Independent held-out flows: mean response and pointwise 95% intervals');fig.savefig(figs/'disturbance.png');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for regime,a in arrays.items():
        s=a['spectra'][:,0].mean(0);positive=s>0;axs[0].semilogy(np.arange(1,len(s)+1)[positive],s[positive],label=regime)
        axs[1].plot(time,a['stretch'].mean(0),label=regime)
        for j,strength in enumerate(p['kick_rms_over_baseline_rms']):axs[2].plot(time,(a['marker'][:,:,j+1]-a['marker'][:,:,0]).mean(0),label=f'{regime}, {strength}')
    axs[0].set(xlabel='Fourier shell |k| in box modes',ylabel='Kinetic energy per shell',title='Spectrum before impulse');axs[0].legend();axs[1].set(xlabel='Time after impulse',ylabel='Enstrophy production',title='Baseline vortex stretching');axs[1].legend();axs[2].set(xlabel='Time after impulse',ylabel='Excess mean log neighbor length',title='Passive-marker deformation');axs[2].legend(fontsize=8);fig.savefig(figs/'flow_scales.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for regime,a in arrays.items():
        for j,strength in enumerate(p['kick_rms_over_baseline_rms']):
            for col,label in enumerate(['production','dissipation']):axs[col].plot(time,a['transfer'][:,:,j,col].mean(0),label=f'{regime}, {strength}')
    for ax,title in zip(axs,['Energy transferred from baseline','Viscous perturbation-energy loss']):ax.set(xlabel='Time after impulse',ylabel='Integral / initial perturbation energy',title=title);ax.legend(fontsize=8)
    fig.savefig(figs/'energy_budget.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
    for ax,(regime,s) in zip(axs,summary.items()):
        means=s['mean_field_rmse'];values=[means['conventional'],means['present'],means['history'],s['contrasts']['history_vs_shuffled_mean']['null_mean_field_rmse'],means['prior_by_strength']]
        ax.bar(np.arange(5),values,color=['#aaa','#287eaf','#d35d24','#9a9a9a','#54634c']);ax.set_xticks(np.arange(5),['Global present','Richer present','With history','Expected shuffled','Training prior'],rotation=30,ha='right');ax.set(ylabel='Mean independent-field RMSE (model time)',title=regime,ylim=(0,1.6*max(values)))
        q=s['contrasts']['history_vs_present'];null=s['contrasts']['history_vs_shuffled_mean']
        ax.text(.02,.96,f'History vs richer present: {q["relative_improvement_percent"]:+.1f}%\n97.5% interval [{q["relative_improvement_interval"][0]:+.1f}, {q["relative_improvement_interval"][1]:+.1f}]%\nHistory vs expected shuffled: {null["relative_improvement_percent"]:+.1f}%\n97.5% interval [{null["relative_improvement_interval"][0]:+.1f}, {null["relative_improvement_interval"][1]:+.1f}]%',transform=ax.transAxes,va='top',fontsize=9)
    fig.suptitle('History does not clearly beat the null controls in either regime')
    fig.savefig(figs/'prediction.png');plt.close(fig)
    save('results/descriptive_scope.json',{'response_intervals':'Descriptive pointwise 95% intervals; no additional confirmatory claims.','primary_intervals':'97.5% per regime; Bonferroni protection for two history-versus-present comparisons.','refinement_scope':'Four prechosen held-out initial fields per regime. Dependent numerical controls, not independent validation samples.'})
if __name__=='__main__':main()
