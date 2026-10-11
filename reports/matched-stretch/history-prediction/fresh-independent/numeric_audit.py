"""Predetermined measurement screens; no statistical outcome selects a grid."""
from pathlib import Path
import json,argparse
import numpy as np
from fit_models import HERE,paths_for,load_models,make_dataset,make_matrices,save_json,sha
from analyze_test import score
from export_models import predict_model


def rms(x):return float(np.sqrt(np.mean(np.asarray(x)**2)))
def relative(a,b):return rms(a-b)/max(rms(b),1e-12)


def datasets(seeds,n,p,dt=None):return make_dataset(paths_for(seeds,n,p,dt=dt),p)


def pair_rows(a,b,p,cutoff,old):
    for key in ['group','anchor','particle']:np.testing.assert_array_equal(a[key],b[key])
    sa=np.array(old['models']['A_full']['continuous']['scale'])
    sh=np.array(old['models']['B_history']['continuous']['scale'])[len(sa):]
    xa=make_matrices(a);xb=make_matrices(b);rows=[]
    for seed in np.unique(a['group']):
        use=a['group']==seed;ea=a['target'][use]>=cutoff;eb=b['target'][use]>=cutoff
        row={'seed':int(seed),'target_RMS':rms(a['target'][use]-b['target'][use]),
            'target_relative_RMS':relative(a['target'][use],b['target'][use]),
            'event_Jaccard':float((ea&eb).sum()/(ea|eb).sum()) if (ea|eb).any() else 1.,
            'gradient_relative_RMS':relative(a['present'][use,3:12],b['present'][use,3:12]),
            'velocity_relative_RMS':relative(a['present'][use,:3],b['present'][use,:3]),
            'present_original_train_standardized_RMS':rms((a['present'][use]-b['present'][use])/sa),
            'history_original_train_standardized_RMS':rms((a['history'][use]-b['history'][use])/sh)}
        for name in ['A_full','B_history']:
            for kind in ['continuous','event']:
                pa=predict_model(old['models'][name][kind],xa[name][use]);pb=predict_model(old['models'][name][kind],xb[name][use])
                row[name+'_'+kind+'_relative_RMS']=relative(pa,pb)
        rows.append(row)
    return rows


def screens_for(rows,p,label,temporal=False,spatial=None):
    lim=p['numeric_screens'];out=[]
    def add(seed,metric,value,bound,direction='max'):
        out.append({'comparison':label,'seed':seed,'metric':metric,'value':value,'bound':bound,
                    'direction':direction,'pass':bool(value<=bound if direction=='max' else value>=bound)})
    for row in rows:
        seed=row['seed']
        add(seed,'target_relative_RMS',row['target_relative_RMS'],lim['half_step_target_relative_RMS_max'] if temporal else lim['target_relative_RMS_max'])
        add(seed,'event_Jaccard',row['event_Jaccard'],lim['event_Jaccard_min'],'min')
        for metric in ['gradient_relative_RMS','velocity_relative_RMS','present_original_train_standardized_RMS','history_original_train_standardized_RMS']:
            add(seed,metric,row[metric],lim[metric+'_max'])
        for name in ['A_full','B_history']:
            for kind in ['continuous','event']:
                key=name+'_'+kind+'_relative_RMS';add(seed,key,row[key],lim['frozen_old_prediction_relative_RMS_max'])
        if temporal:
            ref=next(r['target_RMS'] for r in spatial if r['seed']==seed)
            add(seed,'half_step_target_fraction_of_spatial',row['target_RMS']/max(ref,1e-30),lim['half_step_target_fraction_of_spatial_max'])
    return out


def volumes(seeds,n,p,dt=None):
    out=[]
    for path in paths_for(seeds,n,p,dt=dt):
        meta=json.loads(path.with_suffix('.json').read_text())
        value=max(d['max_det_error'] for d in meta['diagnostics'])
        out.append({'comparison':f'grid{n}','seed':meta['seed'],'metric':'probe_volume_error','value':value,
                    'bound':p['numeric_screens']['probe_volume_error_max'],'direction':'max',
                    'pass':value<=p['numeric_screens']['probe_volume_error_max']})
    return out


def model_screens(a,b,p,manifest,models,label):
    entry=manifest['cases']['primary'];selected={k:models['primary'][k] for k in ['A_full','B_history']}
    ra,pa=score(a,entry,selected);rb,pb=score(b,entry,selected);lim=p['numeric_screens'];out=[];deltas=[]
    denom=entry['models']['A_full']['validation']['macro']['RMSE']
    def add(seed,metric,value,bound):
        out.append({'comparison':label,'seed':seed,'metric':metric,'value':float(value),'bound':bound,'direction':'max','pass':bool(value<=bound)})
    for seed in np.unique(a['group']):
        key=str(int(seed));use=a['group']==seed
        for name in selected:
            add(int(seed),name+'_new_probability_relative_RMS',relative(pa[name+'_probability'][use],pb[name+'_probability'][use]),lim['new_model_probability_relative_RMS_max'])
        da=ra['B_history']['spatial'][key]['mean_AP']-ra['A_full']['spatial'][key]['mean_AP']
        db=rb['B_history']['spatial'][key]['mean_AP']-rb['A_full']['spatial'][key]['mean_AP']
        deltas.append(db-da);add(int(seed),'spatial_AP_contrast_drift',abs(db-da),lim['new_model_spatial_AP_contrast_per_run_drift_max'])
        ca=ra['B_history']['per_run'][key]['RMSE']-ra['A_full']['per_run'][key]['RMSE']
        cb=rb['B_history']['per_run'][key]['RMSE']-rb['A_full']['per_run'][key]['RMSE']
        add(int(seed),'RMSE_contrast_drift_relative_validation_A',abs(cb-ca)/max(denom,1e-30),lim['new_model_RMSE_contrast_relative_validation_A_max'])
    add('macro','spatial_AP_contrast_drift',abs(np.mean(deltas)),lim['new_model_spatial_AP_contrast_macro_drift_max'])
    return out


def audit(n,mode):
    p=json.loads((HERE/'protocol.json').read_text());old=json.loads((HERE.parent/'results'/'primary_models.json').read_text())
    finer={80:112,112:160}[n]
    seeds=p['test_numerical_controls']['seeds'] if mode=='test' else p['pilot']['seeds']
    a=datasets(seeds,n,p);b=datasets(seeds,finer,p);rows=pair_rows(a,b,p,p['pilot']['event_cutoff'],old)
    screens=screens_for(rows,p,f'{n}→{finer}')+volumes(seeds,n,p)+volumes(seeds,finer,p)
    halfseeds=seeds if mode=='test' else p['pilot']['half_step_seeds']
    ta=datasets(halfseeds,n,p);tb=datasets(halfseeds,n,p,p['pilot']['half_step_dt'])
    temporal=pair_rows(ta,tb,p,p['pilot']['event_cutoff'],old)
    screens+=screens_for(temporal,p,f'{n} half-step',temporal=True,spatial=rows)
    screens+=volumes(halfseeds,n,p,p['pilot']['half_step_dt'])
    if mode=='pilot':
        previous={80:56,112:80}[n];prev=datasets(seeds,previous,p)
        before=pair_rows(prev,a,p,p['pilot']['event_cutoff'],old)
        screens+=screens_for(before,p,f'{previous}→{n}')
        for row,oldrow in zip(rows,before):
            value=row['target_RMS']/max(oldrow['target_RMS'],1e-30)
            screens.append({'comparison':'adjacent contraction','seed':row['seed'],'metric':'target_RMS_ratio',
                            'value':value,'bound':1.,'direction':'max','pass':value<1.})
    if mode in ['model','test']:
        manifest_path=HERE/(json.loads((HERE/'test-lock.json').read_text())['manifest'] if mode=='test' else f'models-{n}/manifest.json')
        m,models=load_models(manifest_path);cutoff=m['cases']['primary']['train_event_cutoff']
        rows=pair_rows(a,b,p,cutoff,old);temporal=pair_rows(ta,tb,p,cutoff,old)
        # Event sets now use the fresh train-only cutoff as well.
        screens+=screens_for(rows,p,'fresh cutoff spatial')+screens_for(temporal,p,'fresh cutoff temporal',True,rows)
        screens+=model_screens(a,b,p,m,models,'fresh fitted spatial')+model_screens(ta,tb,p,m,models,'fresh fitted temporal')
    result={'mode':mode,'grid':n,'finer_grid':finer,'seeds':seeds,'rows':rows,'temporal':temporal,'screens':screens,
        'protocol_sha256':sha(HERE/'protocol.json'),'old_models_sha256':sha(HERE.parent/'results'/'primary_models.json'),
        'all_screens_pass':all(s['pass'] for s in screens)}
    save_json(HERE/('test-numerical-audit.json' if mode=='test' else f'{mode}-audit-{n}.json'),result)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--grid',type=int,choices=[80,112],required=True)
    ap.add_argument('--mode',choices=['pilot','model','test'],required=True);args=ap.parse_args()
    r=audit(args.grid,args.mode);print(f"{sum(s['pass'] for s in r['screens'])}/{len(r['screens'])} screens pass")
