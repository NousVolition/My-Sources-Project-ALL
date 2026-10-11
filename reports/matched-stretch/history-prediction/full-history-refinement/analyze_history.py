"""Score frozen coefficients; compare complete histories and paired metric drift."""
from pathlib import Path
import sys
import json
import csv
import numpy as np
HERE=Path(__file__).resolve().parent; BASE=HERE.parent
sys.path.insert(0,str(HERE))
from run_history import sha, save_json, sources, source_hash_matches
from features import dataset, shuffle_blocks, PRESENT_NAMES, HISTORY_NAMES
from export_models import predict_model
from evaluate import summarize, paired


def checked(path,config):
    """Read-only analysis permits newline-equivalent source, never changed numerical text."""
    m=json.loads(path.with_suffix('.json').read_text())
    assert m['config']==config and sha(path)==m['sha256']
    for rel,digest in m['sources'].items():
        source=BASE.parent/'numerics.py' if rel=='numerics.py' else BASE/rel
        assert source_hash_matches(source.read_bytes(),digest),rel
    if 'study_sources' in m:
        assert set(m['study_sources'])==set(sources())
        for rel,digest in m['study_sources'].items():
            assert source_hash_matches((BASE.parent/rel).read_bytes(),digest),rel
    assert m['max_cfl']<.5
    assert max(abs(d['energy_budget_residual']) for d in m['diagnostics'])<.005
    assert max(d['divergence_max'] for d in m['diagnostics'])<1e-10
    with np.load(path,allow_pickle=False) as z:
        assert all(np.isfinite(z[k]).all() for k in z.files)
        np.testing.assert_allclose(z['time'],np.arange(round(config['end']/config['save_dt'])+1)*config['save_dt'],atol=3e-16,rtol=0)
        expected=np.random.default_rng(m['seed']+80000).uniform(-3,3,(config['markers'],3))
        np.testing.assert_array_equal(z['positions'][0],expected)
    return m


def score(d,models,p):
    event=(d['target']>=models['event_target_cutoff']).astype(int)
    out={}; predictions={}
    for name,X in [('A_full',d['present']),('B_history',np.column_stack([d['present'],d['history']]))]:
        m=models['models'][name]
        y=predict_model(m['continuous'],X); prob=predict_model(m['event'],X)
        out[name]=summarize(d['target'],event,y,prob,d['group'],m['decision_threshold'])
        predictions[name+'_prediction']=y; predictions[name+'_probability']=prob
    out['paired']=paired(out['A_full'],out['B_history'],p)
    return out,predictions


def compare_datasets(a,b,pa,pb,models,p):
    for key in ['group','anchor','particle']: np.testing.assert_array_equal(a[key],b[key])
    out=[]; cutoff=models['event_target_cutoff']
    sa=np.array(models['models']['A_full']['continuous']['scale'])
    sh=np.array(models['models']['B_history']['continuous']['scale'])[len(PRESENT_NAMES):]
    for seed in p['test_seeds']:
        ix=b['group']==seed; ya=a['target'][ix]; yb=b['target'][ix]
        ea=ya>=cutoff; eb=yb>=cutoff; union=np.sum(ea|eb)
        rms=float(np.sqrt(np.mean((ya-yb)**2)))
        row={'seed':seed,'target_RMSE':rms,'relative_target_RMSE':rms/float(np.sqrt(np.mean(yb**2))),
             'event_Jaccard':float(np.sum(ea&eb)/union) if union else 1.,'event_flip_fraction':float(np.mean(ea!=eb)),
             'present_standardized_RMS':float(np.sqrt(np.mean(((a['present'][ix]-b['present'][ix])/sa)**2))),
             'history_standardized_RMS':float(np.sqrt(np.mean(((a['history'][ix]-b['history'][ix])/sh)**2))),
             'present_per_feature_standardized_RMS':np.sqrt(np.mean(((a['present'][ix]-b['present'][ix])/sa)**2,axis=0)).tolist(),
             'history_per_feature_standardized_RMS':np.sqrt(np.mean(((a['history'][ix]-b['history'][ix])/sh)**2,axis=0)).tolist()}
        for name in ['A_full','B_history']:
            for kind in ['prediction','probability']:
                key=name+'_'+kind; coarse=pa[key][ix]; fine=pb[key][ix]
                change=float(np.sqrt(np.mean((coarse-fine)**2)))
                row[key+'_RMS']=change; row[key+'_relative_RMS']=change/max(float(np.sqrt(np.mean(fine*fine))),1e-12)
        out.append(row)
    return out


def metric_drift(a,b,metric,p):
    da=np.array(a['paired'][metric]['paired_run_deltas']); db=np.array(b['paired'][metric]['paired_run_deltas'])
    change=db-da; rng=np.random.default_rng(p['bootstrap']['seed'])
    draw=rng.integers(0,len(change),(p['bootstrap']['replicates'],len(change)))
    return {'mean_fine_minus_coarse':float(change.mean()),'per_run_fine_minus_coarse':change.tolist(),
            'ci95':np.quantile(change[draw].mean(1),[.025,.975]).tolist()}


def main():
    p=json.loads((BASE/'protocol.json').read_text()); q=json.loads((HERE/'protocol.json').read_text())
    execution=json.loads((HERE/'execution.json').read_text())
    assert execution['complete']; grids=execution['grids']
    model_path=BASE/'results'/'primary_models.json'; assert sha(model_path)==q['models']['sha256']
    models=json.loads(model_path.read_text())
    assert models['models']['A_full']['feature_names']==PRESENT_NAMES
    assert models['models']['B_history']['feature_names']==PRESENT_NAMES+HISTORY_NAMES
    assert models['event_target_cutoff']==q['primary']['event_cutoff']
    for key in ['anchors','primary_lookback','primary_horizon','gap']:
        expected=q['primary'][{'primary_lookback':'lookback','primary_horizon':'horizon'}.get(key,key)]
        assert p[key]==expected
    result={'protocol':q,'grids':grids,'models_refitted':False,'frozen_models_sha256':sha(model_path),
            'analysis_sha256':sha(Path(__file__)),'grid_models':{},'spatial':{},'temporal':{},'files':[],
            'screens':[],'history_names':HISTORY_NAMES,'present_names':PRESENT_NAMES,'frozen_history_controls':{}}
    data={}; predictions={}; saved={}
    for n in grids:
        paths=[HERE/rec['path'] for rec in execution['records'] if rec['kind']=='main' and rec['config']['n']==n]
        assert len(paths)==8
        for path in paths:
            m=checked(path,{**q['base'],'n':n})
            with np.load(path,allow_pickle=False) as z:
                div=np.trace(z['gradient'],axis1=-2,axis2=-1)
                result['files'].append({'path':str(path.relative_to(BASE)),'sha256':m['sha256'],'seed':m['seed'],
                    'config':m['config'],'max_CFL':m['max_cfl'],'interpolated_divergence_RMS':float(np.sqrt(np.mean(div*div))),
                    'max_probe_det_error':max(row['max_det_error'] for row in m['diagnostics']),
                    'max_energy_residual':max(abs(row['energy_budget_residual']) for row in m['diagnostics'])})
        d=dataset(paths,p); data[n]=d
        result['grid_models'][str(n)],predictions[n]=score(d,models,p)
        for key in ['target','present','history','position','group','anchor','particle']:saved[f'n{n}_{key}']=d[key]
        for key,value in predictions[n].items():saved[f'n{n}_{key}']=value
    limit=q['screens']; finest=grids[-1]; previous=grids[-2]
    def screen(kind,identifier,value,threshold,comparison='<='):
        ok=value>=threshold if comparison=='>=' else value<threshold if comparison=='<' else value<=threshold
        result['screens'].append({'kind':kind,'identifier':identifier,'value':float(value),'limit':threshold,'comparison':comparison,'pass':bool(ok)})
    prior_rows=None
    for lo,hi in zip(grids[:-1],grids[1:]):
        rows=compare_datasets(data[lo],data[hi],predictions[lo],predictions[hi],models,p)
        drift={metric:metric_drift(result['grid_models'][str(lo)],result['grid_models'][str(hi)],metric,p) for metric in ['RMSE','AUPRC_AP']}
        result['spatial'][f'{lo}-{hi}']={'per_run':rows,'paired_metric_drift':drift}
        for i,row in enumerate(rows):
            key=f's{row["seed"]}/{lo}-{hi}'
            if prior_rows is not None: screen('target_difference_contracts',key,row['target_RMSE']/prior_rows[i]['target_RMSE'],1.,'<')
            if hi==finest:
                screen('target_relative_RMS',key,row['relative_target_RMSE'],limit['finest_pair_target_relative_RMSE_max'])
                screen('event_Jaccard',key,row['event_Jaccard'],limit['finest_pair_event_Jaccard_min'],'>=')
                for name in ['A_full','B_history']:
                    for kind in ['prediction','probability']:
                        screen('prediction_relative_RMS',key+'/'+name+'/'+kind,row[name+'_'+kind+'_relative_RMS'],limit['finest_pair_prediction_relative_RMSE_max'])
        if hi==finest:
            for metric,label in [('RMSE','RMSE'),('AUPRC_AP','AP')]:
                d=drift[metric]
                screen('paired_'+label+'_macro_drift',f'{lo}-{hi}',abs(d['mean_fine_minus_coarse']),limit['paired_'+label+'_difference_macro_drift_max'])
                for seed,change in zip(p['test_seeds'],d['per_run_fine_minus_coarse']):
                    screen('paired_'+label+'_per_run_drift',f's{seed}/{lo}-{hi}',abs(change),limit['paired_'+label+'_difference_per_run_drift_max'])
        prior_rows=rows
    for n in [g for g in q['half_step_grids'] if g in grids]:
        paths=[HERE/rec['path'] for rec in execution['records'] if rec['kind']=='half' and rec['config']['n']==n]
        assert len(paths)==2
        for path in paths: checked(path,{**q['base'],'n':n,'dt':q['half_step_dt']})
        d=dataset(paths,p); evaluated,pred=score(d,models,p)
        # The control subset contains only the two prechosen whole runs.
        ix=np.isin(data[n]['group'],q['half_step_seeds']); sub={k:v[ix] for k,v in data[n].items()}
        subpred={k:v[ix] for k,v in predictions[n].items()}; subp={**p,'test_seeds':q['half_step_seeds']}
        rows=compare_datasets(sub,d,subpred,pred,models,subp)
        main,_=score(sub,models,p); drift={metric:metric_drift(main,evaluated,metric,p) for metric in ['RMSE','AUPRC_AP']}
        result['temporal'][str(n)]={'per_run':rows,'paired_metric_drift':drift}
        if n==finest:
            spatial=result['spatial'][f'{previous}-{finest}']['per_run']
            for row in rows:
                base=next(r for r in spatial if r['seed']==row['seed'])
                screen('time_target_fraction',f's{row["seed"]}',row['target_RMSE']/max(base['target_RMSE'],1e-12),limit['half_step_target_fraction_of_spatial_max'])
            for metric,label in [('RMSE','RMSE'),('AUPRC_AP','AP')]:
                for seed,change in zip(q['half_step_seeds'],drift[metric]['per_run_fine_minus_coarse']):
                    screen('time_paired_'+label+'_drift',f's{seed}',abs(change),limit['half_step_paired_'+label+'_difference_abs_max'])
    d=data[finest]; m=models['models']['B_history']; e=(d['target']>=models['event_target_cutoff']).astype(int)
    controls={'reversed':d['reverse'],'irrelevant':d['irrelevant'],
              'row_permuted':shuffle_blocks(d['history'],d['group'],d['anchor'],718)}
    for label,H in controls.items():
        X=np.column_stack([d['present'],H]); pred=predict_model(m['continuous'],X); prob=predict_model(m['event'],X)
        result['frozen_history_controls'][label]=summarize(d['target'],e,pred,prob,d['group'],m['decision_threshold'])
    assert sha(model_path)==q['models']['sha256']
    result['all_screens_pass']=all(s['pass'] for s in result['screens'])
    result['failures']=[s for s in result['screens'] if not s['pass']]
    save_json(HERE/'results.json',result); np.savez_compressed(HERE/'predictions-and-features.npz',**saved)
    with (HERE/'per-run-comparisons.csv').open('w',newline='') as f:
        fields=['comparison','seed','target_RMSE','relative_target_RMSE','event_Jaccard','event_flip_fraction','present_standardized_RMS','history_standardized_RMS']
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore'); w.writeheader()
        for label,rec in result['spatial'].items():w.writerows({'comparison':label,**row} for row in rec['per_run'])
    plot(result)
    print(json.dumps({'grids':grids,'screens':len(result['screens']),'passed':sum(s['pass'] for s in result['screens']),
                      'failures':result['failures'],'finest_paired':result['grid_models'][str(finest)]['paired']},indent=2))


def plot(r):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    grids=r['grids']; fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    hi=grids[1:]
    for seed in range(8):
        rows=[r['spatial'][f'{a}-{b}']['per_run'][seed] for a,b in zip(grids[:-1],grids[1:])]
        axes[0,0].plot(hi,[100*x['relative_target_RMSE'] for x in rows],'o-',alpha=.65)
        axes[0,1].plot(hi,[x['event_Jaccard'] for x in rows],'o-',alpha=.65)
    axes[0,0].axhline(1,color='black',linestyle='--'); axes[0,0].set(title='Full-history target refinement',ylabel='Relative target RMS (%)')
    axes[0,1].axhline(.95,color='black',linestyle='--');axes[0,1].set(title='Strong-event sets across grids',ylabel='Intersection / union')
    for ax,metric,title in [(axes[1,0],'RMSE','B minus A error; negative favors history'),(axes[1,1],'AUPRC_AP','B minus A average precision; positive favors history')]:
        vals=[r['grid_models'][str(n)]['paired'][metric] for n in grids]; means=np.array([v['mean_B_minus_A'] for v in vals]); ci=np.array([v['ci95'] for v in vals])
        ax.errorbar(grids,means,yerr=np.stack([means-ci[:,0],ci[:,1]-means]),fmt='o-',capsize=4);ax.axhline(0,color='black',linestyle='--');ax.set(title=title)
    for ax in axes.flat:ax.set_xlabel('Finer grid size' if ax in axes[0] else 'Grid size');ax.grid(alpha=.2)
    fig.suptitle('Full-history refinement: eight reused flows, frozen predictors')
    fig.savefig(HERE/'convergence.png',dpi=160);fig.savefig(HERE/'convergence.svg');plt.close(fig)


if __name__=='__main__':main()
