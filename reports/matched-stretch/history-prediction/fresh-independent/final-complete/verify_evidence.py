"""Read-only final evidence reconstruction; never simulate, fit or alter study files."""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,hashlib,datetime,argparse,tempfile
import numpy as np
HERE=Path(__file__).resolve().parent
F=HERE.parent
sys.path.insert(0,str(F))
import fit_models as fm
import record_history as rh
import numeric_audit as na
from evaluate import summarize
from fresh_statistics import spatial_scores
from analyze_test import evaluate_case

def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(exact=False):
    tol=0. if exact else 1e-12
    diffs=[]; manual_diffs=[]; probability_diffs=[]
    def same(a,b,path='root'):
        if isinstance(a,dict):
            assert a.keys()==b.keys(),path
            for k in a: same(a[k],b[k],path+'.'+str(k))
        elif isinstance(a,(list,tuple)):
            assert len(a)==len(b),path
            for i,(x,y) in enumerate(zip(a,b)): same(x,y,path+'.'+str(i))
        elif isinstance(a,(float,np.floating)):
            err=float(abs(a-b)); assert np.isfinite(a) and np.isfinite(b) and err<=tol,(path,a,b,err)
            diffs.append(err)
        else: assert a==b,(path,a,b)
    def arrays(a,b,label):
        assert a.shape==b.shape,label
        err=float(np.max(abs(a-b))) if a.size else 0.
        assert np.isfinite(a).all() and np.isfinite(b).all() and err<=tol,(label,err)
        diffs.append(err)
        return err
    p=read(F/'protocol.json');lock=read(F/'test-lock.json')
    completion=read(F/'completion.json');auto=read(F/'verification.json')
    r=read(F/'results.json');audit=read(F/'test-numerical-audit.json')
    assert completion['complete'] and auto['verified']
    for key,name in [('results_sha256','results.json'),('verification_sha256','verification.json'),('report_sha256','README.md')]:
        assert completion[key]==sha(F/name)
    assert sha(F/'protocol.json')==lock['protocol_sha256']==r['protocol_sha256']==auto['protocol_sha256']
    mp=F/lock['manifest']; assert sha(mp)==lock['manifest_sha256']==r['manifest_sha256']
    assert sha(F/'model-audit-112.json')==lock['model_numeric_audit_sha256']
    assert sha(F/'test-lock.json')==r['test_lock_sha256']
    publication=read(F/'model-lock-publication.json')
    assert publication['checks_success'] and publication['test_lock_sha256']==sha(F/'test-lock.json')
    assert publication['commit']=='7bc2059a899552cc882b4cb81d27d0e2b6b115c5'
    m,models=fm.load_models(mp)
    assert m['training_only_scalers'] and m['no_train_validation_refit'] and m['test_absent_at_freeze'] and lock['test_absent_at_lock']
    inv=read(F/'pilot-complete/source-inventory.json')
    for group,base in [('fresh_runtime',F.parent.parent),('fresh_implementation',F),('predecessor',F.parent)]:
        for rel,digest in inv[group].items(): assert rh.source_hash_matches((base/rel).read_bytes(),digest),rel
    sets=[set(p[k]) for k in ['train_seeds','validation_seeds','test_seeds']]+[set(p['pilot']['seeds'])]
    assert all(not sets[i]&sets[j] for i in range(4) for j in range(i+1,4))
    records=read(F/'execution.json')['records']
    assert len(records)==146 and len({x['name'] for x in records})==146<=p['max_new_simulations']
    expected={}
    for seeds,n,dt,nu in [(p['test_seeds'],112,.0025,.02),(p['test_numerical_controls']['seeds'],160,.0025,.02),
                          (p['test_numerical_controls']['seeds'],112,.00125,.02)]+[
                          (p['viscosity_controls']['seeds'],112,.0025,nu) for nu in p['viscosity_controls']['nu']]:
        for seed in seeds:
            config={**p['base'],'n':n,'dt':dt,'nu':nu}
            expected[rh.name_for(seed,config)]=(seed,config)
    new=[x for x in records if x['kind']=='untouched-tests-and-controls']
    assert {x['name'] for x in new}==set(expected) and len(new)==76
    archives=[];diagnostics={'max_cfl':0.,'max_energy_residual':0.,'max_divergence':0.,'max_volume_error':0.}
    for row in records:
        path=F/'recorded-data'/(row['name']+'.npz')
        meta=rh.validate(path,row['config']); assert meta['sha256']==row['sha256']
        assert meta['study_sources']==lock['sources']
        if row['name'] in expected:
            seed,config=expected[row['name']]
            assert meta['seed']==seed and meta['config']==config
            assert lock['locked_at_UTC']<meta['started_UTC']<meta['completed_UTC']
            assert meta['started_UTC']>='2026-10-15T06:55:12'
            diagnostics['max_cfl']=max(diagnostics['max_cfl'],meta['max_cfl'])
            diagnostics['max_energy_residual']=max(diagnostics['max_energy_residual'],max(abs(d['energy_budget_residual']) for d in meta['diagnostics']))
            diagnostics['max_divergence']=max(diagnostics['max_divergence'],max(d['divergence_max'] for d in meta['diagnostics']))
            diagnostics['max_volume_error']=max(diagnostics['max_volume_error'],max(d['max_det_error'] for d in meta['diagnostics']))
        archives.append({'path':path.relative_to(F).as_posix(),'sha256':sha(path),'metadata_sha256':sha(path.with_suffix('.json'))})
    print('146 archives/configurations/sources/times/initial markers/tangents verified; 76 exact preregistered new configurations',flush=True)
    capture=[];original_save=na.save_json
    try:
        na.save_json=lambda path,value:capture.append(value)
        reconstructed=na.audit(112,'test')
    finally: na.save_json=original_save
    assert len(capture)==1
    same(reconstructed,audit,'test_numeric_audit')
    limits=p['numeric_screens']
    for s in audit['screens']:
        metric=s['metric']
        if metric=='target_relative_RMS': key='half_step_target_relative_RMS_max' if ('temporal' in s['comparison'] or 'half-step' in s['comparison']) else 'target_relative_RMS_max'
        elif metric=='event_Jaccard': key='event_Jaccard_min'
        elif metric.endswith('_new_probability_relative_RMS'): key='new_model_probability_relative_RMS_max'
        elif metric=='spatial_AP_contrast_drift': key='new_model_spatial_AP_contrast_macro_drift_max' if s['seed']=='macro' else 'new_model_spatial_AP_contrast_per_run_drift_max'
        elif metric=='RMSE_contrast_drift_relative_validation_A': key='new_model_RMSE_contrast_relative_validation_A_max'
        elif metric.startswith(('A_full_','B_history_')): key='frozen_old_prediction_relative_RMS_max'
        else: key=metric+'_max'
        assert s['bound']==limits[key],(metric,key)
        assert s['pass']==(s['value']>=s['bound'] if s['direction']=='min' else s['value']<=s['bound'])
    assert len(audit['screens'])==108 and all(s['pass'] for s in audit['screens'])==audit['all_screens_pass']==r['numerically_audited']
    print('108/108 test screens reconstructed; all actual protocol bounds checked',flush=True)
    def independent_pair(a,b,expected_pair):
        keys=sorted(a['per_run'],key=int);assert keys==expected_pair['group_order']
        quantities={key:np.array([b['per_run'][g][key]-a['per_run'][g][key] for g in keys]) for key in a['macro']}
        quantities['spatial_AP']=np.array([b['spatial'][g]['mean_AP']-a['spatial'][g]['mean_AP'] for g in keys])
        quantities['relative_RMSE_reduction']=np.array([1-b['per_run'][g]['RMSE']/a['per_run'][g]['RMSE'] for g in keys])
        draw=np.random.default_rng(p['bootstrap']['seed']).integers(0,len(keys),(p['bootstrap']['replicates'],len(keys)))
        for metric,values in quantities.items():
            value={'mean':float(values.mean()),'ci95':np.quantile(values[draw].mean(1),[.025,.975]).tolist(),
                   'per_run':values.tolist(),'independent_runs':len(keys)}
            same(value,expected_pair[metric],'paired.'+metric)
    paths=fm.paths_for(p['test_seeds'],112,p)
    assert r['test_files']=={x.name:sha(x) for x in paths}
    case_counts={};no_event={}
    for name,case in m['cases'].items():
        assert set(case['models'])==set(fm.model_names(name=='primary'))
        d=fm.make_dataset(paths,p,case['case']['options']);xs=fm.make_matrices(d)
        with np.load(F/(name+'-predictions.npz'),allow_pickle=False) as saved:
            for key in ['target','group','anchor','particle','position']: arrays(d[key],saved[key],name+'.'+key)
            event=(saved['target']>=case['train_event_cutoff']).astype(int)
            np.testing.assert_array_equal(saved['event'],event)
            assert set(np.unique(saved['group']))==set(p['test_seeds'])
            rebuilt={}
            for model,pair in models[name].items():
                key=model[5:] if model.startswith('tree_') else model
                pred=pair[0].predict(xs[key]);prob=pair[1].predict_proba(xs[key])[:,1]
                arrays(pred,saved[model+'_prediction'],name+'.'+model+'.prediction')
                probability_diffs.append(arrays(prob,saved[model+'_probability'],name+'.'+model+'.probability'))
                info=summarize(saved['target'],event,saved[model+'_prediction'],saved[model+'_probability'],saved['group'],case['models'][model]['decision_threshold'])
                info['spatial']=spatial_scores(event,saved[model+'_probability'],saved['group'],saved['anchor'])
                same(info,r['cases'][name]['models'][model],name+'.'+model+'.metrics')
                rebuilt[model]=info
                # Independent AP formula checks ties and the no-positive convention at every run/anchor.
                for g,item in info['spatial'].items():
                    for row in item['anchors']:
                        use=(saved['group']==int(g))&(saved['anchor']==row['anchor'])
                        y=event[use];q=saved[model+'_probability'][use];order=np.argsort(-q,kind='stable')
                        yy=y[order];qq=q[order];ends=np.r_[np.flatnonzero(np.diff(qq)),len(qq)-1]
                        tp=np.cumsum(yy)[ends];precision=tp/(ends+1)
                        ap=float(np.sum(np.diff(np.r_[0,tp])*precision)/y.sum()) if y.sum() else 0.
                        err=abs(ap-row['AP']);assert err<=1e-14,(name,model,g,row['anchor'],err);manual_diffs.append(err)
            expected_result=r['cases'][name]
            assert expected_result['n_rows']==len(d['target']) and expected_result['n_independent_runs']==64
            assert expected_result['cutoff']==case['train_event_cutoff']
            for model in rebuilt:
                if model!='A_full': independent_pair(rebuilt['A_full'],rebuilt[model],expected_result['paired'][model+'_vs_A_full'])
            for pair_name,a,b in [('B_vs_matched_current','A_matched_current','B_history'),('tree_B_vs_tree_A','tree_A_full','tree_B_history'),('B_weak_vs_A_weak','A_weak','B_weak')]:
                if a in rebuilt: independent_pair(rebuilt[a],rebuilt[b],expected_result['paired'][pair_name])
            no_event[name]=sum(row['events']==0 for v in rebuilt['A_full']['spatial'].values() for row in v['anchors'])
        case_counts[name]=len(models[name]);print(name+': saved rows, frozen predictions, all metrics and whole-run intervals verified',flush=True)
    for nu in p['viscosity_controls']['nu']:
        vp=fm.paths_for(p['viscosity_controls']['seeds'],112,p,nu=nu)
        value,_=evaluate_case(vp,p,m['cases']['primary'],{k:models['primary'][k] for k in ['A_full','B_history']})
        value['dependent_variants_of_test_groups']=True
        same(value,r['viscosity'][str(nu)],'viscosity.'+str(nu))
    primary=r['cases']['primary']['paired']['B_history_vs_A_full']
    assert primary['spatial_AP']['ci95'][1]<p['primary_endpoint']['useful_absolute_AP_gain']
    assert primary['relative_RMSE_reduction']['ci95'][0]>p['secondary_continuous']['useful_gain']
    assert r['primary_conclusion']=='useful-size benefit ruled out for the predeclared comparison; small positive effect remains possible'
    assert r['secondary_continuous_conclusion']=='useful benefit supported for the predeclared comparison'
    assert sha(F/'results.json')==completion['results_sha256'] and sha(mp)==lock['manifest_sha256']
    return {'verified_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_checks_pass':True,
        'mode':'exact Windows' if exact else 'cross-platform absolute tolerance 1e-12',
        'archives_checked':146,'new_archives_checked':76,'independent_primary_tests':64,'cases':case_counts,
        'model_pairs':sum(case_counts.values()),'test_numeric_screens_passed':108,
        'reconstruction_max_absolute_difference':max(diffs),'manual_AP_formula_max_difference':max(manual_diffs),
        'frozen_probability_max_absolute_difference':max(probability_diffs),
        'no_event_anchors':no_event,'new_archive_diagnostics':diagnostics,
        'all_preregistered_metrics_controls_and_whole_run_intervals_reproduced':True,
        'primary':primary,'protocol_sha256':sha(F/'protocol.json'),'results_sha256':sha(F/'results.json'),
        'model_manifest_sha256':sha(mp),'test_lock_sha256':sha(F/'test-lock.json'),
        'test_audit_sha256':sha(F/'test-numerical-audit.json'),'verification_code_sha256':sha(Path(__file__)),
        'archives':archives}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--exact',action='store_true');ap.add_argument('--output',type=Path)
    args=ap.parse_args();result=verify(args.exact)
    if args.output:
        with args.output.open('x',encoding='utf-8',newline='\n') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['archives','primary']}),flush=True)
