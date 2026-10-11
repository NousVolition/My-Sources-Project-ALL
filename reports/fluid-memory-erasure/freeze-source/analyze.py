"""Whole-flow inference, frozen prediction models and numerical gates."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
from run import ROOT,P,digest,dump,filename
AMENDMENT=json.loads((ROOT/'production-amendment.json').read_text()) if (ROOT/'production-amendment.json').exists() else {}
N_MAIN=AMENDMENT.get('production_n',P['n'])


def load_run(folder,seed,n=None,dt=None,nu=None,amplitude=None,end=None,backend='gpu',keys=None):
    cfg=dict(seed=seed,n=N_MAIN if n is None else n,dt=P['dt'] if dt is None else dt,
             nu=P['nu'] if nu is None else nu,amplitude=P['probe_velocity_rms'] if amplitude is None else amplitude,
             end=max(P['horizons']) if end is None else end,backend=backend)
    stem=filename(**cfg);path=Path(folder)/(stem+'.npz');mpath=path.with_suffix('.json')
    meta=json.loads(mpath.read_text())
    if meta['raw_sha256']!=digest(path):raise RuntimeError('Data hash mismatch: '+str(path))
    if cfg!=meta['config']:raise RuntimeError('Config mismatch')
    for source,sha in meta['sources'].items():
        if digest(ROOT/source)!=sha:raise RuntimeError('Executed source changed: '+source)
    with np.load(path,allow_pickle=False) as archive:
        arrays={key:archive[key] for key in archive.files if keys is None or key in keys}
    return arrays,meta


def branch_index(meta,arm,history,kicked):
    return next(i for i,b in enumerate(meta['branches']) if b['arm']==arm and b['history']==history and b['kicked']==bool(kicked))


def ftle(arrays,horizon):
    ix=np.flatnonzero(np.isclose(arrays['time'],horizon,atol=1e-9))[0]
    return np.log(np.linalg.svd(arrays['tangent'][ix],compute_uv=False)[...,0])/horizon


def response_arrays(arrays,meta,horizon):
    f=ftle(arrays,horizon)
    return {arm:np.stack([f[branch_index(meta,arm,h,1)]-f[branch_index(meta,arm,h,0)] for h in range(2)])
            for arm in dict.fromkeys(b['arm'] for b in meta['branches'])}


def causal_metrics(arrays,meta,horizon):
    r=response_arrays(arrays,meta,horizon)
    contrast={arm:float(np.sqrt(np.mean((value[0]-value[1])**2))) for arm,value in r.items()}
    phases=[contrast['scramble'+str(i)] for i in range(P['phase_scramble_replicates'])]
    delta=contrast['retained']-np.mean(phases)
    reset=max(float(abs(r['reset'][0]-r['reset'][1]).max()),
              float(abs(arrays['start_fields'][branch_index(meta,'reset',0,0)]-arrays['start_fields'][branch_index(meta,'reset',1,0)]).max()))
    # Marker-array resetting is only an observation reset; fluid branches are identical.
    row=dict(primary_attenuation=float(delta),retained_history_contrast=contrast['retained'],
                scrambled_history_contrast=float(np.mean(phases)),scramble_replicates=phases,
                sham_history_contrast=contrast['sham'],sham_attenuation=float(contrast['retained']-contrast['sham']),
                reset_discrepancy=reset,native_history_contrast=contrast['native'],
                retained_signed_response_mean=float(np.mean(r['retained'][0]-r['retained'][1])))
    ri=np.flatnonzero(np.isclose(arrays['response_time'],horizon))[0]
    weight=np.where(arrays['modes'][:,2]==0,1.,2.)
    field_norm=lambda h:np.sqrt(np.sum(abs(h)**2*weight,axis=(-2,-1)))
    field_metrics={}
    for arm in r:
        # Response archive pair order corresponds to unprobed branch index / 2.
        pair=np.stack([arrays['response_fields'][ri,branch_index(meta,arm,h,0)//2] for h in range(2)])
        field_metrics[arm]=dict(gain=float(np.mean(field_norm(pair))/meta['config']['amplitude']),
                               history_difference=float(field_norm(pair[0]-pair[1])/meta['config']['amplitude']))
    row['field_response']=field_metrics
    return row


def matched_sham_study(folder,matched_folder):
    results={};numerical=[];matching=[];rows=[]
    for horizon in P['horizons']:
        values=[];trans=[]
        for seed in P['test_seeds']:
            primary,meta=load_run(folder,seed)
            by_grid=[]
            for n in [N_MAIN,max(P['grids'])]:
                stem=filename(seed,n,P['dt'],P['nu'],P['probe_velocity_rms'],max(P['horizons']),'gpu')+'-matched-sham'
                mp=Path(matched_folder)/(stem+'.json');m=json.loads(mp.read_text());rp=mp.with_suffix('.npz')
                if digest(rp)!=m['raw_sha256']:raise RuntimeError('Matched-sham data hash mismatch')
                with np.load(rp) as z:arrays={k:z[k] for k in z.files}
                f=ftle(arrays,horizon);response=np.stack([f[1]-f[0],f[3]-f[2]])
                magnitude=float(np.sqrt(np.mean((response[0]-response[1])**2)))
                by_grid.append((response,magnitude,m['matching']))
                if horizon==P['primary_horizon']:matching.append(dict(seed=seed,n=n,**m['matching']))
            original=causal_metrics(primary,meta,horizon)
            delta=by_grid[0][1]-original['scrambled_history_contrast'];values.append(delta);trans.append(by_grid[0][1])
            rows.append(dict(seed=seed,horizon=horizon,matched_sham_contrast=by_grid[0][1],scrambled_contrast=original['scrambled_history_contrast'],matched_minus_scrambled=delta))
            if horizon==P['primary_horizon']:
                numerical.append(dict(seed=seed,response_ftle_rmse=float(np.sqrt(np.mean((by_grid[0][0]-by_grid[1][0])**2))),
                                      contrast_difference=by_grid[0][1]-by_grid[1][1],translation_length_difference=by_grid[0][2]['length']-by_grid[1][2]['length']))
        results[str(horizon)]=dict(matched_minus_scrambled=interval(values),matched_sham_history_contrast=interval(trans))
    return dict(horizons=results,per_flow=rows,matching=matching,numerical=numerical,
                interpretation='An exploratory strength-matched structured displacement. Translation also changes location relative to the probe/background; a comparison must not be called proof of a unique memory carrier.')


def interval(values,reps=None,seed=None):
    values=np.asarray(values,float)
    reps=P['inference']['bootstrap_replicates'] if reps is None else reps
    rng=np.random.default_rng(P['inference']['bootstrap_seed'] if seed is None else seed)
    bootstrap=values[rng.integers(0,len(values),(reps,len(values)))].mean(axis=1)
    return dict(mean=float(values.mean()),ci95=np.quantile(bootstrap,[.025,.975]).tolist(),
                n_independent_flows=len(values),per_flow=values.tolist(),sd_between_flows=float(values.std(ddof=1)))


def spectral_distance(a,ka,b,kb):
    """Full relative L2, including fine modes absent from the coarse field."""
    lookup={tuple(k):i for i,k in enumerate(kb)}
    indexes=np.array([lookup[tuple(k)] for k in ka])
    weighta=np.where(ka[:,2]==0,1.,2.);weightb=np.where(kb[:,2]==0,1.,2.)
    eb=np.sum(abs(b)**2*weightb,axis=(-2,-1))
    delta=np.sum(abs(a-b[...,indexes])**2*weighta,axis=(-2,-1))
    missing=np.ones(len(kb),bool);missing[indexes]=False
    delta+=np.sum(abs(b[...,missing])**2*weightb[missing],axis=(-2,-1))
    return np.sqrt(delta/np.maximum(eb,1e-30))


def numerical_comparison(coarse,fine,horizon):
    a,ma=coarse;b,mb=fine
    fa,fb=ftle(a,horizon),ftle(b,horizon)
    ra,rb=response_arrays(a,ma,horizon),response_arrays(b,mb,horizon)
    arms=[v for v in ra if v!='reset']
    response=np.concatenate([ra[v]-rb[v] for v in arms])
    ix=np.flatnonzero(np.isclose(a['response_time'],horizon))[0]
    iy=np.flatnonzero(np.isclose(b['response_time'],horizon))[0]
    distances=spectral_distance(a['response_fields'][ix],a['modes'],b['response_fields'][iy],b['modes'])
    metricsa=causal_metrics(a,ma,horizon);metricsb=causal_metrics(b,mb,horizon)
    return dict(ftle_rmse=float(np.sqrt(np.mean((fa-fb)**2))),
                response_ftle_rmse=float(np.sqrt(np.mean(response**2))),
                max_probe_response_field_relative_L2=float(distances.max()),
                primary_effect_difference=metricsa['primary_attenuation']-metricsb['primary_attenuation'])


def numerical_audit(folder,pilot=False):
    seeds=P['pilot_seeds'] if pilot else P['test_seeds']
    grid=[];time=[];maxima={};source_sets=set();raw=[]
    all_meta=[]
    for seed in seeds:
        runs=[load_run(folder,seed,n=n) for n in P['grids']]
        per=[]
        for a,b in zip(runs[:-1],runs[1:]):
            row=dict(seed=seed,coarse_n=a[1]['config']['n'],fine_n=b[1]['config']['n'])
            row.update(numerical_comparison(a,b,P['primary_horizon']));grid.append(row);per.append(row)
        if not pilot:
            for n in [N_MAIN]+([P['grids'][-1]] if seed in seeds[:2] else []):
                time.append(dict(seed=seed,n=n,**numerical_comparison(load_run(folder,seed,n=n),load_run(folder,seed,n=n,dt=P['dt']/2),P['primary_horizon'])))
        all_meta.extend(r[1] for r in runs)
    # Audit every completed configuration, including dependent sensitivity controls.
    for path in Path(folder).glob('s*.json'):
        meta=json.loads(path.read_text());raw_path=path.with_suffix('.npz')
        if digest(raw_path)!=meta['raw_sha256']:raise RuntimeError('Bad SHA256')
        all_meta.append(meta)
        raw.append(dict(file=raw_path.name,sha256=meta['raw_sha256'],bytes=raw_path.stat().st_size,config=meta['config']))
        source_sets.add(json.dumps(meta['sources'],sort_keys=True))
    for meta in all_meta:
        for key,value in meta['maxima'].items():maxima[key]=max(maxima.get(key,0),abs(value))
    g=P['numerical_gates']
    finest=[r for r in grid if r['fine_n']==max(P['grids'])]
    decreases=[]
    for seed in seeds:
        a,b=[r for r in grid if r['seed']==seed]
        decreases.append(all(b[key]<=max(g['tiny_difference_floor'],g['required_grid_difference_decrease_factor']*a[key])
                             for key in ['ftle_rmse','response_ftle_rmse','max_probe_response_field_relative_L2']))
    gates=dict(divergence=maxima['divergence_rms']<g['divergence_rms'],gradient_trace=maxima['gradient_trace']<g['marker_gradient_trace'],
               energy_budget=maxima['energy_budget_relative']<g['relative_energy_budget'],cfl=maxima['max_cfl']<g['cfl'],
               tail=maxima['tail_enstrophy']<g['tail_enstrophy_fraction'],volume=maxima['tangent_volume_error']<g['tangent_determinant_error'],
               finest_ftle=all(r['ftle_rmse']<g['finest_pair_ftle_rmse'] for r in finest),
               finest_response=all(r['response_ftle_rmse']<g['finest_pair_response_ftle_rmse'] for r in finest),
               finest_response_fields=all(r['max_probe_response_field_relative_L2']<g['finest_pair_probe_response_field_relative_L2'] for r in finest),
               three_grid_decrease=all(decreases),single_executed_source=len(source_sets)==1)
    if time:gates['time_response']=all(r['response_ftle_rmse']<g['time_response_ftle_rmse'] for r in time)
    return dict(scope='pilot' if pilot else 'test-flow grid/time plus all-recording integrity',grid=grid,time=time,maxima=maxima,gates=gates,
                all_passed=all(gates.values()),verified_recordings=raw,count=len(raw))


def wrap(x):return (x+np.pi)%(2*np.pi)-np.pi


def neighbor_offsets(x):
    c=P['centers'];nn=P['neighbor_probes_per_center']
    return wrap(x[c:].reshape(c,nn,3)-x[:c,None,:])


def geometry(e):
    r=np.linalg.norm(e,axis=-1);cov=np.einsum('cki,ckj->cij',e,e)/e.shape[1]
    tri=np.triu_indices(3)
    return np.column_stack([r.mean(1),r.std(1),cov[:,tri[0],tri[1]],np.linalg.eigvalsh(cov)])


def features(past_time,past_positions,current_velocity,current_gradient):
    """Accepts past observations only; no target, future state or tangent argument."""
    if not np.isclose(past_time[-1],P['prepare_end'],atol=1e-12) or not np.all(np.diff(past_time)>0):
        raise ValueError('Strictly ordered past-only observations ending at 0.4 required')
    x=past_positions[-1];e=neighbor_offsets(x);r=np.linalg.norm(e,axis=-1);c=P['centers']
    v=current_velocity;g=current_gradient
    pairv=v[c:].reshape(c,6,3)-v[:c,None]
    rate=np.sum(pairv*e,axis=-1)/np.maximum(r*r,1e-30)
    strain=(g+g.swapaxes(-1,-2))/2;rotation=(g-g.swapaxes(-1,-2))/2
    geo=geometry(e)
    current=np.column_stack([v[:c],g.reshape(c,9),np.linalg.eigvalsh(strain),np.linalg.norm(strain,axis=(-2,-1)),
                              np.linalg.norm(rotation,axis=(-2,-1)),geo,rate.mean(1),rate.std(1),rate.min(1),rate.max(1)])
    history=[]
    for lag in [.1,.2,.3]:
        ix=np.flatnonzero(np.isclose(past_time,past_time[-1]-lag,atol=1e-9))[0]
        old=neighbor_offsets(past_positions[ix]);ro=np.linalg.norm(old,axis=-1)
        logs=np.log(r/np.maximum(ro,1e-30))/lag
        cosine=np.sum(e*old,axis=-1)/np.maximum(r*ro,1e-30)
        history.append(np.column_stack([logs.mean(1),logs.std(1),logs.min(1),logs.max(1),cosine.mean(1),(geo-geometry(old))/lag]))
    return current,np.column_stack(history)


def dataset(folder,seeds,horizon,n=None):
    aa,bb,yy,groups=[],[],[],[]
    for seed in seeds:
        z,meta=load_run(folder,seed,n=n,keys=['time','tangent','prep_time','prep_positions','prep_velocity','prep_gradient'])
        targets=ftle(z,horizon);resp=response_arrays(z,meta,horizon)
        for history in range(2):
            a,b=features(z['prep_time'],z['prep_positions'][:,history],z['prep_velocity'][-1,history],z['prep_gradient'][-1,history])
            aa.append(a);bb.append(b)
            yy.append(np.column_stack([targets[branch_index(meta,'native',history,0)],resp['native'][history]]))
            groups.extend([seed]*len(a))
    return dict(A=np.concatenate(aa),H=np.concatenate(bb),Y=np.concatenate(yy),groups=np.array(groups))


def expansion(a):
    i,j=np.triu_indices(a.shape[1]);return np.column_stack([a,a[:,i]*a[:,j]])


def fit_model(train,val,column,family):
    am=train['A'].mean(0);astd=train['A'].std(0);astd[astd<1e-10]=1
    hm=train['H'].mean(0);hstd=train['H'].std(0);hstd[hstd<1e-10]=1
    rng=np.random.default_rng(703)
    triplets=rng.integers(0,train['A'].shape[1],(train['H'].shape[1],3))
    def transform(d):
        a=(d['A']-am)/astd;cur=expansion(a);h=(d['H']-hm)/hstd
        if family in ['history','shuffled_history','across_flow_history']:return np.column_stack([cur,h])
        if family=='capacity_current':return np.column_stack([cur,np.prod(a[:,triplets],axis=-1)])
        return cur
    tr,va=transform(train),transform(val)
    mean=tr.mean(0);std=tr.std(0);std[std<1e-10]=1
    tr=(tr-mean)/std;va=(va-mean)/std;y=train['Y'][:,column];ym=float(y.mean())
    with threadpool_limits(limits=1):u,s,vt=svd(tr,full_matrices=False,check_finite=False)
    uy=u.T@(y-ym);best=None
    for alpha in P['prediction']['alphas']:
        coef=vt.T@(s/(s*s+alpha)*uy)
        score=float(np.mean((va@coef+ym-val['Y'][:,column])**2))
        if best is None or score<best[0]:best=score,alpha,coef
    return dict(family=family,column=column,alpha=best[1],validation_rmse=float(np.sqrt(best[0])),A_mean=am.tolist(),A_scale=astd.tolist(),
                H_mean=hm.tolist(),H_scale=hstd.tolist(),triplets=triplets.tolist(),expanded_mean=mean.tolist(),expanded_scale=std.tolist(),
                coefficient=best[2].tolist(),intercept=ym,input_count=len(best[2]),training_seeds=P['train_seeds'],validation_seeds=P['validation_seeds'])


def predict_model(model,data):
    a=(data['A']-np.array(model['A_mean']))/np.array(model['A_scale']);cur=expansion(a)
    if model['family'] in ['history','shuffled_history','across_flow_history']:
        h=(data['H']-np.array(model['H_mean']))/np.array(model['H_scale']);x=np.column_stack([cur,h])
    elif model['family']=='capacity_current':x=np.column_stack([cur,np.prod(a[:,np.array(model['triplets'])],axis=-1)])
    else:x=cur
    return ((x-np.array(model['expanded_mean']))/np.array(model['expanded_scale']))@np.array(model['coefficient'])+model['intercept']


def shuffled(d,across=False):
    out={k:v.copy() for k,v in d.items()};rng=np.random.default_rng(905)
    if across:
        seeds=np.unique(d['groups']);other=np.roll(seeds,1)
        for s,os in zip(seeds,other):out['H'][d['groups']==s]=d['H'][d['groups']==os]
    else:
        for seed in np.unique(d['groups']):
            ix=np.flatnonzero(d['groups']==seed);out['H'][ix]=d['H'][rng.permutation(ix)]
    return out


def rmse_by_flow(data,pred,column):
    return np.array([np.sqrt(np.mean((pred[data['groups']==s]-data['Y'][data['groups']==s,column])**2)) for s in np.unique(data['groups'])])


def prediction_study(folder,output):
    scores={};saved={};model_path=output/'frozen-models.json'
    models=json.loads(model_path.read_text()) if model_path.exists() else {}
    for horizon in P['horizons']:
        train=dataset(folder,P['train_seeds'],horizon);val=dataset(folder,P['validation_seeds'],horizon)
        # Freeze models to disk BEFORE loading test observations or targets.
        for column,target in enumerate(P['prediction']['targets']):
            key=f'{target}:T{horizon:g}';models.setdefault(key,{})
            for family in ['current','history','capacity_current','shuffled_history','across_flow_history']:
                if family in models[key]:
                    if models[key][family]['training_seeds']!=P['train_seeds'] or models[key][family]['validation_seeds']!=P['validation_seeds']:
                        raise RuntimeError('Frozen model split mismatch')
                    continue
                tr,va=train,val
                if family=='shuffled_history':tr,va=shuffled(train),shuffled(val)
                if family=='across_flow_history':tr,va=shuffled(train,True),shuffled(val,True)
                models[key][family]=fit_model(tr,va,column,family)
        dump(output/'frozen-models.json',models)
        test=dataset(folder,P['test_seeds'],horizon)
        for column,target in enumerate(P['prediction']['targets']):
            key=f'{target}:T{horizon:g}';scores[key]={};rmse={}
            for family,model in models[key].items():
                te=shuffled(test) if family=='shuffled_history' else shuffled(test,True) if family=='across_flow_history' else test
                pred=predict_model(model,te);errors=rmse_by_flow(test,pred,column);rmse[family]=errors
                scores[key][family]=dict(**interval(errors),alpha=model['alpha'],validation_rmse=model['validation_rmse'],input_count=model['input_count'])
                saved[key.replace(' ','_')+'_'+family]=pred
            for control in ['current','capacity_current','shuffled_history','across_flow_history']:
                scores[key]['history_minus_'+control]=interval(rmse['history']-rmse[control])
            scores[key]['relative_rmse_benefit_percent']=float(100*(rmse['current'].mean()-rmse['history'].mean())/rmse['current'].mean())
            refined=dataset(folder,P['test_seeds'],horizon,n=max(P['grids']))
            refined_errors={}
            for family in ['current','history','capacity_current']:
                pr=predict_model(models[key][family],refined);refined_errors[family]=rmse_by_flow(refined,pr,column)
                saved[key.replace(' ','_')+'_refined_'+family]=pr
            scores[key]['frozen_finer_grid']={family:interval(error) for family,error in refined_errors.items()}
            scores[key]['frozen_finer_grid']['history_minus_current']=interval(refined_errors['history']-refined_errors['current'])
            scores[key]['frozen_finer_grid']['history_minus_capacity_current']=interval(refined_errors['history']-refined_errors['capacity_current'])
    np.savez_compressed(output/'heldout-predictions.npz',**saved)
    return scores


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');ap.add_argument('--out',type=Path,default=ROOT/'results');ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    audit=numerical_audit(args.data,args.pilot);dump(args.out/'numerical-audit.json',audit)
    if args.pilot:
        print(json.dumps(dict(all_passed=audit['all_passed'],gates=audit['gates'],grid=audit['grid']),indent=2));return
    metrics={};rows=[]
    for horizon in P['horizons']:
        row=[dict(seed=seed,horizon=horizon,**causal_metrics(*load_run(args.data,seed),horizon)) for seed in P['test_seeds']]
        rows.extend(row)
        metrics[str(horizon)]={key:interval([r[key] for r in row]) for key in ['primary_attenuation','retained_history_contrast','scrambled_history_contrast','sham_attenuation','reset_discrepancy','retained_signed_response_mean']}
    primary=metrics[str(P['primary_horizon'])]['primary_attenuation']
    finest=[r for r in audit['grid'] if r['fine_n']==max(P['grids'])]
    allgrid=[r for r in audit['grid'] if r['coarse_n']>=N_MAIN]
    numerical_size=float(np.mean([abs(r['primary_effect_difference']) for r in finest]))
    coarse_bias=float(np.mean([sum(abs(r['primary_effect_difference']) for r in allgrid if r['seed']==s) for s in P['test_seeds']]))
    claim_gate=coarse_bias < P['numerical_gates']['numerical_uncertainty_fraction_of_claimed_primary_effect']*abs(primary['mean'])
    sensitivity=[]
    for meta_path in args.data.glob('s*.json'):
        meta=json.loads(meta_path.read_text());cfg=meta['config']
        if cfg['seed'] in P['test_seeds'] and (cfg['end']>.8 or cfg['nu']!=P['nu'] or cfg['amplitude']!=P['probe_velocity_rms']):
            z,m=load_run(args.data,**{k:v for k,v in cfg.items() if k!='seed'},seed=cfg['seed'])
            horizons=[t for t in P['horizons']+P['long_horizons'] if t<=cfg['end']]
            for t in horizons:sensitivity.append(dict(config=cfg,horizon=t,**causal_metrics(z,m,t)))
    result=dict(primary_horizon=P['primary_horizon'],primary=primary,horizons=metrics,per_flow=rows,sensitivity=sensitivity,
                finest_grid_effect_change_mean_absolute=numerical_size,coarse_to_finest_effect_change_bound=coarse_bias,
                numerical_effect_gate=claim_gate,numerical_screens_pass=audit['all_passed'],
                attenuation_supported=bool(primary['ci95'][0]>0 and audit['all_passed'] and claim_gate),
                interpretation='Phase organization candidate only. Same low modes and scalar spectra do not equal a complete present state. Reset identity is expected, not evidence of a new physical memory law.')
    dump(args.out/'causal-results.json',result)
    matched_folder=ROOT/'matched-sham-data'
    if (matched_folder/'execution.json').exists():
        execution=json.loads((matched_folder/'execution.json').read_text())
        if execution['completed']==execution['planned']:
            dump(args.out/'matched-sham-results.json',matched_sham_study(args.data,matched_folder))
    prediction=prediction_study(args.data,args.out);dump(args.out/'prediction-results.json',prediction)
    dump(args.out/'analysis-provenance.json',dict(sources={p.name:digest(p) for p in ROOT.glob('*.py')},protocol_sha256=digest(ROOT/'protocol.json'),raw_count=audit['count']))
    print(json.dumps(dict(primary=primary,numerical_screens_pass=audit['all_passed'],numerical_effect_gate=claim_gate,attenuation_supported=result['attenuation_supported']),indent=2))

if __name__=='__main__':main()
