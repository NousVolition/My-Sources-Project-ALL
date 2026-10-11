"""Whole-flow inference for common pulses, forcing feedback and numerical controls.

No incomplete matrix is silently analyzed. Layouts, markers and resolutions
are paired conditions, never independent statistical samples.
"""
import json
import numpy as np
from final_experiment import ROOT,P,dump,digest
from analysis import read,filename,response,measurements,comparison,interval,effective_regimes,distance

METRICS=['retained','scrambled','attenuation','matched_translation','specificity','reset']
def path(folder,seed,regime,n,dt,layout=None):
    middle=f'{regime}-{layout}' if layout is not None else regime
    return ROOT/folder/f'{middle}-s{seed}-n{n}-dt{dt:g}.npz'

def common_response(a,m,t=None):
    t=P['horizon'] if t is None else t
    at=np.flatnonzero(np.isclose(a['time'],t))[0]
    if 'base_tangent' in a:
        kicked=a['tangent'][at];base=a['base_tangent'][at]
    else:
        kicked=a['tangent'][at,1::2];base=a['tangent'][at,::2]
    ftle=lambda f:np.log(np.linalg.svd(f,compute_uv=False)[...,0])/t
    r=ftle(kicked)-ftle(base)
    return {arm:r[2*j:2*j+2] for j,arm in enumerate(P['arms'])}

def common_measure(a,m,t=None):
    r=common_response(a,m,t)
    values={arm:float(np.sqrt(np.mean((v[0]-v[1])**2))) for arm,v in r.items()}
    scrambled=float(np.mean([values['scramble'+str(j)] for j in range(P['scrambles'])]))
    return dict(retained=values['retained'],scrambled=scrambled,attenuation=values['retained']-scrambled,
        matched_translation=values['matched_translation'],specificity=values['matched_translation']-scrambled,reset=values['reset'])

def common_fields(a):
    return a['end_fields']-a['base_end_fields'] if 'base_end_fields' in a else a['end_fields'][1::2]-a['end_fields'][::2]

def common_comparison(low,high):
    a,ma=read(low,['time','tangent','base_tangent','end_fields','base_end_fields','modes','pulse'])
    b,mb=read(high,['time','tangent','base_tangent','end_fields','base_end_fields','modes','pulse'])
    ra,rb=common_response(a,ma),common_response(b,mb)
    diff=np.stack([ra[k]-rb[k] for k in P['arms']])
    field=distance(common_fields(a),a['modes'],common_fields(b),b['modes'])
    pulse=distance(a['pulse'],a['modes'],b['pulse'],b['modes'])
    return dict(response_ftle_rmse=float(np.sqrt(np.mean(diff**2))),response_field_relative=float(field.max()),
        pulse_relative=float(pulse),effect_difference=common_measure(a,ma)['attenuation']-common_measure(b,mb)['attenuation'],
        constraint_rank_low=ma.get('audit',{}).get('constraint_rank'),constraint_rank_high=mb.get('audit',{}).get('constraint_rank'))

def basic_screens(m):
    g=P['gates'];d=m['maxima']
    return {**{k:d[k]<g[k] for k in ['divergence_rms','gradient_trace','energy_budget_relative','tail_enstrophy','tangent_volume_error']},'cfl':d['max_cfl']<g['cfl']}

def endpoint_resolution(a,regime):
    """Post-hoc screen all available start/end fields, not just one preparation."""
    k2=np.sum(a['modes']**2,axis=1);w=np.where(a['modes'][:,2]==0,1.,2.);nu=P['regimes'][regime]['nu'];kmax=int(abs(a['modes']).max())
    eta=[]
    for key in ['start_fields','end_fields','base_start_fields','base_end_fields']:
        if key in a:
            diss=nu*np.sum(abs(a[key])**2*k2*w,axis=(-2,-1))
            eta.extend((kmax*(nu**3/np.maximum(diss,1e-300))**.25).tolist())
    return dict(min_endpoint_kmax_eta=float(min(eta)),passed=regime=='smooth' or min(eta)>=P['gates']['kmax_eta'],scope='Endpoints only; not a bound between recorded times')

def summarize(rows,methods):
    return {method:{metric:interval([r[method][metric] for r in rows]) for metric in METRICS} for method in methods}

def claim(low,high,changes):
    same_sign=low['mean']*high['mean']>0
    excluded=lambda x:x['ci95'][0]>0 or x['ci95'][1]<0
    mean_change=float(np.mean(abs(np.array(changes))))
    scale=abs(high['mean'])
    stable=mean_change<.1*scale
    return dict(ci_excludes_zero_on_both=excluded(low) and excluded(high),same_sign=bool(same_sign),
        mean_absolute_grid_change=mean_change,relative_grid_change=mean_change/max(scale,1e-30),effect_resolution_passed=bool(same_sign and excluded(low) and excluded(high) and stable))

def main_analysis():
    groups={};per=[];grid=[];steps=[];numerical=[];claims={};time=[]
    for regime,cfg in effective_regimes().items():
        seeds=P['smooth_seeds'] if regime=='smooth' else P['forced_seeds']
        for layout in P['layouts']:
            rows=[[],[]]
            for seed in seeds:
                files=[filename(seed,regime,layout,n) for n in cfg['grids']]
                values=[]
                for level,file in enumerate(files):
                    a,m=read(file);v=measurements(a,m);rows[level].append(v);values.append(v)
                    numeric=dict(config=m['config'],screens=basic_screens(m),work_error=m['work_audit']['equal_work_error'],rms_error=m['work_audit']['equal_rms_error'],endpoint_resolution=endpoint_resolution(a,regime),field_stats=m['field_stats'],work_audit=m['work_audit'],matching=m['matching'])
                    numerical.append(numeric)
                    if level==1:
                        for horizon in a['time'][1:]:
                            r=response(a,m,float(horizon))
                            mags={method:{arm:float(np.sqrt(np.mean((z[0]-z[1])**2))) for arm,z in arms.items()} for method,arms in r.items()}
                            time.append(dict(seed=seed,regime=regime,layout=layout,time=float(horizon),values=mags))
                per.append(dict(seed=seed,regime=regime,layout=layout,coarse=values[0],fine=values[1]))
                grid.append(dict(seed=seed,regime=regime,layout=layout,**comparison(*files)))
                if seed in seeds[:2]:
                    half=filename(seed,regime,layout,cfg['grids'][-1],cfg['dt']/2)
                    a,m=read(half);numerical.append(dict(config=m['config'],screens=basic_screens(m),work_error=m['work_audit']['equal_work_error'],rms_error=m['work_audit']['equal_rms_error'],endpoint_resolution=endpoint_resolution(a,regime)))
                    steps.append(dict(seed=seed,regime=regime,layout=layout,**comparison(files[1],half)))
            key=regime+'/'+layout;groups[key]={label:summarize(v,P['probe_methods']) for label,v in zip(['coarse','fine'],rows)}
            claims[key]={method:claim(groups[key]['coarse'][method]['attenuation'],groups[key]['fine'][method]['attenuation'],[a[method]['attenuation']-b[method]['attenuation'] for a,b in zip(*rows)]) for method in P['probe_methods']}
    interactions={}
    for regime in effective_regimes():
        base=[r for r in per if r['regime']==regime and r['layout']=='coincident']
        for layout in ['writer_moved','probe_moved']:
            other=[r for r in per if r['regime']==regime and r['layout']==layout]
            interactions[regime+'/'+layout]={method:interval([b['fine'][method]['attenuation']-a['fine'][method]['attenuation'] for a,b in zip(base,other)]) for method in P['probe_methods']}
    method_interactions={key:{method:interval([r['fine'][method]['attenuation']-r['fine']['fixed_rms']['attenuation'] for r in per if r['regime']+'/'+r['layout']==key]) for method in P['probe_methods'][1:]} for key in groups}
    regime_interactions={layout:{method:interval([next(r for r in per if r['seed']==seed and r['regime']=='forced_nu002' and r['layout']==layout)['fine'][method]['attenuation']-next(r for r in per if r['seed']==seed and r['regime']=='forced_nu004' and r['layout']==layout)['fine'][method]['attenuation'] for seed in P['forced_seeds']]) for method in P['probe_methods']} for layout in P['layouts']}
    result=dict(groups=groups,per_flow=per,grid=grid,half_step=steps,numerical=numerical,claims=claims,layout_interactions=interactions,method_interactions=method_interactions,regime_interactions=regime_interactions,time_course=time)
    dump(ROOT/'results'/'main-analysis.json',result);return result

def common_analysis():
    groups={};per=[];grid=[];steps=[];numerical=[];claims={}
    for regime,cfg in effective_regimes().items():
        seeds=P['smooth_seeds'] if regime=='smooth' else P['forced_seeds']
        for layout in P['layouts'] if regime=='smooth' else ['coincident']:
            rows=[[],[]]
            for seed in seeds:
                files=[path('common-probe-data',seed,regime,n,cfg['dt'],layout) for n in cfg['grids']]
                vals=[]
                for level,file in enumerate(files):
                    a,m=read(file);v=common_measure(a,m);rows[level].append(v);vals.append(v)
                    numerical.append(dict(config=m['config'],screens=basic_screens(m),audit=m['audit'],endpoint_resolution=endpoint_resolution(a,regime)))
                per.append(dict(seed=seed,regime=regime,layout=layout,coarse=vals[0],fine=vals[1]));grid.append(dict(seed=seed,regime=regime,layout=layout,**common_comparison(*files)))
                if (regime=='smooth' and seed in seeds[:2]) or (regime!='smooth' and seed==seeds[0]):
                    half=path('common-probe-data',seed,regime,cfg['grids'][-1],cfg['dt']/2,layout)
                    a,m=read(half);numerical.append(dict(config=m['config'],screens=basic_screens(m),audit=m['audit'],endpoint_resolution=endpoint_resolution(a,regime)))
                    steps.append(dict(seed=seed,regime=regime,layout=layout,**common_comparison(files[1],half)))
            key=regime+'/'+layout
            groups[key]={level:{metric:interval([v[metric] for v in rs]) for metric in METRICS} for level,rs in zip(['coarse','fine'],rows)}
            claims[key]=claim(groups[key]['coarse']['attenuation'],groups[key]['fine']['attenuation'],[a['attenuation']-b['attenuation'] for a,b in zip(*rows)])
    interactions={}
    for layout in ['writer_moved','probe_moved']:
        a=[r for r in per if r['regime']=='smooth' and r['layout']=='coincident'];b=[r for r in per if r['regime']=='smooth' and r['layout']==layout]
        interactions[layout]=interval([y['fine']['attenuation']-x['fine']['attenuation'] for x,y in zip(a,b)])
    result=dict(groups=groups,per_flow=per,grid=grid,half_step=steps,numerical=numerical,claims=claims,layout_interactions=interactions)
    dump(ROOT/'results'/'common-analysis.json',result);return result

def unforced_analysis():
    results={};numeric=[]
    for regime,cfg in effective_regimes().items():
        if regime=='smooth':continue
        seeds=P['forced_seeds'][:2];groups={};paired=[];grids=[];steps=[]
        for method in ['individual','common']:
            rows=[[],[]]
            for seed in seeds:
                vals=[];files=[]
                for level,n in enumerate(cfg['grids']):
                    folder='unforced-readout-data' if method=='individual' else 'common-unforced-data';file=path(folder,seed,regime,n,cfg['dt']);files.append(file);a,m=read(file)
                    values=measurements(a,m) if method=='individual' else {'same_pulse':common_measure(a,m)}
                    rows[level].append(values);vals.append(values);numeric.append(dict(config=m['config'],kind=method,screens=basic_screens(m),endpoint_resolution=endpoint_resolution(a,regime)))
                parent=filename(seed,regime,'coincident',cfg['grids'][-1]) if method=='individual' else path('common-probe-data',seed,regime,cfg['grids'][-1],cfg['dt'],'coincident')
                a,m=read(parent);on=measurements(a,m) if method=='individual' else {'same_pulse':common_measure(a,m)}
                paired.append(dict(seed=seed,method=method,forcing_off=vals[1],forcing_on=on))
                cmp=comparison if method=='individual' else common_comparison
                grids.append(dict(seed=seed,method=method,**cmp(*files)))
                if seed==seeds[0]:
                    half=path('unforced-readout-data' if method=='individual' else 'common-unforced-data',seed,regime,cfg['grids'][-1],cfg['dt']/2)
                    a,m=read(half);numeric.append(dict(config=m['config'],kind=method,screens=basic_screens(m),endpoint_resolution=endpoint_resolution(a,regime)))
                    steps.append(dict(seed=seed,method=method,**cmp(files[1],half)))
            methods=P['probe_methods'] if method=='individual' else ['same_pulse']
            groups[method]={level:summarize(rs,methods) for level,rs in zip(['coarse','fine'],rows)}
            groups[method]['paired_off_minus_on']={m:interval([r['forcing_off'][m]['attenuation']-r['forcing_on'][m]['attenuation'] for r in paired if r['method']==method]) for m in methods}
        results[regime]=dict(groups=groups,paired=paired,grid=grids,half_step=steps,n_flows=2,exploratory=True)
    result=dict(regimes=results,numerical=numeric);dump(ROOT/'results'/'unforced-analysis.json',result);return result

def controls():
    translation=[];extra=[]
    for seed in P['smooth_seeds'][:2]:
        src=filename(seed,'smooth','coincident',32);a,ma=read(src,['time','tangent']);b,mb=read(src.with_name(src.stem+'-translated.npz'),['time','tangent'])
        ra,rb=response(a,ma),response(b,mb);err=max(float(abs(ra[m][arm]-rb[m][arm]).max()) for m in ra for arm in P['arms'])
        translation.append(dict(seed=seed,response_error=err,passed=err<1e-10))
    cfg=effective_regimes()['forced_nu004']
    for seed in P['forced_seeds'][:2]:
        a=filename(seed,'forced_nu004','coincident',cfg['grids'][-1]);b=filename(seed,'forced_nu004','coincident',96)
        ca=path('common-probe-data',seed,'forced_nu004',cfg['grids'][-1],cfg['dt'],'coincident');cb=path('common-probe-data',seed,'forced_nu004',96,cfg['dt'],'coincident')
        extra.append(dict(seed=seed,main=comparison(a,b),common=common_comparison(ca,cb)))
    cpu=[]
    for folder in ['independent-full-data','independent-common-data']:
        for f in sorted((ROOT/folder).glob('*.npz')):
            a,m=read(f,['modes']);cpu.append(dict(folder=folder,file=f.name,**m))
    if len(cpu)!=12:raise RuntimeError('Incomplete independent solver matrix')
    result=dict(translation=translation,finest_moderate_audits=extra,independent_solver=cpu)
    dump(ROOT/'results'/'controls-analysis.json',result);return result

def all_results():
    for name in ['gpu-pipeline-complete.json','cpu-pipeline-complete.json']:
        if not (ROOT/'results'/name).exists():raise RuntimeError('Production incomplete: '+name)
    a=main_analysis();b=common_analysis();c=unforced_analysis();d=controls();g=P['gates']
    checks={}
    for label,result in [('main',a),('common',b)]:
        for regime in effective_regimes():
            nums=[r for r in result['numerical'] if r['config']['regime']==regime]
            grids=[r for r in result['grid'] if r['regime']==regime];steps=[r for r in result['half_step'] if r['regime']==regime]
            checks[label+'/'+regime]=dict(trajectory=all(all(r['screens'].values()) for r in nums),endpoint_resolution=all(r['endpoint_resolution']['passed'] for r in nums),
                work=all((r.get('work_error',r.get('audit',{}).get('work_error',0)))<g['work_error'] for r in nums),
                rms=all((r.get('rms_error',r.get('audit',{}).get('rms_error',0)))<g['rms_error'] for r in nums),
                grid_ftle=all(r['response_ftle_rmse']<g['response_ftle_grid_rmse'] for r in grids),grid_field=all(r['response_field_relative']<g['response_field_grid_relative'] for r in grids),
                half_step=all(r['response_ftle_rmse']<g['response_ftle_step_rmse'] for r in steps))
    for regime,result in c['regimes'].items():
        nums=[r for r in c['numerical'] if r['config']['regime']==regime]
        checks['unforced/'+regime]=dict(trajectory=all(all(r['screens'].values()) and r['endpoint_resolution']['passed'] for r in nums),grid_ftle=all(r['response_ftle_rmse']<g['response_ftle_grid_rmse'] for r in result['grid']),grid_field=all(r['response_field_relative']<g['response_field_grid_relative'] for r in result['grid']),half_step=all(r['response_ftle_rmse']<g['response_ftle_step_rmse'] for r in result['half_step']))
    checks['controls']=dict(translation=all(r['passed'] for r in d['translation']),independent_solver=all(r['passed'] for r in d['independent_solver']),
        moderate_finest_ftle=all(r[k]['response_ftle_rmse']<g['response_ftle_grid_rmse'] for r in d['finest_moderate_audits'] for k in ['main','common']),moderate_finest_field=all(r[k]['response_field_relative']<g['response_field_grid_relative'] for r in d['finest_moderate_audits'] for k in ['main','common']))
    def regime_control_pass(regime):
        cpu=[r for r in d['independent_solver'] if r['regime']==regime]
        passed=bool(cpu) and all(r['passed'] for r in cpu)
        if regime=='smooth':passed=passed and checks['controls']['translation']
        if regime=='forced_nu004':passed=passed and checks['controls']['moderate_finest_ftle'] and checks['controls']['moderate_finest_field']
        return passed
    for key,v in a['claims'].items():
        regime=key.split('/')[0]
        for method,r in v.items():r['numerically_qualified']=all(checks['main/'+regime].values()) and regime_control_pass(regime);r['physical_claim_passed']=r['effect_resolution_passed'] and r['numerically_qualified']
    for key,r in b['claims'].items():
        regime=key.split('/')[0];r['numerically_qualified']=all(checks['common/'+regime].values()) and regime_control_pass(regime);r['physical_claim_passed']=r['effect_resolution_passed'] and r['numerically_qualified']
    # Save augmented claims too, so report and machine-readable files agree.
    dump(ROOT/'results'/'main-analysis.json',a);dump(ROOT/'results'/'common-analysis.json',b)
    result=dict(checks=checks,all_numerical_checks_passed=all(all(v.values()) for v in checks.values()),main_claims=a['claims'],common_claims=b['claims'],
        independent_statistical_units='Whole fresh initial-flow seeds; layouts/methods/markers/grids remain paired',exploratory_intervals_unadjusted=True,
        inference_safeguard_sha256=digest(ROOT/'inference-safeguard.json'),request_id='2300ecd1-d208-4695-9b39-62b93d705892')
    dump(ROOT/'results'/'extension-summary.json',result);print(json.dumps(result,indent=2));return result
if __name__=='__main__':all_results()
