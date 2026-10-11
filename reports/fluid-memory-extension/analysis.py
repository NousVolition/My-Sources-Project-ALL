"""Independent-flow inference and separate numerical uncertainty checks."""
import json
from pathlib import Path
import numpy as np
from experiment import ROOT,P,digest,dump,sources

def read(path,keys=None):
    path=Path(path);meta=json.loads(path.with_suffix('.json').read_text())
    if any(digest(ROOT/f)!=sha for f,sha in meta['sources'].items()) or digest(path)!=meta['raw_sha256']:raise RuntimeError('Archive/source integrity failure: '+str(path))
    with np.load(path,allow_pickle=False) as z:a={k:z[k] for k in z.files if keys is None or k in keys}
    return a,meta
def idx(meta,arm,hist,method):return next(j for j,b in enumerate(meta['branches']) if b['arm']==arm and b['history']==hist and b['method']==method)
def response(a,m,horizon=None):
    horizon=P['horizon'] if horizon is None else horizon
    at=np.flatnonzero(np.isclose(a['time'],horizon))[0]
    f=np.log(np.linalg.svd(a['tangent'][at],compute_uv=False)[...,0])/horizon
    return {method:{arm:np.stack([f[idx(m,arm,h,method)]-f[idx(m,arm,h,'unprobed')] for h in range(2)]) for arm in P['arms']} for method in P['probe_methods']}
def measurements(a,m):
    r=response(a,m);out={}
    for method,v in r.items():
        mags={arm:float(np.sqrt(np.mean((x[0]-x[1])**2))) for arm,x in v.items()}
        sc=np.mean([mags['scramble'+str(j)] for j in range(P['scrambles'])])
        out[method]=dict(retained=mags['retained'],scrambled=float(sc),attenuation=float(mags['retained']-sc),
                        matched_translation=mags['matched_translation'],specificity=float(mags['matched_translation']-sc),reset=mags['reset'])
    return out
def interval(x):
    x=np.array(x,float);rng=np.random.default_rng(P['bootstrap_seed']);res=x[rng.integers(len(x),size=(P['bootstrap_replicates'],len(x)))].mean(axis=1)
    return dict(mean=float(x.mean()),ci95=np.quantile(res,[.025,.975]).tolist(),n_flows=len(x),sd=float(x.std(ddof=1)),per_flow=x.tolist())
def field_response(a,m):
    return np.stack([a['end_fields'][idx(m,arm,h,method)]-a['end_fields'][idx(m,arm,h,'unprobed')] for method in P['probe_methods'] for arm in P['arms'] for h in range(2)])
def distance(a,ka,b,kb):
    lookup={tuple(k):i for i,k in enumerate(kb)};indexes=np.array([lookup[tuple(k)] for k in ka]);wa=np.where(ka[:,2]==0,1.,2.);wb=np.where(kb[:,2]==0,1.,2.)
    numerator=np.sum(abs(a-b[...,indexes])**2*wa,axis=(-2,-1));missing=np.ones(len(kb),bool);missing[indexes]=False
    numerator+=np.sum(abs(b[...,missing])**2*wb[missing],axis=(-2,-1));denominator=np.sum(abs(b)**2*wb,axis=(-2,-1))
    return np.sqrt(numerator/np.maximum(denominator,1e-30))
def comparison(pa,pb):
    a,ma=read(pa,['time','tangent','end_fields','modes']);b,mb=read(pb,['time','tangent','end_fields','modes'])
    ra,rb=response(a,ma),response(b,mb);diff=np.stack([ra[method][arm]-rb[method][arm] for method in P['probe_methods'] for arm in P['arms']])
    field=distance(field_response(a,ma),a['modes'],field_response(b,mb),b['modes'])
    mm,nn=measurements(a,ma),measurements(b,mb)
    return dict(response_ftle_rmse=float(np.sqrt(np.mean(diff**2))),response_field_relative=float(field.max()),
                primary_effect_difference={method:mm[method]['attenuation']-nn[method]['attenuation'] for method in P['probe_methods']})
def effective_regimes():
    cfg=json.loads(json.dumps(P['regimes']))
    amendment=ROOT/'resolution-amendment.json'
    if amendment.exists():
        v=json.loads(amendment.read_text())
        for key,override in v['overrides'].items():cfg[key].update(override)
    return cfg
def filename(seed,regime,layout,n,dt=None,folder='final-data'):
    dt=effective_regimes()[regime]['dt'] if dt is None else dt
    return ROOT/folder/f'{regime}-{layout}-s{seed}-n{n}-dt{dt:g}-gpu.npz'
def screens(m):
    g=P['gates'];maxima=m['maxima'];out={k:maxima[k]<g[k] for k in ['divergence_rms','gradient_trace','energy_budget_relative','tail_enstrophy','tangent_volume_error']}
    out['cfl']=maxima['max_cfl']<g['cfl'];out['work']=m['work_audit']['equal_work_error']<g['work_error'];out['rms']=m['work_audit']['equal_rms_error']<g['rms_error']
    if m['config']['regime']!='smooth':out['kmax_eta']=m['field_stats']['kmax_eta']>=g['kmax_eta']
    return out
def pilot_audit():
    pairs=[];runs=[]
    for regime in P['regimes']:
        cfg=P['regimes'][regime]
        for seed in P['pilot_seeds']:
            lo,hi=[filename(seed,regime,'coincident',n,cfg['dt'],'fast-data') for n in cfg['grids']]
            if not lo.exists() or not hi.exists():continue
            comp=comparison(lo,hi);pair=dict(seed=seed,regime=regime,**comp)
            pair.pop('primary_effect_difference') # Do not inspect/select on pilot effect sizes.
            pair['grid_pass']=comp['response_ftle_rmse']<P['gates']['response_ftle_grid_rmse'] and comp['response_field_relative']<P['gates']['response_field_grid_relative']
            pairs.append(pair)
            for file in [lo,hi]:
                _,m=read(file,[]);runs.append(dict(config=m['config'],screens=screens(m),maxima=m['maxima'],Re_lambda=m['field_stats']['Re_lambda'],kmax_eta=m['field_stats']['kmax_eta']))
    result=dict(pairs=pairs,runs=runs,all_passed=len(pairs)==len(P['regimes'])*len(P['pilot_seeds']) and all(x['grid_pass'] for x in pairs) and all(all(x['screens'].values()) for x in runs),
                effect_sizes_inspected=False)
    dump(ROOT/'results'/'pilot-audit.json',result);return result
def analyse():
    groups={};per=[];grid=[];step=[];numerical=[]
    for regime,cfg in effective_regimes().items():
        seeds=P['smooth_seeds'] if regime=='smooth' else P['forced_seeds']
        for layout in P['layouts']:
            rows=[]
            for seed in seeds:
                pair=[filename(seed,regime,layout,n) for n in cfg['grids']]
                a,m=read(pair[0],['time','tangent']);v=measurements(a,m);rows.append(v);per.append(dict(seed=seed,regime=regime,layout=layout,measurements=v))
                grid.append(dict(seed=seed,regime=regime,layout=layout,**comparison(*pair)))
                for file in pair:
                    _,meta=read(file,[]);numerical.append(dict(config=meta['config'],screens=screens(meta),maxima=meta['maxima'],work_audit=meta['work_audit'],field_stats=meta['field_stats']))
                if seed in seeds[:2]:step.append(dict(seed=seed,regime=regime,layout=layout,**comparison(pair[1],filename(seed,regime,layout,cfg['grids'][-1],cfg['dt']/2))))
            group={method:{metric:interval([r[method][metric] for r in rows]) for metric in ['retained','scrambled','attenuation','specificity']} for method in P['probe_methods']}
            for method in P['probe_methods'][1:]:group[method]['versus_fixed_rms']=interval([r[method]['attenuation']-r['fixed_rms']['attenuation'] for r in rows])
            groups[regime+'/'+layout]=group
    layout_interactions={}
    for regime in P['regimes']:
        base=[v for v in per if v['regime']==regime and v['layout']=='coincident']
        for layout in ['writer_moved','probe_moved']:
            other=[v for v in per if v['regime']==regime and v['layout']==layout]
            layout_interactions[regime+'/'+layout]={method:interval([b['measurements'][method]['attenuation']-a['measurements'][method]['attenuation'] for a,b in zip(base,other)]) for method in P['probe_methods']}
    g=P['gates'];gates=dict(all_trajectory_screens=all(all(m['screens'].values()) for m in numerical),grid_ftle=all(x['response_ftle_rmse']<g['response_ftle_grid_rmse'] for x in grid),grid_field=all(x['response_field_relative']<g['response_field_grid_relative'] for x in grid),half_step=all(x['response_ftle_rmse']<g['response_ftle_step_rmse'] for x in step),reset=all(r['measurements'][method]['reset']<1e-12 for r in per for method in P['probe_methods']))
    result=dict(groups=groups,per_flow=per,layout_interactions=layout_interactions,grid=grid,half_step=step,numerical=numerical,gates=gates,all_passed=all(gates.values()),independent_flow_unit=True,exploratory_comparisons_unadjusted=True)
    dump(ROOT/'results'/'analysis.json',result);return result
if __name__=='__main__':print(json.dumps(analyse()['gates'],indent=2))
