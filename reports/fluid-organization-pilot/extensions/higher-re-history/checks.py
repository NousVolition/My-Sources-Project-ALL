"""Resolution gates fixed before calibration and independent held-out runs."""
import json
import numpy as np
from common import ROOT,read,save,protocol,digest,sources

def response_subset(z,events):
    ids=[int(np.flatnonzero(z['events']==t)[0]) for t in events]
    return {k:np.asarray(z[k])[ids] for k in ['events','response_time','energy_ratio']}

def numerical(path,meta,p):
    with np.load(path) as z:
        rc=z['response_checks'];cfl=float(z['response_cfl'].max())
        finite=all(np.isfinite(z[k]).all() for k in z.files if z[k].dtype.kind in 'fci')
    rows=meta['baseline'];g=p['gates'];c=meta['config'];cut=int(np.ceil(c['n']/3)-1)*2*np.pi/6
    return dict(finite=bool(finite),cfl=max(cfl,max(r['cfl'] for r in rows))<=g['cfl_max'],diffusion=3*c['nu']*cut**2*max(c['spin_dt'],c['response_dt'])<=g['diffusion_max'],total_budget=max(float(rc[...,2].max()),max(abs(r['balance']) for r in rows))<=g['budget_relative_max'],difference_budget=float(rc[...,3].max())<=g['budget_relative_max'],resolution=min(float(rc[...,0].min()),min(r['keta'] for r in rows))>=g['keta_min'],tail=max(float(rc[...,1].max()),max(r['tail_energy'] for r in rows))<=g['tail_energy_max'],divergence=max(r['divergence'] for r in rows)<=g['divergence_max'])

def screening():
    p=protocol();records={};hashes={};arrays={};pilot=[]
    for case,c in enumerate(p['pilot_cases']):
        name=f'data/pilot/pilot_{case}.json';q=read(name)
        assert q['case']==case and q['config']==c
        assert q['source_sha256']=={**sources(),'pilot.py':digest(ROOT/'pilot.py')}
        # Full Fourier snapshots are retained in the pilot artifact. Lightweight
        # calibration/test bundles contain the verified metadata and responses.
        snapshot=ROOT/f'data/pilot/pilot_{case}.npz'
        if snapshot.exists():assert digest(snapshot)==q['sha256']
        hashes[name]=digest(ROOT/name);records[str(case)]=q['checks'];pilot.append(q)
    for label in ['base','half','fine']:
        name=f'data/pilot/response_{label}';q=read(name+'.json');path=ROOT/(name+'.npz');assert digest(path)==q['sha256']
        assert q['source_sha256']=={**sources(),'pilot.py':digest(ROOT/'pilot.py')}
        for key in ['nu','power','horizon','strengths','impulse_width']:assert q['config'][key]==p[key]
        assert q['config']['n']==(p['fine_n'] if label=='fine' else p['n'])
        assert q['config']['response_dt']==p['response_dt']/(1 if label=='base' else 2)
        assert q['input_sha256']==pilot[p['response_pilot']['case']]['sha256']
        records[label]=q['checks'];hashes[name+'.json']=digest(ROOT/(name+'.json'));hashes[name+'.npz']=digest(path)
        with np.load(path) as z:arrays[label]=z['energy_ratio']
    step=float(np.max(abs(arrays['half']-arrays['base']))/np.max(abs(arrays['base'])))
    grid=float(np.max(abs(arrays['fine']-arrays['half']))/np.max(abs(arrays['half'])))
    pairs=[]
    for i in range(0,len(pilot),2):
        a,b=pilot[i:i+2];row={}
        for key in ['energy','dissipation','re_lambda','integral_length','flux_k4']:
            x=float(np.mean([r[key] for r in a['baseline'] if 20<=r['t']<=44]));y=float(np.mean([r[key] for r in b['baseline'] if 20<=r['t']<=44]));row[key]=dict(base=x,fine=y,relative_change=abs(y-x)/max(abs(x),1e-12))
        pairs.append(dict(seed=a['config']['seed'],statistics=row))
    long_ok=all(q['statistics'][k]['relative_change']<=p['gates']['pilot_grid_statistics_relative_max'] for q in pairs for k in ['energy','dissipation','re_lambda','integral_length'])
    # These pilot means choose a regime, not a confidence interval for population
    # stationarity. That separate gate uses all calibration/test runs.
    re=np.mean([r['re_lambda'] for q in pilot for r in q['baseline'] if 20<=r['t']<=44])
    result=dict(checks=records,step_curve_change=step,grid_curve_change=grid,long_grid_pairs=pairs,pilot_mean_re_lambda=float(re),sha256=hashes,passed=bool(all(all(q.values()) for q in records.values()) and long_ok and re>p['gates']['mean_re_lambda_min'] and step<=p['gates']['step_curve_change'] and grid<=p['gates']['grid_curve_change']))
    return result

def verify_all(write=True):
    p=protocol();checks={};refinements=[]
    for split in ['train','validation','test']:
        for seed in p['splits'][split]:
            path=ROOT/f'data/runs/base_{seed}.npz';q=read(f'data/runs/base_{seed}.json');assert digest(path)==q['sha256'];assert q['source_sha256']==sources();checks[path.name]=numerical(path,q,p)
    for seed in p['refinement_seeds']:
        zs={}
        for label in ['base','half','fine']:
            path=ROOT/f'data/runs/{label}_{seed}.npz';q=read(f'data/runs/{label}_{seed}.json');assert digest(path)==q['sha256'];checks[path.name]=numerical(path,q,p)
            with np.load(path) as z:zs[label]={k:z[k] for k in z.files}
        for label,reference in [('half','base'),('fine','half')]:
            subset=response_subset(zs[reference],zs[label]['events'])
            a=zs[label]['energy_ratio'];b=subset['energy_ratio']
            relative=float(np.max(abs(a-b))/max(float(np.max(abs(b))),1e-12))
            from features import target
            ya=target(zs[label]);yb=target(subset)
            refinements.append(dict(seed=seed,label=label,reference=reference,events=zs[label]['events'].tolist(),max_curve_relative_change=relative,max_target_relative_change=float(np.max(abs(ya-yb)/abs(yb))),passed=relative<=p['gates']['step_curve_change' if label=='half' else 'grid_curve_change']))
    screen=screening();result=dict(numerical=checks,refinements=refinements,screening=screen,all_passed=bool(screen['passed'] and all(all(q.values()) for q in checks.values()) and all(q['passed'] for q in refinements)))
    if write:save('results/numerical_checks.json',result)
    print(json.dumps({'all_passed':result['all_passed'],'runs_checked':len(checks),'refinement_pairs':len(refinements)}));return result

if __name__=='__main__':verify_all()
