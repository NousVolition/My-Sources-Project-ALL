import argparse,json
import numpy as np
from common import ROOT,save,digest,read_protocol,numerical_checks
def pair(a,b):
    with np.load(a) as x,np.load(b) as y:
        ex=x['perturbation_energy']/x['initial_perturbation_energy'];ey=y['perturbation_energy']/y['initial_perturbation_energy']
        assert np.allclose(x['time'],y['time']);assert np.allclose(x['amplitudes'],y['amplitudes'],rtol=.002)
        return {'energy_curve_relative':float(np.max(abs(ex-ey))/max(1,float(abs(ex).max()))),'marker_position_rms':float(np.max(np.sqrt(np.mean(np.sum((x['positions']-y['positions'])**2,axis=-1),axis=(1,2))))),'target_relative':float(np.max(abs(np.trapezoid(ex,x['time'],axis=0)-np.trapezoid(ey,y['time'],axis=0))/np.trapezoid(ex,x['time'],axis=0)))}
def verify(screen=False):
    p=read_protocol(screen);folder=ROOT/'data'/('screen' if screen else 'runs');rows=[];comparisons=[]
    sources={s:digest(ROOT/s) for s in ['simulate.py','physics.py','vendor/numerics.py','common.py']}
    for path in sorted(folder.glob('*.npz')):
        m=json.loads(path.with_suffix('.json').read_text(encoding='utf-8'));assert digest(path)==m['sha256']
        expected=sources.copy()
        if path.stem=='smooth_base_90000' and screen:expected['simulate.py']=digest(ROOT/'screen_initial_simulate.py')
        elif screen and 'cubic' in path.stem:expected['simulate.py']=digest(ROOT/'screen_cubic_simulate.py')
        elif screen and 'cubic' not in path.stem:expected['simulate.py']=digest(ROOT/'screen_trilinear_simulate.py')
        # The audit snapshot used the original common helper before gate-check support.
        if screen and m['sources']['common.py']==digest(ROOT/'screen_initial_common.py'):
            expected['common.py']=digest(ROOT/'screen_initial_common.py')
        assert m['sources']==expected,path.name
        checks=numerical_checks(m,p['gates']);rows.append({'file':path.name,'sha256':digest(path),'checks':checks,**{k:m[k] for k in ['max_cfl','diffusion_number','min_kcut_eta','max_high_fraction','max_energy_balance','max_difference_balance','max_divergence']}})
    for regime in p['regimes']:
        seeds=p['marker_revalidation_seeds'] if screen else p['refinement_seeds']
        for seed in seeds:
            base=folder/f'{regime}_{"cubic" if screen else "base"}_{seed}.npz'
            for label in ['fine','half']:
                other=folder/f'{regime}_{"cubic_" if screen else ""}{label}_{seed}.npz'
                if not other.exists():
                    if label=='half' and screen and seed==90001:continue
                    raise FileNotFoundError(other)
                q=pair(base,other);q.update(regime=regime,seed=seed,label=label)
                q['passed']=q['energy_curve_relative']<=p['gates']['grid_energy_curve_relative' if label=='fine' else 'step_energy_curve_relative'] and q['marker_position_rms']<=p['gates']['grid_marker_position_rms'];comparisons.append(q)
    result={'screen':screen,'runs':rows,'refinements':comparisons,'all_numerical_passed':all(all(q['checks'].values()) for q in rows),'all_refinements_passed':all(q['passed'] for q in comparisons)}
    if screen:
        original=json.loads((ROOT/'results/screen_initial_verification.json').read_text(encoding='utf-8'));result['historical_failed_marker_screen']=original['refinements'];result['historical_energy_refinements_passed']=all(q['energy_curve_relative']<=p['gates']['grid_energy_curve_relative' if q['label']=='fine' else 'step_energy_curve_relative'] for q in original['refinements'])
        for regime in p['regimes']:
            for seed in p['marker_revalidation_seeds']:
                with np.load(folder/f'{regime}_cubic_{seed}.npz') as a,np.load(folder/f'{regime}_revised_{seed}.npz') as b:
                    assert np.array_equal(a['energy'],b['energy']) and np.array_equal(a['perturbation_energy'],b['perturbation_energy']),'Marker-only revision altered the fluid calculation'
        result['marker_revision_leaves_fluid_energy_bit_identical']=True
    save('results/'+('screen' if screen else 'numerical')+'_verification.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='runs'},indent=2))
    if not result['all_numerical_passed'] or not result['all_refinements_passed']:raise ValueError('Predeclared numerical gate failed')
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--screen',action='store_true');verify(ap.parse_args().screen)
