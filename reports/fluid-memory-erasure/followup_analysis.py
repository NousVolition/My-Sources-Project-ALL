"""Secondary erasure/persistence audits; no fitted models or tuned thresholds."""
import json
from pathlib import Path
import numpy as np
from run import ROOT,P,dump,digest,filename
from analyze import BACKEND,N_MAIN,load_run,ftle,causal_metrics,interval,numerical_comparison,spectral_distance,matched_sham_study


def read_aux(path):
    m=json.loads(path.with_suffix('.json').read_text())
    if digest(path)!=m['raw_sha256']:raise RuntimeError('Auxiliary raw hash mismatch')
    for f,sha in m['sources'].items():
        if digest(ROOT/f)!=sha:raise RuntimeError('Auxiliary source changed: '+f)
    with np.load(path,allow_pickle=False) as z:a={k:z[k] for k in z.files}
    if not all(np.isfinite(v).all() for v in a.values()):raise RuntimeError('Nonfinite auxiliary array')
    return a,m


def gates(maxima):
    p=P['numerical_gates']
    mapping={'divergence_rms':'divergence_rms','gradient_trace':'marker_gradient_trace','energy_budget_relative':'relative_energy_budget',
             'max_cfl':'cfl','tail_enstrophy':'tail_enstrophy_fraction','tangent_volume_error':'tangent_determinant_error'}
    return {k:bool(maxima[k]<p[v]) for k,v in mapping.items()}


def main():
    matched=matched_sham_study(ROOT/'data',ROOT/'matched-sham-data')
    dump(ROOT/'results/matched-sham-results.json',matched)
    auxiliary=[];maxima={};matched_time=[];delay_grid=[];delay_time=[];delayed=[]
    for folder in ['matched-sham-data','delayed-probe-data']:
        for path in sorted((ROOT/folder).glob('s*.npz')):
            a,m=read_aux(path)
            for k,v in m['maxima'].items():maxima[k]=max(maxima.get(k,0.),abs(v))
            if folder=='matched-sham-data':
                if max(m['matching']['absolute_matching_error'],m['per_mode_power_error'],m['coarse_coefficient_error'])>1e-10:raise RuntimeError('Matching/invariant error')
            else:
                for t in P['horizons']:delayed.append(dict(config=m['config'],horizon=t,**causal_metrics(a,m,t)))
                if m['config']['wait']==.8:
                    c=m['config'];source,_=load_run(ROOT/'data',c['seed'],n=c['n'],dt=c['dt'])
                    if not np.array_equal(a['waited_fields'],source['end_fields'][::2]):raise RuntimeError('Reused checkpoint not exact')
            auxiliary.append(dict(file=str(path.relative_to(ROOT)),sha256=m['raw_sha256'],config=m['config']))
    matching=json.loads((ROOT/'matched-sham-protocol.json').read_text())
    for seed in matching['half_step_seeds']:
        pairs=[]
        for dt in [P['dt'],P['dt']/2]:
            stem=filename(seed,matching['half_step_n'],dt,P['nu'],P['probe_velocity_rms'],.8,BACKEND)+'-matched-sham.npz'
            pairs.append(read_aux(ROOT/'matched-sham-data'/stem)[0])
        f=[ftle(a,P['primary_horizon']) for a in pairs];r=[np.stack([v[1]-v[0],v[3]-v[2]]) for v in f]
        matched_time.append(dict(seed=seed,response_ftle_rmse=float(np.sqrt(np.mean((r[0]-r[1])**2)))))
    design=json.loads((ROOT/'delayed-probe-protocol.json').read_text())
    def delay(seed,n,dt):
        stem=filename(seed,n,dt,P['nu'],P['probe_velocity_rms'],.8,BACKEND)+'-wait0.8.npz'
        return read_aux(ROOT/'delayed-probe-data'/stem)
    for seed in design['fine_seeds']:
        delay_grid.append(dict(seed=seed,**numerical_comparison(delay(seed,N_MAIN,P['dt']),delay(seed,design['fine_n'],P['dt']),P['primary_horizon'])))
    for seed in design['half_step_seeds']:
        delay_time.append(dict(seed=seed,**numerical_comparison(delay(seed,N_MAIN,P['dt']),delay(seed,N_MAIN,P['dt']/2),P['primary_horizon'])))
    summary={}
    for wait in design['waits']:
        summary[str(wait)]={}
        for t in P['horizons']:
            rows=[r for r in delayed if r['config']['n']==N_MAIN and r['config']['dt']==P['dt'] and r['config']['wait']==wait and r['horizon']==t]
            summary[str(wait)][str(t)]={k:interval([r[k] for r in rows]) for k in ['primary_attenuation','retained_history_contrast','scrambled_history_contrast']}
    p=P['numerical_gates'];passes=gates(maxima)
    passes.update(matched_grid=all(r['response_ftle_rmse']<p['finest_pair_response_ftle_rmse'] for r in matched['numerical']),
                  matched_half_step=all(r['response_ftle_rmse']<p['time_response_ftle_rmse'] for r in matched_time),
                  delayed_grid=all(r['response_ftle_rmse']<p['finest_pair_response_ftle_rmse'] and r['ftle_rmse']<p['finest_pair_ftle_rmse'] and r['max_probe_response_field_relative_L2']<p['finest_pair_probe_response_field_relative_L2'] for r in delay_grid),
                  delayed_half_step=all(r['response_ftle_rmse']<p['time_response_ftle_rmse'] for r in delay_time))
    dump(ROOT/'results/followup-numerical-audit.json',dict(gates=passes,all_passed=all(passes.values()),maxima=maxima,matched_half_step=matched_time,
         delayed_grid=delay_grid,delayed_half_step=delay_time,verified_auxiliary_recordings=auxiliary))
    dump(ROOT/'results/delayed-probe-results.json',dict(horizons=summary,per_flow=delayed,
         interpretation='Four-flow exploratory waiting experiment. A later probe is distinct from a longer response observation. Two-flow grid/time checks do not verify every delayed branch. Matching at erasure need not persist during the wait.'))
    print(json.dumps(dict(auxiliary_recordings=len(auxiliary),numerical_gates=passes),indent=2))


if __name__=='__main__':main()
