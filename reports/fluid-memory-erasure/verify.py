"""Recompute scientific summaries and check every recording without rerunning NS."""
import argparse,json,sys,unittest
from pathlib import Path
import numpy as np
from run import ROOT,P,digest,dump
from analyze import load_run,causal_metrics,branch_index,features,ftle,N_MAIN,numerical_audit


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');ap.add_argument('--out',type=Path,default=ROOT/'results');args=ap.parse_args()
    errors=[];checked=0;checks=[]
    for mp in args.data.glob('s*.json'):
        meta=json.loads(mp.read_text());cfg=meta['config'];z,m=load_run(args.data,**cfg)
        if not all(np.isfinite(value).all() for value in z.values()):errors.append(mp.name+': nonfinite')
        for name,c in m['intervention_checks'].items():
            if max(c.values())>1e-10:errors.append(mp.name+': spectrum/coarse/divergence intervention failure '+name)
        # Recompute current marker-gradient trace and determinant directly from data.
        if np.max(abs(np.trace(z['prep_gradient'],axis1=-2,axis2=-1)))>1e-10:errors.append(mp.name+': trace')
        if abs(np.linalg.det(z['tangent'])-1).max()>P['numerical_gates']['tangent_determinant_error']:errors.append(mp.name+': volume')
        for kicked in [0,1]:
            a,b=[branch_index(m,'reset',history,kicked) for history in [0,1]]
            if max(abs(z['positions'][:,a]-z['positions'][:,b]).max(),abs(z['tangent'][:,a]-z['tangent'][:,b]).max())>1e-11:
                errors.append(mp.name+': common reset identity')
        if cfg['n']==N_MAIN and cfg['seed'] in P['test_seeds'] and cfg['dt']==P['dt'] and cfg['nu']==P['nu'] and cfg['amplitude']==P['probe_velocity_rms'] and cfg['end']==.8:
            order=np.array([7,2,9,0,1,3,4,5,6,8]);indexes=np.concatenate([order,np.concatenate([np.arange(10+6*i,10+6*i+6) for i in order])])
            error=0.
            for history in [0,1]:
                a,h=features(z['prep_time'],z['prep_positions'][:,history],z['prep_velocity'][-1,history],z['prep_gradient'][-1,history])
                ar,hr=features(z['prep_time'],z['prep_positions'][:,history,indexes],z['prep_velocity'][-1,history,indexes],z['prep_gradient'][-1,history,order])
                error=max(error,float(abs(ar-a[order]).max()),float(abs(hr-h[order]).max()))
            relabeled=dict(z);relabeled['tangent']=z['tangent'][:,:,order]
            original=causal_metrics(z,m,P['primary_horizon']);permuted=causal_metrics(relabeled,m,P['primary_horizon'])
            target_error=max(abs(original[k]-permuted[k]) for k in ['primary_attenuation','retained_history_contrast','scrambled_history_contrast','sham_attenuation','retained_signed_response_mean'])
            checks.append(dict(seed=cfg['seed'],actual_history_relabeling_error=error,actual_response_target_relabeling_error=target_error))
            if error>1e-10:errors.append(mp.name+': label invariance')
            if target_error>1e-10:errors.append(mp.name+': target label invariance')
        checked+=1
    for folder in [ROOT/'pilot-data',ROOT/'matched-sham-data',ROOT/'delayed-probe-data']:
        for mp in folder.glob('s*.json'):
            m=json.loads(mp.read_text());rp=mp.with_suffix('.npz')
            if digest(rp)!=m['raw_sha256']:errors.append(str(rp)+': hash')
    duration_checks=[]
    for seed in P['long_horizon_seeds']:
        short,_=load_run(args.data,seed,keys=['positions','tangent','response_fields'])
        long,_=load_run(args.data,seed,end=max(P['long_horizons']),keys=['positions','tangent','response_fields'])
        error=max(float(abs(short[k]-long[k][:len(short[k])]).max()) for k in short)
        duration_checks.append(dict(seed=seed,prefix_error=error))
        if error>1e-11:errors.append('Long-run prefix differs: '+str(seed))
    freeze=json.loads((args.out/'model-freeze-provenance.json').read_text())
    if digest(args.out/'frozen-models.json')!=freeze['models_sha256']:errors.append('Frozen model changed')
    for f,sha in freeze['source_sha256'].items():
        candidates=[ROOT/'freeze-source'/f,ROOT/f]
        if not any(p.exists() and digest(p)==sha for p in candidates):errors.append('Missing exact fitting source: '+f)
    suite=unittest.defaultTestLoader.loadTestsFromNames(['test_solver','test_analysis','test_matched_sham'])
    runner=unittest.TextTestRunner(verbosity=1);result=runner.run(suite)
    summary=dict(all_passed=not errors and result.wasSuccessful(),errors=errors,production_configurations_verified=checked,
                 analytic_and_analysis_tests=result.testsRun,failed_tests=len(result.failures)+len(result.errors),actual_label_controls=checks,
                 duration_prefix_controls=duration_checks,source_sha256={p.name:digest(p) for p in ROOT.glob('*.py')},protocol_sha256=digest(ROOT/'protocol.json'),
                 note='Integrity, invariance and numerical screens are distinct from statistical evidence. Verification does not establish continuum convergence globally.')
    args.out.mkdir(parents=True,exist_ok=True);dump(args.out/'verification.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','actual_label_controls']},indent=2))
    if not summary['all_passed']:raise SystemExit(1)

if __name__=='__main__':main()
