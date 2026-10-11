"""Select feasible production settings using numerical errors, not effects."""
import json,time
from production_experiment import ROOT,P,run_one,dump,digest
from analysis import read,comparison,screens,filename
from independent_solver import IndependentSolver

def qualify():
    audit=json.loads((ROOT/'results'/'pilot-audit.json').read_text())
    if len(audit['pairs'])!=6:raise RuntimeError('Original pilot matrix incomplete')
    reference_screens=[r for r in audit['runs'] if r['config']['n']==max(P['regimes'][r['config']['regime']]['grids'])]
    if not all(all(r['screens'].values()) for r in reference_screens):raise RuntimeError('Finest pilot trajectory screens failed; further refinement required')
    pairs=[];rows=[];equivalence=[]
    for regime,n in [('forced_nu004',48),('forced_nu002',80)]:
        for seed in P['pilot_seeds']:
            run_one(seed,regime,'coincident',n)
            trial=filename(seed,regime,'coincident',n,P['regimes'][regime]['dt'])
            reference=filename(seed,regime,'coincident',64 if n==48 else 96,P['regimes'][regime]['dt'],'fast-data')
            comp=comparison(trial,reference);comp.pop('primary_effect_difference')
            _,meta=read(trial,[]);s=screens(meta)
            passed=all(s.values()) and comp['response_ftle_rmse']<P['gates']['response_ftle_grid_rmse'] and comp['response_field_relative']<P['gates']['response_field_grid_relative']
            pairs.append(dict(seed=seed,regime=regime,n=n,comparison=comp,screens=s,passed=passed));print('qualification '+json.dumps(pairs[-1]),flush=True)
    # Verify changed batch scheduling against an already completed case.
    regime='forced_nu004';seed=P['pilot_seeds'][0];run_one(seed,regime,'coincident',64)
    a=filename(seed,regime,'coincident',64,P['regimes'][regime]['dt']);b=filename(seed,regime,'coincident',64,P['regimes'][regime]['dt'],'fast-data')
    comp=comparison(a,b);comp.pop('primary_effect_difference');equivalence.append(dict(comparison=comp,passed=comp['response_ftle_rmse']<1e-10 and comp['response_field_relative']<1e-10))
    result=dict(all_passed=all(r['passed'] for r in pairs+equivalence),pairs=pairs,batch_equivalence=equivalence,effect_sizes_inspected=False,original_pilot_all_passed=audit['all_passed'],reference_screens_passed=True)
    dump(ROOT/'results'/'production-qualification.json',result)
    if not result['all_passed']:raise RuntimeError('Smaller grids/batching not qualified; do not publish a validated effect')
    dump(ROOT/'resolution-amendment.json',dict(date='2026-10-11',before_production_effects=True,overrides={'smooth':{'dt':.005},'forced_nu004':{'grids':[48,64]},'forced_nu002':{'grids':[80,96]}},
         reason='Runtime/memory feasibility and a stronger margin at lower viscosity: select 48/64 for nu=0.04 after direct numerical qualification, and promote the nu=0.02 baseline to 80 with 96 refinement because its first 64/96 pilot is close to the FTLE error limit. No tolerances changed and no causal effect sizes inspected. Halve smooth step for CFL margin. Shared spinup remains on n=96; first-two forced coincident flows at nu=0.04 also have extra 96 and common-probe audits.',
         original_protocol_sha256=digest(ROOT/'protocol.json'),qualification_sha256=digest(ROOT/'results'/'production-qualification.json'),extra_finest_audits=2))
    return result

if __name__=='__main__':qualify()
