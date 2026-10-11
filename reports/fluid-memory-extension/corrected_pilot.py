"""Requalify exact integer-wave implementation; no history effects inspected."""
import json
from final_experiment import ROOT,P,run_one,dump,digest
from analysis import read,filename,comparison,screens

def qualify():
    pairs=[]
    # Start with the regime that failed; never relax its thresholds.
    designs=[('forced_nu002',[80,96]),('smooth',[32,48]),('forced_nu004',[48,64])]
    for regime,grids in designs:
        for seed in P['pilot_seeds']:
            for n in grids:run_one(seed,regime,'coincident',n,.005 if regime=='smooth' else .002)
            dt=.005 if regime=='smooth' else .002
            files=[filename(seed,regime,'coincident',n,dt) for n in grids]
            comp=comparison(*files);comp.pop('primary_effect_difference')
            ss=[screens(read(f,[])[1]) for f in files]
            passed=all(all(s.values()) for s in ss) and comp['response_ftle_rmse']<P['gates']['response_ftle_grid_rmse'] and comp['response_field_relative']<P['gates']['response_field_grid_relative']
            pairs.append(dict(seed=seed,regime=regime,grids=grids,comparison=comp,screens=ss,passed=passed))
            result=dict(all_passed=all(r['passed'] for r in pairs),complete=len(pairs)==6,pairs=pairs,effect_sizes_inspected=False,sources={f:digest(ROOT/f) for f in ['corrected_pilot.py','corrected_solver.py','corrected_observation.py','final_experiment.py']})
            dump(ROOT/'results'/'corrected-pilot-qualification.json',result)
            print('corrected pilot '+json.dumps(pairs[-1]),flush=True)
    if not result['all_passed']:raise RuntimeError('Corrected pilots still fail; physical interpretations withheld')
    dump(ROOT/'resolution-amendment.json',dict(date='2026-10-11',before_production_effects=True,overrides={'smooth':{'dt':.005},'forced_nu004':{'grids':[48,64]},'forced_nu002':{'grids':[80,96]}},
         reason='New corrected integer-wave solver and spinup. All six two-grid pilots pass the unchanged gates after fixing the grid-dependent reset degeneracy; use 48/64 moderate and 80/96 lower viscosity. No pilot effect sizes inspected. Extra n=96 moderate audits retained.',
         correction_sha256=digest(ROOT/'integer-wave-correction.json'),qualification_sha256=digest(ROOT/'results'/'corrected-pilot-qualification.json')))
    return result
if __name__=='__main__':qualify()
