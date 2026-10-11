"""Frozen production matrix; prerequisites prevent release on failed pilot gates."""
import argparse,json
from final_experiment import ROOT,P,run_one,dump,digest
from analysis import effective_regimes
from common_probe import run_one as common
from unforced_readout import run_one as unforced
from common_unforced import run_one as common_off

def tasks():
    main=[];extra=[];common_tasks=[];off=[]
    for regime,cfg in effective_regimes().items():
        seeds=P['smooth_seeds'] if regime=='smooth' else P['forced_seeds']
        for seed in seeds:
            for layout in P['layouts']:
                for n in cfg['grids']:main.append(dict(seed=seed,regime=regime,layout=layout,n=n,dt=cfg['dt']))
                if seed in seeds[:2]:main.append(dict(seed=seed,regime=regime,layout=layout,n=cfg['grids'][-1],dt=cfg['dt']/2))
                if regime=='smooth' or layout=='coincident':
                    for n in cfg['grids']:common_tasks.append(dict(seed=seed,regime=regime,layout=layout,n=n,dt=cfg['dt']))
                    if seed in seeds[:2] and (regime=='smooth' or seed==seeds[0]):common_tasks.append(dict(seed=seed,regime=regime,layout=layout,n=cfg['grids'][-1],dt=cfg['dt']/2))
                if regime!='smooth' and layout=='coincident' and seed in seeds[:2]:
                    for n in cfg['grids']:off.append(dict(seed=seed,regime=regime,n=n,dt=cfg['dt']))
                    if seed==seeds[0]:off.append(dict(seed=seed,regime=regime,n=cfg['grids'][-1],dt=cfg['dt']/2))
    for seed in P['smooth_seeds'][:2]:extra.append(dict(seed=seed,regime='smooth',layout='coincident',n=32,dt=effective_regimes()['smooth']['dt'],translation=[.37,-.29,.61]))
    finest=[dict(seed=seed,regime=regime,layout='coincident',n=96,dt=effective_regimes()[regime]['dt']) for regime in ['forced_nu004','forced_nu002'] if max(effective_regimes()[regime]['grids'])<96 for seed in P['forced_seeds'][:2]]
    return dict(main=main,resolution_audit=finest,translation=extra,common_probe=common_tasks,common_resolution=finest,unforced_readout=off,common_unforced_readout=off)

def execute(stage):
    plan=tasks();path=ROOT/'production-plan.amended.json'
    if path.exists():
        if json.loads(path.read_text())['tasks']!=plan:raise RuntimeError('Frozen plan mismatch')
    else:dump(path,dict(protocol_sha256=digest(ROOT/'protocol.json'),tasks=plan,counts={k:len(v) for k,v in plan.items()},production_effects_inspected=False))
    audit=json.loads((ROOT/'results'/'pilot-audit.json').read_text())
    fast=json.loads((ROOT/'results'/'fast-observation-validation.json').read_text())
    equivalence=json.loads((ROOT/'results'/'accelerated-pilot-equivalence.json').read_text())
    qualification=json.loads((ROOT/'results'/'corrected-pilot-qualification.json').read_text())
    if not(qualification['all_passed'] and qualification['complete'] and json.loads((ROOT/'results'/'corrected-validation.json').read_text())['all_passed']):raise RuntimeError('Production withheld: corrected pilot settings not qualified')
    completed=[];method={'main':run_one,'resolution_audit':run_one,'translation':run_one,'common_probe':common,'common_resolution':common,'unforced_readout':unforced,'common_unforced_readout':common_off}[stage]
    for task in plan[stage]:
        result=method(**task);completed.append(result['config'])
        dump(ROOT/'results'/('execution-'+stage+'.json'),dict(stage=stage,planned=len(plan[stage]),completed=len(completed),configs=completed,production_plan_sha256=digest(path),controller_sha256=digest(ROOT/'production.py')))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['main','resolution_audit','translation','common_probe','common_resolution','unforced_readout','common_unforced_readout'],default='main');a=ap.parse_args();execute(a.stage)
