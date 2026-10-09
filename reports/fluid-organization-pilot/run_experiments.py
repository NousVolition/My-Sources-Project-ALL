"""Run gates, then the planned ensemble. Cached results require matching hashes."""
import argparse
import json
from pathlib import Path
from core import ROOT,run,save_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['screen','ensemble','regimes','refine','all'],default='screen')
    ap.add_argument('--data',type=Path,default=ROOT/'data');a=ap.parse_args()
    p=json.loads((ROOT/'protocol.json').read_text());base=p['base'];d=a.data
    if a.stage in ['screen','all']:
        paths=[]
        for name,change in [('screen',{}),('screen_dt',{'dt':.005}),('screen_grid',{'n':38})]:
            paths.append(run(6000,{**base,**change},d,name))
        import numpy as np
        meta=[json.loads(x.with_suffix('.json').read_text()) for x in paths]
        gains=[r['diagnostics'][-1]['perturbation_energy']/r['diagnostics'][0]['perturbation_energy'] for r in meta]
        budget=max(abs(row[k]) for m in meta for row in m['diagnostics'] for k in ['baseline_energy_residual','disturbed_energy_residual','perturbation_budget_residual'])
        arrays=[np.load(x) for x in paths]
        marker=float(np.sqrt(np.mean(np.sum((arrays[2]['positions_perturbed'][-1]-arrays[0]['positions_perturbed'][-1])**2,axis=1))))
        checks={'max_budget_residual':budget,'dt_gain_relative_change':abs(gains[1]/gains[0]-1),'grid_gain_relative_change':abs(gains[2]/gains[0]-1),'grid_marker_rms':marker}
        checks['passed']=budget<1e-4 and checks['dt_gain_relative_change']<.005 and checks['grid_gain_relative_change']<.02 and marker<.005
        save_json(d/'screen_checks.json',checks);print(checks,flush=True)
        if not checks['passed']:raise RuntimeError('Screen failed; do not expand')
    if a.stage!='screen':
        if not (d/'screen_checks.json').exists() or not json.loads((d/'screen_checks.json').read_text())['passed']:
            raise RuntimeError('Run and pass --stage screen first')
    if a.stage in ['ensemble','all']:
        for seed in p['train_seeds']+p['validation_seeds']+p['test_seeds']:run(seed,base,d/'ensemble')
    if a.stage in ['regimes','all']:
        for label,change in p['regimes'].items():
            for seed in p['regime_seeds']:run(seed,{**base,**change},d/'regimes',f'{label}_{seed}')
    if a.stage in ['refine','all']:
        for seed in p['refinement_seeds']:
            for name,change in [('grid',{'n':38}),('dt',{'dt':.005})]:
                run(seed,{**base,**change},d/'refinement',f'{name}_{seed}')
        for name,change in [('fast_water',{'nu':.01,'speed':2.}),('water',{'nu':.01}),('air',{'nu':.15})]:
            for suffix,tweak in [('grid',{'n':38,'dt':.005}),('dt',{'dt':.005})]:
                run(7000,{**base,**change,**tweak},d/'refinement',f'{name}_{suffix}_7000')


if __name__=='__main__':main()
