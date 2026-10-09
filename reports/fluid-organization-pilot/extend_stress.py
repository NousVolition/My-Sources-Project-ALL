"""Targeted continuation after the 26-to-38 fast-water marker screen failed.

These are dependent sensitivity runs (same seed 7000), not new statistical runs.
38->50->62 isolates spatial refinement with dt=.005 fixed throughout.
"""
import json
import numpy as np
from core import ROOT,run,wrap,save_json


def main():
    p=json.loads((ROOT/'protocol.json').read_text());c={**p['base'],'nu':.01,'speed':2.,'dt':.005}
    out=ROOT/'data/stress_extension'
    for n in [50,62]:run(7000,{**c,'n':n},out,f'fast_water_n{n}_7000')
    paths=[ROOT/'data/refinement/fast_water_grid_7000.npz']+[out/f'fast_water_n{n}_7000.npz' for n in [50,62]]
    rec=[]
    for a,b in zip(paths[:-1],paths[1:]):
        za=np.load(a);zb=np.load(b);ma=json.loads(a.with_suffix('.json').read_text());mb=json.loads(b.with_suffix('.json').read_text())
        gain=lambda m:m['diagnostics'][-1]['perturbation_energy']/m['diagnostics'][0]['perturbation_energy']
        change=float(np.sqrt(np.mean(np.sum(wrap(zb['positions_perturbed'][-1]-za['positions_perturbed'][-1])**2,axis=1))))
        rec.append({'from_n':ma['config']['n'],'to_n':mb['config']['n'],'dt':.005,'marker_rms_change':change,
                    'gain_relative_change':abs(gain(mb)/gain(ma)-1),'fine_gain':gain(mb),'marker_screen_passed':change<.005})
    save_json(ROOT/'results/stress_extension.json',{'seed':7000,'comparisons':rec,'interpretation':'One-seed spatial sensitivity only; no evidence that every stress realization is resolved.'})
    print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
