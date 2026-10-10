"""Recheck the three historical summary tables without evolving a fluid field."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def audit(folder):
    names = ['bound-check.json', 'ratio-rose.json', 'pressure-at-peak.json']
    data = {n: json.loads((folder / 'results' / n).read_text()) for n in names}
    bound = data['bound-check.json']['n48']
    rows = bound['series']
    times = [r['t'] for r in rows]
    w = [r['vortMax'] for r in rows]
    assert all(b > a for a, b in zip(times, times[1:]))
    intervals = [dict(t0=times[i], t1=times[i+1], W0=w[i], W1=w[i+1],
                      trapezoid=(times[i+1]-times[i])*(w[i+1]+w[i])/2)
                 for i in range(len(rows)-1)]
    integral = math.fsum(r['trapezoid'] for r in intervals)
    # Trapezoidal integral is linear in its time coordinates when W is fixed.
    # Include endpoint errors too, yielding a conservative conditional bound.
    coeff = [-(w[0]+w[1])/2] + [(w[i-1]-w[i+1])/2 for i in range(1,len(w)-1)] + [(w[-2]+w[-1])/2]
    rounding_bound = 0.00005 * math.fsum(abs(c) for c in coeff)
    ratio = data['ratio-rose.json']['series']
    high = max(ratio, key=lambda r:r['ratio'])
    first = ratio[0]
    assert all(math.isfinite(r[k]) for r in ratio for k in ['t','biggest','average','ratio'])
    residual = max(abs(r['ratio'] - r['biggest']/r['average']) for r in ratio)
    assert residual < 1e-12
    pressure = data['pressure-at-peak.json']
    repair = json.loads((folder/'corrected_v1/verification.json').read_text())
    hashes = {n: hashlib.sha256((folder/'corrected_v1'/n).read_bytes()).hexdigest()
              for n in repair['source_sha256']}
    assert hashes == repair['source_sha256']
    return {
        'scope':'Historical table arithmetic and current repair provenance; zero new fluid trajectories',
        'source_sha256':{n:hashlib.sha256((folder/'results'/n).read_bytes()).hexdigest() for n in names},
        'repaired_source_hashes_match_saved_verification':True,
        'pressure':{'top_level_keys':list(pressure), 'historical_split_reconstructed':False,
                    'needed':'Original velocity fields at the stated times, viscosity, domain/grid, producing solver/filter and exact time-step record; or the complete original replay package. A sign flip alone is insufficient.'},
        'integral':{'saved_row_trapezoid':integral, 'reported':bound['integral'],
                    'difference':integral-bound['integral'],
                    'difference_percent_of_reported':100*(integral/bound['integral']-1),
                    'maximum_time_rounding_effect_if_nearest_4_decimals_and_fixed_W':rounding_bound,
                    'rounding_alone_explains_difference':abs(integral-bound['integral'])<=rounding_bound,
                    'intervals':intervals,
                    'historical_internal_step_integral_reconstructed':False,
                    'needed':'Quadrature formula plus every unrounded time and maximum-vorticity value used in that quadrature. The saved-row trapezoid is a different calculation and does not replace the original value.'},
        'ratio':{'maximum_saved_ratio_row':high, 'maximum_saved_vorticity_row':max(ratio,key=lambda r:r['biggest']),
                 'arithmetic_max_absolute_residual':residual,
                 'at_ratio_maximum_change_from_initial_percent':{k:100*(high[k]/first[k]-1) for k in ['biggest','average','ratio']},
                 'trajectory_reconstructed':False,
                 'needed':'Grid, domain, viscosity, complete initial field or construction and parameters, force, solver/filter revision, exact step schedule, and average/norm definition.'},
    }


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    result=audit(args.folder)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['scope','repaired_source_hashes_match_saved_verification']}))
