"""Check supplied JSON arithmetic and optionally rebuild only the initial fields.

This does not rerun a trajectory or certify the reported long-run measurements.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import solver

HERE=Path(__file__).resolve().parent


def numeric_values(obj):
    if isinstance(obj,dict):
        for value in obj.values():
            yield from numeric_values(value)
    elif isinstance(obj,list):
        for value in obj:
            yield from numeric_values(value)
    elif isinstance(obj,(float,int)):
        yield obj


def series_summary(rows):
    assert len(rows)>1 and rows[0]['t']==0
    assert all(b['t']>a['t'] for a,b in zip(rows,rows[1:]))
    out={'samples':len(rows),'end_time':rows[-1]['t']}
    if 'peak' in rows[0]:
        out.update(start_peak=rows[0]['peak'],end_peak=rows[-1]['peak'],
                   largest_saved_peak=max(r['peak'] for r in rows))
    return out


def check_results(rebuild_starts=False):
    received={p.name:json.loads(p.read_text(encoding='utf-8')) for p in sorted((HERE/'results').glob('*.json'))}
    manifest=json.loads((HERE/'provenance.json').read_text(encoding='utf-8'))
    expected={Path(r['published_path']).name for r in manifest['source_files']+manifest.get('supplemental_files',[])
              if r['published_path'].startswith('results/')}
    assert set(received)==expected
    assert all(math.isfinite(v) for data in received.values() for v in numeric_values(data))
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((HERE/'results').glob('*.json'))}
    series={}
    for name in ('tube-finer.json','two-tubes.json'):
        data=received[name]
        rows=data['series']
        series[name]=series_summary(rows)
        assert math.isclose(data['startPeak'],rows[0]['peak'],rel_tol=1e-12)
        assert data['highPeak']>=max(r['peak'] for r in rows)-1e-12
        series[name]['energy_decreases_at_saved_samples']=all(b['energy']<=a['energy'] for a,b in zip(rows,rows[1:]))
    for name in ('spin-ratio.json','spin-ratio-64.json','spin-ratio-low-nu.json','ratio-rose.json'):
        data=received[name]
        rows=data if isinstance(data,list) else data['series']
        series[name]=series_summary(rows)
        assert all(r['average']>0 and math.isclose(r['ratio'],r['biggest']/r['average'],rel_tol=1e-12) for r in rows)
        series[name].update(start_ratio=rows[0]['ratio'],end_ratio=rows[-1]['ratio'],
                            monotone_decrease=all(b['ratio']<=a['ratio'] for a,b in zip(rows,rows[1:])))
        if name=='ratio-rose.json':
            series[name]['highest_ratio_row']=max(rows,key=lambda r:r['ratio'])
            series[name]['highest_vorticity_row']=max(rows,key=lambda r:r['biggest'])
            series[name]['required_linear_constant_on_samples']=max(r['ratio'] for r in rows)
            series[name]['required_log_constant_on_samples']=max(r['ratio']/(1+math.log1p(r['average'])) for r in rows)
    rows=received['stretch-at-max.json']
    series['stretch-at-max.json']=series_summary(rows)
    series['stretch-at-max.json'].update(positive_stretch_samples=sum(r['stretch']>0 for r in rows),
        negative_stretch_samples=sum(r['stretch']<0 for r in rows),
        all_smoothing_negative=all(r['smoothing']<0 for r in rows),
        initial_max=rows[0]['biggest'],final_max=rows[-1]['biggest'],
        largest_saved_max=max(r['biggest'] for r in rows))
    for name,rows in received['pressure-at-peak.json'].items():
        series['pressure-'+name]=series_summary(rows)
        assert all(math.isclose(r['sum'],r['carry']+r['pressure']+r['smoothing'],abs_tol=1e-10) for r in rows)
        series['pressure-'+name].update(negative_pressure_samples=sum(r['pressure']<0 for r in rows),
                                       positive_pressure_samples=sum(r['pressure']>0 for r in rows),
                                       all_smoothing_negative=all(r['smoothing']<0 for r in rows))
    data=received['same-spot.json']
    rows=data['series']
    series['same-spot.json']=series_summary(rows)
    assert data['samples']==len(rows)
    assert all(math.isclose(r['ahead'],r['smoothAtFastest']-r['steepAtFastest'],abs_tol=1e-12) for r in rows)
    assert data['smoothAhead']==sum(r['ahead']>0 for r in rows)
    bound=received['bound-check.json']['n48']
    rows=bound['series']
    series['bound-check.json']=series_summary(rows)
    assert bound['tEnd']==rows[-1]['t'] and bound['integral']==rows[-1]['integral']
    assert bound['vortHigh']==max(r['vortMax'] for r in rows)
    assert all(b['integral']>=a['integral'] for a,b in zip(rows,rows[1:]))
    trapezoid=sum((b['t']-a['t'])*(a['vortMax']+b['vortMax'])/2 for a,b in zip(rows,rows[1:]))
    series['bound-check.json'].update(reported_integral=bound['integral'],
                                     trapezoid_from_saved_samples=trapezoid)
    rate=received['stretch-rate.json']
    rows=rate['series']
    series['stretch-rate.json']=series_summary(rows)
    assert rate['integral']==rows[-1]['integral']
    series['stretch-rate.json'].update(reported_signed_integral=rate['integral'],
        trapezoid_from_saved_samples=sum((b['t']-a['t'])*(a['rate']+b['rate'])/2 for a,b in zip(rows,rows[1:])),
        initial_rate=rows[0]['rate'],final_rate=rows[-1]['rate'])
    other=received['stretch-at-max.json'][0]
    assert math.isclose(rows[0]['rate'],other['stretch']/other['biggest'],rel_tol=1e-12)
    for name in ('high-fraction.json','matched-strain.json','matched-continue.json'):
        rows=received[name]['series']
        series[name]=series_summary(rows)
        assert all(r['biggest']>0 and math.isclose(r['fraction'],r['rate']/r['biggest'],rel_tol=1e-12) for r in rows)
        series[name].update(initial_max=rows[0]['biggest'],final_max=rows[-1]['biggest'],
            highest_saved_max=max(r['biggest'] for r in rows),
            min_fraction=min(r['fraction'] for r in rows),max_fraction=max(r['fraction'] for r in rows),
            negative_fraction_samples=sum(r['fraction']<0 for r in rows))
    continuation=received['matched-continue.json']
    assert continuation['series'][0]==received['matched-strain.json']['series'][0]
    assert continuation['high']>=max(r['biggest'] for r in continuation['series'])
    series['matched-continue.json'].update(reported_high=continuation['high'],
        reported_stop_reason=continuation['reason'],
        initial_row_matches_matched_strain=True,
        reported_high_present_in_saved_rows=any(r['biggest']==continuation['high'] for r in continuation['series']))
    data=received['matched-nostop.json']
    rows=data['series']
    assert all(r['biggest']>0 for r in rows)
    assert rows[0]['biggest']==continuation['series'][0]['biggest']
    series['matched-nostop.json']=series_summary(rows)
    series['matched-nostop.json'].update(initial_max=rows[0]['biggest'],final_max=rows[-1]['biggest'],
        highest_saved_row=max(rows,key=lambda r:r['biggest']),
        saved_decreases=sum(b['biggest']<a['biggest'] for a,b in zip(rows,rows[1:])),
        reported_stop_reason=data['reason'],initial_max_matches_matched_continue=True)
    starts=[]
    if rebuild_starts:
        for n,name in ((48,'spin-ratio.json'),(64,'spin-ratio-64.json')):
            u,ops,report=solver.build(n,psi=solver.pulled_tube,nu=.002)
            kx,ky,kz,_,keep,_=ops
            curl=[solver.deriv(u[2],ky,keep)-solver.deriv(u[1],kz,keep),
                  solver.deriv(u[0],kz,keep)-solver.deriv(u[2],kx,keep),
                  solver.deriv(u[1],kx,keep)-solver.deriv(u[0],ky,keep)]
            mag2=sum(c*c for c in curl)
            maximum=float(np.sqrt(mag2).max())
            rms=float(np.sqrt(mag2.mean()))
            saved=received[name][0]
            assert math.isclose(maximum,saved['biggest'],rel_tol=1e-10)
            assert math.isclose(rms,saved['average'],rel_tol=1e-10)
            pressure_start=received['pressure-at-peak.json']['n'+str(n)][0]
            assert math.isclose(report['peak'],pressure_start['peak'],rel_tol=1e-10)
            if n==64:
                assert math.isclose(report['energy'],received['tube-finer.json']['series'][0]['energy'],rel_tol=1e-10)
            starts.append({'n':n,'reproduced_peak':report['peak'],'reproduced_vorticity_max':maximum,
                           'reproduced_vorticity_rms':rms,'steps_advanced':0})
    return {'scope':'Received JSON arithmetic and starting values only; no long trajectory replay.',
            'solver_sha256':hashlib.sha256((HERE/'solver.py').read_bytes()).hexdigest(),
            'all_numeric_values_finite':True,'sha256':hashes,'series':series,'rebuilt_starts':starts,
            'pressure_status':'Historical values supplied alongside the uncorrected diagnostic. Preserved as received; not corrected result data.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--check-starts',action='store_true')
    args=parser.parse_args()
    result=check_results(args.check_starts)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,allow_nan=False))
