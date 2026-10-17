"""Check available fresh grid pairs without inspecting memory effect sizes."""
import datetime,json
import numpy as np
from final_experiment import ROOT,P,dump,digest
from analysis import read,filename,response,field_response,distance,effective_regimes

def check_available():
    rows=[]
    for regime,cfg in effective_regimes().items():
        seeds=P['smooth_seeds'] if regime=='smooth' else P['forced_seeds']
        for seed in seeds:
            for layout in P['layouts']:
                lo,hi=[filename(seed,regime,layout,n) for n in cfg['grids']]
                if not all(f.with_suffix('.json').exists() for f in [lo,hi]):continue
                a,ma=read(lo,['time','tangent','end_fields','modes']);b,mb=read(hi,['time','tangent','end_fields','modes'])
                ra,rb=response(a,ma),response(b,mb)
                diff=np.stack([ra[method][arm]-rb[method][arm] for method in P['probe_methods'] for arm in P['arms']])
                rmse=float(np.sqrt(np.mean(diff**2)))
                relative=float(distance(field_response(a,ma),a['modes'],field_response(b,mb),b['modes']).max())
                rows.append(dict(seed=seed,regime=regime,layout=layout,response_ftle_rmse=rmse,response_field_relative=relative,
                    passed=rmse<P['gates']['response_ftle_grid_rmse'] and relative<P['gates']['response_field_grid_relative'],
                    source_raw_sha256=[ma['raw_sha256'],mb['raw_sha256']]))
    result=dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),complete_study=False,
        scope='Available paired-grid numerical agreement only; no history-contrast or erasure effect estimates inspected',
        sources={f:digest(ROOT/f) for f in ['numerical_progress.py','analysis.py']},pairs=rows)
    dump(ROOT/'results'/'numerical-progress.json',result)
    print(json.dumps(dict(checked_pairs=len(rows),failed_pairs=[r for r in rows if not r['passed']])),flush=True)
    return result

if __name__=='__main__':check_available()
