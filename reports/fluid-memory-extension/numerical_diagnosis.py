"""Compare numerical conditions only, without examining history-effect contrasts."""
import json
import numpy as np
from analysis import read,filename,response,field_response,distance
from experiment import ROOT,P,dump

def diagnose():
    out=[]
    for seed in P['pilot_seeds']:
        fine=filename(seed,'forced_nu002','coincident',96,.002,'fast-data')
        b,mb=read(fine)
        for n,folder in [(64,'fast-data'),(80,'production-data')]:
            a,ma=read(filename(seed,'forced_nu002','coincident',n,.002,folder))
            ra,rb=response(a,ma),response(b,mb)
            rows=[]
            dr=distance(field_response(a,ma),a['modes'],field_response(b,mb),b['modes']).reshape(3,7,2)
            for j,method in enumerate(P['probe_methods']):
                for k,arm in enumerate(P['arms']):
                    rows.append(dict(method=method,arm=arm,response_ftle_rmse=float(np.sqrt(np.mean((ra[method][arm]-rb[method][arm])**2))),field_relative=float(dr[j,k].max())))
            start=distance(a['start_fields'],a['modes'],b['start_fields'],b['modes'])
            coeff=abs(a['pulse_coefficients']-b['pulse_coefficients'])
            out.append(dict(seed=seed,n=n,reference_n=96,start_field_relative_max=float(start.max()),pulse_coefficient_absolute_max=float(coeff.max()),rows=rows))
    result=dict(effect_sizes_inspected=False,numerical_comparisons=out)
    dump(ROOT/'results'/'numerical-diagnosis.json',result)
    print(json.dumps(result,indent=2))
if __name__=='__main__':diagnose()
