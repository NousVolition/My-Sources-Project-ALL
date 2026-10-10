"""Reproduce the direct-phase precision failure on two validation starts."""
from pathlib import Path
import csv,json
import numpy as np
from scipy.integrate import solve_ivp
from model import integrate,integrate_state,to_state,from_state,stable_rhs,voltage
ROOT=Path(__file__).resolve().parent


def main():
    original=json.loads((ROOT/'audit/initial_independent_validation.json').read_text())
    cases=[r for r in original['independent_solver_runs'] if r['alpha']==2 and r['bias']==1.005]
    rows=[];records={}
    for j,row in enumerate(cases):
        p0=np.array(row['initial']);s0,sign=to_state(p0)
        ref=solve_ivp(lambda t,y:stable_rhs(y,1.005,2.),(0,400),s0,method='DOP853',rtol=1e-12,atol=1e-13,max_step=.1,t_eval=np.arange(401))
        pr=from_state(ref.y.T,sign);vr=voltage(pr[200],pr[400],200)
        records[f'reference_{j}']=pr;records[f'log_state_{j}']=ref.y.T
        for dt in [.05,.025,.0125]:
            _,p,_=integrate(p0,1.005,2.,400,dt,stride=1.)
            vd=voltage(p[200],p[400],200)
            rows.append({'case':j,'dt':dt,'direct_v1':float(vd[0]),'direct_v2':float(vd[1]),
                         'reference_v1':float(vr[0]),'reference_v2':float(vr[1]),
                         'total_error':float(abs((vd-vr).sum())),
                         'smallest_log_separation_coordinate':float(ref.y[1].min()),
                         'largest_log_separation_coordinate':float(ref.y[1].max())})
            records[f'direct_{j}_{dt}']=p
    with (ROOT/'results/precision_diagnostic.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    np.savez_compressed(ROOT/'data/precision_diagnostic.npz',time=np.arange(401),**records)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
