"""Resolve strong-drive failures without replacing the original attempts."""
from pathlib import Path
from dataclasses import asdict, replace
import json, hashlib, time
import numpy as np
from clay_hug import Parameters, integrate, object_outline
from run_stress import summarize

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'

def main():
    sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['clay_hug.py','run_stress.py','run_followup.py']}
    rows=[]
    def run(name,p,end,dt):
        started=time.perf_counter()
        mesh,t,y,f=integrate(p,end=end,dt=dt)
        row=dict(name=name,group='followup',parameters=asdict(p),initial=[1,1,1],end=end,dt=dt,
                 failure=f,seconds=time.perf_counter()-started,metrics=summarize(mesh,t,y) if len(t)>1 else None)
        np.savez_compressed(DATA/(name+'.npz'),time=t,state=y,triangles=mesh.tri,initial_positions=mesh.initial,outline=object_outline(p))
        rows.append(row)
        (DATA/'followup.json').write_text(json.dumps(dict(sources=sources,runs=rows),indent=2,allow_nan=False)+'\n')
        print(name,'end',t[-1],'failure',f,'energy error',row['metrics']['total_budget_error'],flush=True)
        return t,y,row
    p=Parameters(preload=.35,alpha=2)
    for feedback in [False,True]:
        run('resolved-alpha2-'+('two' if feedback else 'one'),replace(p,feedback=feedback),16,.00025)
    a=run('resolved-step-strong-base',p,4,.00025)
    b=run('resolved-step-strong-fine',p,4,.000125)
    scale=np.maximum(1,np.max(abs(b[1][:,:-4]),axis=0))
    error=float(np.max(abs(a[1][:,:-4]-b[1][:,:-4])/scale))
    check=dict(name='strong-resolved',scaled_state_error=error,passed=bool(error<1e-4 and a[2]['failure'] is None and b[2]['failure'] is None))
    (DATA/'followup-checks.json').write_text(json.dumps([check],indent=2)+'\n')
    print(check,flush=True)
    if not check['passed']:raise SystemExit(1)

if __name__=='__main__':main()
