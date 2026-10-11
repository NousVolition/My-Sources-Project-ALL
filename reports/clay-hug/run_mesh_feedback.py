"""Check whether the reciprocal response persists on a finer material mesh."""
from pathlib import Path
from dataclasses import asdict,replace
import hashlib,json,time
import numpy as np
from clay_hug import Parameters,integrate,object_outline
from run_stress import summarize
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'

def main():
    sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['clay_hug.py','run_stress.py','run_mesh_feedback.py']}
    rows=[];base=Parameters(preload=.35,angles=48,layers=5)
    for feedback in [False,True]:
        p=replace(base,feedback=feedback);name='feedback-mesh48-'+('two' if feedback else 'one');start=time.perf_counter()
        m,t,y,f=integrate(p,end=16,dt=.0005)
        row=dict(name=name,group='mesh-feedback',parameters=asdict(p),initial=[1,1,1],end=16,dt=.0005,
                 failure=f,seconds=time.perf_counter()-start,metrics=summarize(m,t,y))
        np.savez_compressed(DATA/(name+'.npz'),time=t,state=y,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p))
        rows.append(row)
        (DATA/'feedback-mesh.json').write_text(json.dumps(dict(sources=sources,runs=rows),indent=2,allow_nan=False)+'\n')
        print(name,'completed',t[-1],'failure',f,'energy error',row['metrics']['total_budget_error'],flush=True)
    if any(r['failure'] or max(r['metrics']['total_budget_error'],r['metrics']['material_budget_error'])>1e-3 for r in rows):raise SystemExit(1)

if __name__=='__main__':main()
