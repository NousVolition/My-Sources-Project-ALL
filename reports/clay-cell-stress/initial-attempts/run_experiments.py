"""Reproducible paired experiments; save every attempt, including failure."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,hashlib
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
from cell_hug import Parameters,Mesh,object_outline,geometry
from contact_hug import Mesh as OldMesh
from implicit_solver import integrate
from baseline_metrics import summarize
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
SOURCES=['clay_hug.py','contact_hug.py','cell_hug.py','implicit_solver.py','run_experiments.py','validate_cells.py','baseline_metrics.py']
def hashes():return {s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES}
def run(name,p,end=1,initial=(-8.,8.,27.),method='Radau',max_step=.02,rtol=2e-7,atol=2e-9,old=False,state=None,group='stress',record_step=.04):
    path=DATA/(name+'.json')
    if path.exists():
        row=json.loads(path.read_text());assert row['source_hashes']==hashes(),name+' source changed; do not silently reuse results'
        raw=np.load(DATA/(name+'.npz'));return row,raw['state'][-1]
    m,t,y,f,audit=integrate(p,end=end,max_step=max_step,initial=initial,method=method,rtol=rtol,atol=atol,mesh_type=OldMesh if old else Mesh,state=state,record_step=record_step)
    metrics=summarize(m,t,y) if len(t)>1 else None
    geoms=[geometry(m,m.unpack(v)[0]) for v in y]
    if metrics:metrics.update({k:min(g[k] for g in geoms) for k in ['min_quality','min_quality_ratio','min_altitude_ratio']})
    row=dict(name=name,group=group,parameters=asdict(p),initial=list(initial),requested_end=end,completed_end=float(t[-1]),failure=f,audit=audit,metrics=metrics,old_material=old,record_step=record_step,source_hashes=hashes())
    np.savez_compressed(DATA/(name+'.npz'),time=t,state=y,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p))
    path.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
    print(name,'t',round(t[-1],5),'sec',round(audit['seconds'],2),'calls',audit['rhs_calls'],'restarts',audit['restarted_trials'],'J',round(audit['minimum_swept_area_ratio'],5),'q',round(audit['minimum_swept_quality_ratio'],5),'failure',f,flush=True)
    return row,y[-1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
    assert json.loads((DATA/'validation.json').read_text())['all_passed']
    base=Parameters(preload=.35,alpha=16,rho=160,viscosity=.15)
    if args.pilot:
        for old in [True,False]:
            for feedback in [False,True]:run(('old' if old else 'new')+'-pilot-'+('two' if feedback else 'one'),replace(base,feedback=feedback),end=.2,old=old,group='pilot')
        return
    configs=[('prior-extreme',{},1),('harder32',dict(alpha=32,rho=256,viscosity=.075),2),('harder64',dict(alpha=64,rho=300,viscosity=.0375),1)]
    for label,changes,end in configs:
        for fb in [False,True]:run(label+'-'+('two' if fb else 'one'),replace(base,feedback=fb,**changes),end=end)
    for fb in [False,True]:run('old-method-'+('two' if fb else 'one'),replace(base,feedback=fb),old=True,group='method-control')
    for a in [4,16,32]:
        for fb in [False,True]:run('amplitude'+str(a)+'-'+('two' if fb else 'one'),replace(base,alpha=a,rho=28,viscosity=1.5,feedback=fb),end=4)
    for shape in ['circle','lobed']:
        for fb in [False,True]:run(shape+'-harder-'+('two' if fb else 'one'),replace(base,alpha=32,rho=256,viscosity=.075,object_shape=shape,feedback=fb),end=1)
    for cap in [.0025,.04]:
        for fb in [False,True]:run('capacity'+str(cap)+'-'+('two' if fb else 'one'),replace(base,lorenz_capacity=cap,feedback=fb))
    for strength in [25,400]:
        run('strength'+str(strength),replace(base,area_strength=strength,shape_strength=strength,feedback=False),group='sensitivity')
    for floor in [.15,.35]:
        run('floor'+str(floor),replace(base,area_floor=floor,feedback=False),group='sensitivity')
    for i,initial in enumerate([(1,1,1),(-8,8,27),(1.000001,1,1)]):
        for fb in [False,True]:run('long-s'+str(i)+'-'+('two' if fb else 'one'),replace(base,alpha=8,rho=80,viscosity=.75,feedback=fb),end=24,initial=initial,group='long',max_step=.04)
    for label,p in [('prior-extreme',base),('harder64',replace(base,alpha=64,rho=300,viscosity=.0375))]:
        for method in ['Radau','BDF']:
            for fb in [False,True]:run(label+'-'+method.lower()+'-fine-'+('two' if fb else 'one'),replace(p,feedback=fb),method=method,max_step=.01,rtol=2e-9,atol=2e-11,group='convergence')
    for fb in [False,True]:run('mesh48-hard-'+('two' if fb else 'one'),replace(base,angles=48,layers=5,feedback=fb),end=1,group='mesh')
    for n,l in [(24,3),(48,5),(72,7)]:run('steady-mesh'+str(n),replace(base,angles=n,layers=l,alpha=0,rho=28,viscosity=1.5),end=6,group='forming',max_step=.04)
    for plastic in [False,True]:
        p=replace(base,alpha=0,rho=28,viscosity=1.5,plastic=plastic)
        _,state=run('form-'+str(plastic),p,end=6,group='memory',max_step=.04)
        run('release-'+str(plastic),replace(p,preload=0),end=12,state=state,group='memory',max_step=.04)
    run('old-form',replace(base,alpha=0,rho=28,viscosity=1.5),end=6,old=True,group='forming',max_step=.04)
if __name__=='__main__':main()
