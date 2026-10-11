"""Frozen final protocol with obstacle, wall and self-contact protection."""
from pathlib import Path
from dataclasses import asdict,replace
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,hashlib,time
import numpy as np
from contact_hug import Parameters,Mesh,integrate,object_outline
from baseline_metrics import summarize
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def write(name,data):
    (DATA/name).write_text(json.dumps(data,indent=2,allow_nan=False,default=lambda x:x.item())+'\n')
def weighted_motion(m,t,y,start=None):
    start=min(4,t[-1]/2) if start is None else start
    idx=np.flatnonzero(t>=start)[0]
    return float(np.sqrt(max(0,(y[-1,-3]-y[idx,-3])/(m.p.viscosity*sum(m.mass)*(t[-1]-t[idx])))))
def measurements(m,t,y):
    metrics=summarize(m,t,y);metrics['material_motion_rms']=weighted_motion(m,t,y)
    metrics['max_pressure']=float(np.max(abs(m.p.alpha/10*y[:,-7])))
    metrics['peak_contact_fraction']=max(m.diagnostics(s)['contact_fraction'] for s in y)
    return metrics
def configurations():
    base=Parameters(preload=.35);out=[]
    def add(key,p,initial=(-8.,8.,27.),end=6,group='paired',**kw):
        out.append(dict(name=key,parameters=asdict(p),initial=initial,end=end,group=group,dt=kw.pop('dt',.004),**kw))
    def pair(key,p,**kw):
        for f in [False,True]:add(key+('-two' if f else '-one'),replace(p,feedback=f),**kw)
    pair('repair-seed0',replace(base,alpha=2),initial=(1.,1.,1.),end=8)
    pair('amplitude4',replace(base,alpha=4))
    pair('amplitude8',replace(base,alpha=8),end=8)
    pair('amplitude16',replace(base,alpha=16),end=4)
    pair('rho160',replace(base,rho=160,alpha=4),end=4)
    pair('fast',replace(base,viscosity=.15,alpha=4),end=2)
    pair('stiff',replace(base,bulk=50,contact=100000,alpha=4),end=2)
    for c in [.0025,.04]:pair('capacity'+str(c),replace(base,lorenz_capacity=c,alpha=4))
    for shape in ['circle','lobed']:pair(shape,replace(base,object_shape=shape,alpha=4))
    for seed,initial in enumerate([(1.,1.,1.),(-8.,8.,27.),(15.,-10.,40.)]):
        pair('long-seed'+str(seed),replace(base,alpha=4),initial=initial,end=24,group='long')
    pair('mesh48',replace(base,angles=48,layers=5,alpha=4),end=4,group='mesh')
    for n,l in [(24,3),(48,5),(72,7)]:add('steady-mesh'+str(n),replace(base,alpha=0,angles=n,layers=l),end=6,group='steady')
    for s in [8,16,32]:add('lobed-sub'+str(s),replace(base,alpha=0,object_shape='lobed',contact_subdivisions=s),end=4,group='contact-resolution')
    pair('adversarial',replace(base,alpha=16,rho=160,viscosity=.15),end=1,dt=.2,group='adversarial')
    for key,p in [('strong',replace(base,alpha=8)),('highrho',replace(base,alpha=4,rho=160))]:
        for name,dt,rtol,atol in [('base',.004,2e-7,2e-9),('fine',.002,2e-9,2e-11)]:
            add('refine-'+key+'-'+name,p,end=2,dt=dt,rtol=rtol,atol=atol,group='convergence')
    pair('zero-port',replace(base,alpha=0),end=2,group='control')
    return out
def run(config):
    started=time.perf_counter();p=Parameters(**config['parameters'])
    m,t,y,f,a=integrate(p,end=config['end'],dt=config['dt'],initial=config['initial'],rtol=config.get('rtol',2e-7),atol=config.get('atol',2e-9))
    row=dict(**config,failure=f,seconds=time.perf_counter()-started,audit=a,metrics=measurements(m,t,y),
        model_sha256=hashlib.sha256((ROOT/'contact_hug.py').read_bytes()).hexdigest())
    np.savez_compressed(DATA/(config['name']+'.npz'),time=t,state=y,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p))
    write(config['name']+'.json',row);return row
def main():
    assert json.loads((DATA/'validation.json').read_text())['all_passed'];configs=configurations()
    write('protocol.json',dict(configurations=configs,sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['clay_hug.py','contact_hug.py','baseline_metrics.py','validate_contact.py','run_experiments.py']},
        targets=dict(overlap=0,proper_crossings=0,swept_clearance_positive=True,swept_self_clearance_positive=True,energy_budget=1e-3,short_scaled_state_error=1e-4),record_step=.04,
        scope='Dimensionless 2-D annular material, fixed object and square, frictionless self-contact. Not calibrated clay.'))
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(run,c):c['name'] for c in configs if not (DATA/(c['name']+'.json')).exists()}
        for future in as_completed(futures):
            r=future.result();print(r['name'],'finished',r['failure'],'self gap',r['audit']['minimum_swept_self_clearance'],'budget',r['metrics']['total_budget_error'],flush=True)
    write('runs.json',[json.loads((DATA/(c['name']+'.json')).read_text()) for c in configs]);print('Saved',len(configs),'configurations',flush=True)
if __name__=='__main__':main()
