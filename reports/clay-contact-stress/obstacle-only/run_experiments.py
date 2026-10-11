"""All experiments are specified before running; retain failures and raw states."""
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

def configurations():
    base=Parameters(preload=.35);out=[]
    def add(key,p,initial=(-8.,8.,27.),end=16,group='paired',**kw):
        out.append(dict(name=key,parameters=asdict(p),initial=initial,end=end,group=group,dt=kw.pop('dt',.004),**kw))
    def pair(key,p,**kw):
        for f in [False,True]:add(key+('-two' if f else '-one'),replace(p,feedback=f),**kw)
    pair('moderate',base,initial=(1.,1.,1.))
    pair('repair-seed0',replace(base,alpha=2),initial=(1.,1.,1.))
    pair('repair-seed1',replace(base,alpha=2))
    for alpha in [4,8,16]:pair('amplitude'+str(alpha),replace(base,alpha=alpha),end=8 if alpha==16 else 16)
    pair('rho100',replace(base,rho=100,alpha=2))
    pair('rho160',replace(base,rho=160,alpha=4))
    pair('fast',replace(base,viscosity=.15,alpha=4),end=8)
    pair('stiff',replace(base,bulk=50,contact=100000,alpha=4),end=8)
    for c in [.0025,.04]:pair('capacity'+str(c),replace(base,lorenz_capacity=c,alpha=4))
    for shape in ['circle','lobed']:pair(shape,replace(base,object_shape=shape,alpha=4))
    pair('third-seed',replace(base,alpha=4),initial=(15.,-10.,40.))
    for seed,initial in enumerate([(1.,1.,1.),(-8.,8.,27.),(15.,-10.,40.)]):
        pair('long-seed'+str(seed),replace(base,alpha=4),initial=initial,end=64,group='long')
    pair('mesh48',replace(base,angles=48,layers=5,alpha=4),group='mesh')
    for n,l in [(24,3),(48,5),(72,7)]:
        add('steady-mesh'+str(n),replace(base,alpha=0,angles=n,layers=l),end=12,group='steady')
    for s in [8,16,32]:add('lobed-sub'+str(s),replace(base,alpha=0,object_shape='lobed',contact_subdivisions=s),end=8,group='contact-resolution')
    pair('adversarial',replace(base,alpha=16,rho=160,viscosity=.15),end=2,dt=.2,group='adversarial')
    for key,p in [('strong',replace(base,alpha=8)),('highrho',replace(base,alpha=4,rho=160))]:
        for name,dt,rtol,atol in [('base',.004,2e-7,2e-9),('fine',.002,2e-9,2e-11)]:
            add('refine-'+key+'-'+name,p,end=4,dt=dt,rtol=rtol,atol=atol,group='convergence')
    pair('zero-port',replace(base,alpha=0),end=4,group='control')
    return out

def run(config):
    started=time.perf_counter();p=Parameters(**config['parameters'])
    m,t,y,f,a=integrate(p,end=config['end'],dt=config['dt'],initial=config['initial'],
        rtol=config.get('rtol',2e-7),atol=config.get('atol',2e-9))
    row=dict(**config,failure=f,seconds=time.perf_counter()-started,audit=a,metrics=summarize(m,t,y),
             model_sha256=hashlib.sha256((ROOT/'contact_hug.py').read_bytes()).hexdigest())
    row['metrics']['max_pressure']=float(np.max(abs(p.alpha/10*y[:,-7])))
    row['metrics']['peak_contact_fraction']=max(m.diagnostics(s)['contact_fraction'] for s in y)
    np.savez_compressed(DATA/(config['name']+'.npz'),time=t,state=y,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p))
    write(config['name']+'.json',row)
    return row

def main():
    DATA.mkdir(exist_ok=True)
    assert json.loads((DATA/'validation.json').read_text())['all_passed']
    configs=configurations()
    write('protocol.json',dict(configurations=configs,sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['clay_hug.py','contact_hug.py','baseline_metrics.py','validate_contact.py','run_experiments.py']},
        targets=dict(overlap=0,swept_clearance_positive=True,energy_budget=1e-3,short_scaled_state_error=1e-4),
        scope='Dimensionless 2-D annular material, fixed rigid object, no friction, fracture or calibrated clay data.'))
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(run,c):c['name'] for c in configs if not (DATA/(c['name']+'.json')).exists()}
        for future in as_completed(futures):
            r=future.result();print(r['name'],'finished',r['failure'],'clearance',r['audit']['minimum_swept_clearance'],'budget',r['metrics']['total_budget_error'],flush=True)
    rows=[json.loads((DATA/(c['name']+'.json')).read_text()) for c in configs];write('runs.json',rows)
    print('Saved',len(rows),'configurations',flush=True)
if __name__=='__main__':main()
