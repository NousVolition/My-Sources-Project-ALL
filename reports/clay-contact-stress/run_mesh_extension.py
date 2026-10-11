"""Continue both finer-mesh runs beyond the time of the discovered fold."""
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,time,hashlib,sys
import numpy as np
from contact_hug import Parameters,integrate,object_outline
from run_experiments import ROOT,DATA,write,measurements
def run(direction):
    parent='mesh48-'+direction;name='extended-'+parent
    if (DATA/(name+'.json')).exists():return json.loads((DATA/(name+'.json')).read_text())
    while not (DATA/(parent+'.json')).exists():
        if '--wait' not in sys.argv:raise FileNotFoundError('Run the main protocol before extending '+parent)
        time.sleep(5)
    row=json.loads((DATA/(parent+'.json')).read_text());raw=np.load(DATA/(parent+'.npz'))
    assert row['failure'] is None and abs(raw['time'][-1]-4)<1e-9
    p=Parameters(**row['parameters']);started=time.perf_counter()
    m,t,y,f,a=integrate(p,end=2,dt=.004,state=raw['state'][-1])
    y[:,-4:]+=raw['state'][-1,-4:]
    tt=np.r_[raw['time'],t[1:]+4];yy=np.vstack([raw['state'],y[1:]])
    audit=dict(a)
    for k in ['minimum_step','minimum_accepted_clearance','minimum_swept_clearance','minimum_swept_area_ratio','minimum_swept_self_clearance']:audit[k]=min(a[k],row['audit'][k])
    for k in ['accepted_half_steps','rejected_trials','domain_rejections','forward_work_step_trapezoid','return_work_step_trapezoid']:audit[k]=a[k]+row['audit'][k]
    if f:f=dict(f,time=f['time']+4)
    result=dict(name=name,group='mesh-extension',parameters=asdict(p),initial=row['initial'],end=6,dt=.004,failure=f,
        seconds=time.perf_counter()-started,source_run=parent,source_npz_sha256=hashlib.sha256((DATA/(parent+'.npz')).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((ROOT/'contact_hug.py').read_bytes()).hexdigest(),audit=audit,metrics=measurements(m,tt,yy),
        note='Full six-unit trajectory assembled from the first four units plus a two-unit autonomous continuation. Not an independent initial-state replicate.')
    np.savez_compressed(DATA/(name+'.npz'),time=tt,state=yy,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p));write(name+'.json',result)
    return result
def main():
    write('mesh-extension-protocol.json',dict(end=6,continuation_start=4,continuation_duration=2,configurations=['extended-mesh48-one','extended-mesh48-two'],
        sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['contact_hug.py','run_experiments.py','run_mesh_extension.py']},
        scope='Two finer-mesh continuations to test the repaired folding interval.'))
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(run,d) for d in ['one','two']]
        for future in as_completed(futures):
            r=future.result();print(r['name'],r['failure'],'self gap',r['audit']['minimum_swept_self_clearance'],flush=True)
    write('mesh-extension.json',dict(runs=[json.loads((DATA/('extended-mesh48-'+d+'.json')).read_text()) for d in ['one','two']]))
if __name__=='__main__':main()
