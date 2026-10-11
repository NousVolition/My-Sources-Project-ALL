"""Recompute result metrics, provenance and topology from stored trajectories."""
from pathlib import Path
from dataclasses import asdict
import hashlib,json
import numpy as np
from clay_hug import Parameters, Mesh
from run_stress import summarize

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def read(name):return json.loads((DATA/name).read_text())

def boundary_crossings(x,n):
    ids=np.r_[np.arange(n),np.arange(len(x)-n,len(x))]
    next_ids=np.r_[(np.arange(n)+1)%n,len(x)-n+(np.arange(n)+1)%n]
    a=x[ids];b=x[next_ids];u=b-a
    cross=lambda v,w:v[...,0]*w[...,1]-v[...,1]*w[...,0]
    v=a[None,:,:]-a[:,None,:];other=u[None,:,:]
    den=cross(u[:,None,:],other)
    safe=np.where(abs(den)>1e-12,den,1)
    s=cross(v,other)/safe;t=cross(v,u[:,None,:])/safe
    hit=(abs(den)>1e-12)&(s>1e-8)&(s<1-1e-8)&(t>1e-8)&(t<1-1e-8)
    return int(np.count_nonzero(np.triu(hit,1)))

def main():
    validation=read('validation.json');assert validation['all_passed']
    protocol=read('protocol.json');follow=read('followup.json');mesh_follow=read('feedback-mesh.json')
    hashes={**protocol['sources'],**follow['sources'],**mesh_follow['sources']}
    for name,digest in hashes.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    rows=read('runs.json')+follow['runs']+mesh_follow['runs'];audit=[]
    for r in rows:
        a=np.load(DATA/(r['name']+'.npz'));p=Parameters(**r['parameters']);m=Mesh(p)
        t=a['time'];y=a['state'];assert np.isfinite(y).all() and np.all(np.diff(t)>0),r['name']
        assert np.array_equal(m.tri,a['triangles']) and np.array_equal(m.initial,a['initial_positions'])
        assert y.shape==(len(t),m.size)
        assert np.min(y[:,m.position_size:m.position_size+m.internal_size])>0,r['name']
        if not r['failure']:assert abs(t[-1]-r['end'])<1e-9
        actual=summarize(m,t,y) if len(t)>1 else None
        def compare(actual,stored):
            if isinstance(actual,dict):
                for k in actual:compare(actual[k],stored[k])
            else:assert np.isclose(actual,stored,rtol=1e-10,atol=1e-10),(r['name'],actual,stored)
        if actual:compare(actual,r['metrics'])
        losses=y[:,-3:-1]
        assert np.min(np.diff(losses,axis=0),initial=0)>-1e-7,r['name']
        crossing=max(boundary_crossings(m.unpack(s)[0],p.angles) for s in y)
        audit.append(dict(name=r['name'],samples=len(t),boundary_crossings=crossing,
                          metrics_recomputed=actual is not None,completed=r['failure'] is None))
        print(r['name'],'audited','boundary crossings',crossing,flush=True)
    checks=read('solver-checks.json')+read('followup-checks.json')
    accepted=['ellipse','lobed','independent-DOP853','strong-resolved']
    assert all(next(c for c in checks if c['name']==n)['passed'] for n in accepted)
    for name,base,fine in [('ellipse','step-ellipse-base','step-ellipse-fine'),('lobed','step-lobed-base','step-lobed-fine'),
                           ('strong-resolved','resolved-step-strong-base','resolved-step-strong-fine')]:
        a=np.load(DATA/(base+'.npz'));b=np.load(DATA/(fine+'.npz'))
        assert np.allclose(a['time'],b['time'])
        scale=np.maximum(1,np.max(abs(b['state'][:,:-4]),axis=0))
        error=float(np.max(abs(a['state'][:,:-4]-b['state'][:,:-4])/scale))
        assert np.isclose(error,next(c for c in checks if c['name']==name)['scaled_state_error'])
    independent=np.load(DATA/'independent.npz');fine=np.load(DATA/'step-ellipse-fine.npz')
    scale=np.maximum(1,np.max(abs(independent['state'][:,:-4]),axis=0))
    error=float(np.max(abs(independent['state'][:,:-4]-fine['state'][:,:-4])/scale))
    assert np.isclose(error,next(c for c in checks if c['name']=='independent-DOP853')['scaled_state_error'])
    analysis=read('analysis.json')
    assert analysis['source_sha256']==hashlib.sha256((ROOT/'build_report.py').read_bytes()).hexdigest()
    report=dict(configurations=len(rows),completed=sum(r['completed'] for r in audit),
                stopped=sum(not r['completed'] for r in audit),all_saved_metrics_recomputed=True,
                initial_validation_passed=True,accepted_solver_checks_passed=True,
                boundary_crossing_runs=[r['name'] for r in audit if r['boundary_crossings']],runs=audit,
                note='Boundary intersections checked at stored times, not every internal RK stage; positive cell areas and simple boundaries do not replace a full self-contact method.')
    (DATA/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='runs'},indent=2))

if __name__=='__main__':main()
