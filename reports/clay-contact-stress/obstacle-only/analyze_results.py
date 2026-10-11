"""Recompute completed results from saved states, never assume success."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib,json,sys
import numpy as np
from scipy.integrate import trapezoid
from contact_hug import *
from baseline_metrics import summarize
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'

def audit(row):
    m=Mesh(Parameters(**row['parameters']));raw=np.load(DATA/(row['name']+'.npz'));t=raw['time'];y=raw['state']
    assert np.isfinite(y).all() and np.all(np.diff(t)>0)
    assert np.min(np.diff(y[:,-3:-1],axis=0))>=-1e-10 # Dissipation cannot decrease.
    assert row['model_sha256']==hashlib.sha256((ROOT/'contact_hug.py').read_bytes()).hexdigest()
    values=summarize(m,t,y)
    error=max(abs(values[k]-v)/max(1,abs(v)) for k,v in row['metrics'].items() if k in values and isinstance(v,(int,float)))
    n=m.p.angles;crossings=0;minimum_gap=float('inf')
    # Every saved frame: all edge probes plus all pairs of nonadjacent mesh edges.
    aa,bb=np.triu_indices(len(m.edges),1)
    keep=(m.i[aa]!=m.i[bb])&(m.i[aa]!=m.j[bb])&(m.j[aa]!=m.i[bb])&(m.j[aa]!=m.j[bb]);aa=aa[keep];bb=bb[keep]
    def cross(a,b):return a[:,0]*b[:,1]-a[:,1]*b[:,0]
    f=np.linspace(0,1,65)
    for state in y:
        x=m.unpack(state)[0];a=x[m.i];b=x[m.j]
        points=((1-f[None,:,None])*a[:,None]+f[None,:,None]*b[:,None]).reshape(-1,2)
        minimum_gap=min(minimum_gap,float(np.min(object_gap(points,m.p)[0])))
        u,v,w,z=a[aa],b[aa],a[bb],b[bb]
        c1=cross(v-u,w-u);c2=cross(v-u,z-u);c3=cross(z-w,u-w);c4=cross(z-w,v-w)
        crossings+=int(np.count_nonzero((c1*c2<-1e-20)&(c3*c4<-1e-20)))
    accepted=row['audit']
    passed=(row['failure'] is None and abs(t[-1]-row['end'])<1e-9 and error<1e-12 and minimum_gap>=0 and crossings==0
        and accepted['minimum_swept_clearance']>0 and accepted['minimum_swept_area_ratio']>.01
        and values['maximum_wall_overlap']==0 and max(values['total_budget_error'],values['material_budget_error'],values['lorenz_budget_error'])<1e-3)
    result=dict(name=row['name'],recomputed_metric_error=error,minimum_probe_radial_gap=minimum_gap,
        proper_edge_crossings=crossings,passed=bool(passed),raw_sha256=hashlib.sha256((DATA/(row['name']+'.npz')).read_bytes()).hexdigest(),
        auditor_sha256=hashlib.sha256((ROOT/'analyze_results.py').read_bytes()).hexdigest())
    cache=DATA/'audits';cache.mkdir(exist_ok=True);(cache/(row['name']+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    return result

def audited(row):
    cache=DATA/'audits'/(row['name']+'.json')
    if cache.exists():
        result=json.loads(cache.read_text())
        if result['raw_sha256']==hashlib.sha256((DATA/(row['name']+'.npz')).read_bytes()).hexdigest() and result['auditor_sha256']==hashlib.sha256((ROOT/'analyze_results.py').read_bytes()).hexdigest():return result
    return audit(row)

def statistics(row,start=4,end=None):
    m=Mesh(Parameters(**row['parameters']));raw=np.load(DATA/(row['name']+'.npz'));t=raw['time'];y=raw['state']
    use=(t>=start)&(t<=(row['end'] if end is None else end));t=t[use];y=y[use]
    mean=lambda a:float(trapezoid(a,t)/(t[-1]-t[0]))
    ly=y[:,-7:-4];zmean=mean(ly[:,2]);motion=np.array([np.sum(m.rhs(tt,s)[:m.position_size]**2)/m.n for tt,s in zip(t,y)])
    imprint=np.array([m.diagnostics(s)['imprint'] for s in y])
    return dict(lorenz_rms=float(np.sqrt(mean(np.sum(ly[:,:2]**2,axis=1)))),z_std=float(np.sqrt(max(0,mean((ly[:,2]-zmean)**2)))),
        mesh_motion_rms=float(np.sqrt(mean(motion))),mean_imprint=mean(imprint))

def main():
    if '--audit-ready' in sys.argv:
        rows=[json.loads(f.read_text()) for f in DATA.glob('*.json') if f.name not in ['runs.json','all-runs.json']]
        rows=[r for r in rows if isinstance(r,dict) and 'metrics' in r and 'name' in r]
        with ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(audited,rows))
        print('Audited',len(results),'completed recordings; failed:',[r['name'] for r in results if not r['passed']],flush=True)
        if not all(r['passed'] for r in results):raise SystemExit(1)
        return
    for file in ['protocol.json','sensitivity-protocol.json']:
        protocol=json.loads((DATA/file).read_text())
        for name,digest in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    rows=json.loads((DATA/'runs.json').read_text())+json.loads((DATA/'sensitivity.json').read_text())['runs'];lookup={r['name']:r for r in rows}
    (DATA/'all-runs.json').write_text(json.dumps(rows,indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=4) as pool:audits=list(pool.map(audited,rows))
    comparisons=[]
    for r in rows:
        if not r['name'].endswith('-one'):continue
        key=r['name'][:-4];other=lookup[key+'-two']
        a=r['metrics'];b=other['metrics']
        changes={k:100*(b[k]/a[k]-1) if abs(a[k])>1e-14 else None for k in ['lorenz_rms','z_std','mesh_motion_rms','mean_imprint']}
        comparisons.append(dict(key=key,group=r['group'],changes_percent=changes,
            one_final_radial_gap_rms=a['final']['surface_gap_rms'],two_final_radial_gap_rms=b['final']['surface_gap_rms'],
            one_forward=r['audit']['forward_work_step_trapezoid'],one_return=r['audit']['return_work_step_trapezoid'],
            two_forward=other['audit']['forward_work_step_trapezoid'],two_return=other['audit']['return_work_step_trapezoid']))
    long=[]
    for seed in range(3):
        key='long-seed'+str(seed)
        for start in [4,24,44]:
            a=statistics(lookup[key+'-one'],start,start+20);b=statistics(lookup[key+'-two'],start,start+20)
            long.append(dict(seed=seed,start=start,end=start+20,changes_percent={k:100*(b[k]/a[k]-1) for k in a}))
    old=json.loads((ROOT.parent/'clay-hug/data/followup.json').read_text())['runs']
    baseline=[]
    for direction in ['one','two']:
        r=next(v for v in old if v['name']=='resolved-alpha2-'+direction);new=lookup['repair-seed0-'+direction]
        baseline.append(dict(direction=direction,old_overlap=r['metrics']['maximum_penetration'],new_overlap=new['metrics']['maximum_penetration'],
            new_certified_minimum_clearance=new['audit']['minimum_accepted_clearance']))
    # Same continuous model: gross work is auxiliary step trapezoid, check refinement.
    gross_checks=[]
    for key in ['strong','highrho']:
        a=lookup['refine-'+key+'-base']['audit'];b=lookup['refine-'+key+'-fine']['audit']
        gross_checks.append(dict(case=key,relative_difference={k:abs(a[k]-b[k])/max(1,abs(b[k])) for k in ['forward_work_step_trapezoid','return_work_step_trapezoid']}))
    steady=[]
    for n in [24,48,72]:
        r=lookup['steady-mesh'+str(n)];steady.append(dict(angles=n,layers=r['parameters']['layers'],final_imprint=r['metrics']['final']['imprint'],gap=r['metrics']['final']['surface_gap_rms']))
    solver=json.loads((DATA/'solver-checks.json').read_text())
    result=dict(configurations=len(rows),finished=sum(r['failure'] is None for r in rows),all_audits_passed=all(a['passed'] for a in audits),
        all_independent_checks_passed=solver['all_passed'],audit=audits,comparisons=comparisons,long_blocks=long,baseline=baseline,
        gross_work_refinement=gross_checks,steady_mesh=steady,
        minimum_swept_clearance=min(r['audit']['minimum_swept_clearance'] for r in rows),
        maximum_budget_error=max(max(r['metrics'][k] for k in ['total_budget_error','material_budget_error','lorenz_budget_error']) for r in rows),
        failures=[r['name'] for r in rows if r['failure']],source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['contact_hug.py','analyze_results.py','independent_checks.py']})
    (DATA/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['audit','comparisons','long_blocks','source_sha256']},indent=2),flush=True)
    if not result['all_audits_passed'] or not solver['all_passed']:raise SystemExit(1)
if __name__=='__main__':main()
