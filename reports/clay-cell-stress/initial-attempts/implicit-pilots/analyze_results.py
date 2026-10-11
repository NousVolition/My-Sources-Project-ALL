"""Audit every recording and quantify matched finite-window effects."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import trapezoid
from cell_hug import Parameters,Mesh,geometry,object_gap
from contact_hug import Mesh as OldMesh
from baseline_metrics import summarize
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def rows():
    return [r for p in sorted(DATA.glob('*.json')) if isinstance(r:=json.loads(p.read_text()),dict) and 'source_hashes' in r and 'name' in r]
def mesh(row):return (OldMesh if row['old_material'] else Mesh)(Parameters(**row['parameters']))
def motion(row,start=0,end=None):
    raw=np.load(DATA/(row['name']+'.npz'));t=raw['time'];y=raw['state'];end=t[-1] if end is None else end
    indices=np.where((t>=start-1e-10)&(t<=end+1e-10))[0];a,b=indices[0],indices[-1]
    duration=t[b]-t[a];m=mesh(row)
    if duration<=0:return None
    z=y[indices,-5];tt=t[indices];zmean=trapezoid(z,tt)/duration
    return dict(material_motion_rms=float(np.sqrt(max(0,(y[b,-3]-y[a,-3])/(m.p.viscosity*np.sum(m.mass)*duration)))),z_std=float(np.sqrt(max(0,trapezoid((z-zmean)**2,tt)/duration))))
def state_error(a,b):
    ra=np.load(DATA/(a+'.npz'));rb=np.load(DATA/(b+'.npz'))
    # Refinement comparisons require exactly matching stored times and exclude integrated bookkeeping.
    assert np.allclose(ra['time'],rb['time'],atol=1e-10,rtol=0)
    x=ra['state'][:,:-4];y=rb['state'][:,:-4];scale=np.maximum(1,np.max(abs(y),axis=0))
    return float(np.max(abs(x-y)/scale))
def main():
    allrows=rows();lookup={r['name']:r for r in allrows};audits=[]
    for row in allrows:
        raw=np.load(DATA/(row['name']+'.npz'));t,y=raw['time'],raw['state'];m=mesh(row)
        assert np.isfinite(y).all() and np.all(np.diff(t)>0)
        for n,h in row['source_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h
        crossings=0;probe=float('inf');aa,bb=np.triu_indices(len(m.edges),1)
        keep=(m.i[aa]!=m.i[bb])&(m.i[aa]!=m.j[bb])&(m.j[aa]!=m.i[bb])&(m.j[aa]!=m.j[bb]);aa,bb=aa[keep],bb[keep]
        cross=lambda a,b:a[:,0]*b[:,1]-a[:,1]*b[:,0]
        f=np.linspace(0,1,33)
        for s in y:
            x=m.unpack(s)[0];a,b=x[m.i],x[m.j]
            pts=((1-f[None,:,None])*a[:,None]+f[None,:,None]*b[:,None]).reshape(-1,2)
            probe=min(probe,float(np.min(object_gap(pts,m.p)[0])))
            u,v,w,z=a[aa],b[aa],a[bb],b[bb]
            crossings+=int(np.sum((cross(v-u,w-u)*cross(v-u,z-u)<-1e-20)&(cross(z-w,u-w)*cross(z-w,v-w)<-1e-20)))
        audit=row['audit'];metrics=row['metrics'];complete=row['failure'] is None and abs(t[-1]-row['requested_end'])<1e-9
        budget=max(metrics[k] for k in ['material_budget_error','lorenz_budget_error','total_budget_error']) if metrics else None
        positive=audit['minimum_swept_clearance']>0 and audit['minimum_swept_self_clearance']>0 and audit['minimum_swept_area_ratio']>(.01 if row['old_material'] else m.p.area_floor) and audit['minimum_swept_quality_ratio']>(0 if row['old_material'] else m.p.shape_floor)
        passed=complete and positive and probe>=0 and crossings==0 and budget is not None and budget<1e-3 and metrics['maximum_wall_overlap']==0 and np.all(np.diff(y[:,-3:-1],axis=0)>=-1e-9)
        audits.append(dict(name=row['name'],completed=complete,passed=bool(passed),geometry_passed=bool(positive and probe>=0 and crossings==0),minimum_probe_gap=probe,proper_edge_crossings=crossings,maximum_budget_error=budget,raw_sha256=hashlib.sha256((DATA/(row['name']+'.npz')).read_bytes()).hexdigest()))
    pairs=[];blocks=[]
    for name,r in lookup.items():
        if not name.endswith('-one') or name[:-4]+'-two' not in lookup:continue
        b=lookup[name[:-4]+'-two']
        if r['failure'] or b['failure']:continue
        start=r['requested_end']/2;a_stats=motion(r,start);b_stats=motion(b,start)
        changes={k:100*(b_stats[k]/a_stats[k]-1) if a_stats[k]>1e-14 else None for k in a_stats}
        pairs.append(dict(key=name[:-4],start=start,end=r['requested_end'],one=a_stats,two=b_stats,changes_percent=changes,
            one_gap=r['metrics']['final']['surface_gap_rms'],two_gap=b['metrics']['final']['surface_gap_rms'],
            forward=b['audit']['forward_work_step_trapezoid'],returned=b['audit']['return_work_step_trapezoid']))
        if r['group']=='long':
            for start in [4,10,16]:
                a_stats=motion(r,start,start+6);b_stats=motion(b,start,start+6)
                blocks.append(dict(key=name[:-4],start=start,end=start+6,changes_percent={k:100*(b_stats[k]/a_stats[k]-1) for k in a_stats}))
    checks=[]
    for label in ['prior-extreme','harder64']:
        for method in ['radau','bdf']:
            for direction in ['one','two']:
                a=label+'-'+direction;b=label+'-'+method+'-fine-'+direction
                if a in lookup and b in lookup and not lookup[a]['failure'] and not lookup[b]['failure']:
                    error=state_error(a,b);checks.append(dict(name=b,reference=a,error=error,target=1e-4,passed=error<1e-4))
    for direction in ['one','two']:
        a='old-method-'+direction;b='prior-extreme-'+direction
        if a in lookup and not lookup[a]['failure']:
            error=state_error(a,b);checks.append(dict(name='old-law-method-control-'+direction,error=error,target=1e-4,passed=error<1e-4))
    prior=[]
    for direction in ['one','two']:
        raw=np.load(DATA/('prior-extreme-'+direction+'.npz'));oldrow=json.loads((DATA/('prior-extreme-'+direction+'.json')).read_text())
        m=OldMesh(Parameters(**oldrow['parameters']));values=[geometry(m,m.unpack(s)[0]) for s in raw['state']]
        prior.append(dict(direction=direction,min_quality_ratio=min(v['min_quality_ratio'] for v in values),min_altitude_ratio=min(v['min_altitude_ratio'] for v in values),minimum_saved_area_ratio=min(v['min_area_ratio'] for v in values),minimum_swept_area_ratio=oldrow['audit']['minimum_swept_area_ratio']))
    out=dict(configurations=len(allrows),completed=sum(a['completed'] for a in audits),audits=audits,pairs=pairs,long_blocks=blocks,solver_checks=checks,prior_baseline=prior,
        failures=[dict(name=r['name'],failure=r['failure'],completed_end=r['completed_end'],requested_end=r['requested_end']) for r in allrows if r['failure']],
        maximum_completed_budget_error=max((a['maximum_budget_error'] for a in audits if a['completed']),default=None),
        all_completed_audits_passed=all(a['passed'] for a in audits if a['completed']),all_solver_checks_passed=bool(checks) and all(c['passed'] for c in checks),
        notes='Finite windows and three starts provide observed variation, not confidence intervals, global stability or experimental calibration.')
    (DATA/'all-runs.json').write_text(json.dumps(allrows,indent=2)+'\n');(DATA/'analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['audits','pairs','long_blocks']},indent=2),flush=True)
if __name__=='__main__':main()
