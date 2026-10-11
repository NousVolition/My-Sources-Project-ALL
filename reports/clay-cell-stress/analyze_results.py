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
from source_versions import resolve_source
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
def prefix_error(a,b):
    ra=np.load(DATA/(a+'.npz'));rb=np.load(DATA/(b+'.npz'))
    ids=[int(np.argmin(abs(ra['time']-t))) for t in rb['time']]
    keep=abs(ra['time'][ids]-rb['time'])<1e-10
    x=ra['state'][np.array(ids)[keep],:-4];y=rb['state'][keep,:-4]
    return float(np.max(abs(x-y)/np.maximum(1,np.max(abs(y),axis=0)))),float(rb['time'][keep][-1])
def main():
    analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    allrows=rows();lookup={r['name']:r for r in allrows};audits=[]
    for row in allrows:
        raw=np.load(DATA/(row['name']+'.npz'));t,y=raw['time'],raw['state'];m=mesh(row)
        assert np.isfinite(y).all() and np.all(np.diff(t)>0)
        for n,h in row['source_hashes'].items():resolve_source(n,h)
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
        audit=row['audit'];metrics=row['metrics'];complete=bool(row['failure'] is None and abs(t[-1]-row['requested_end'])<1e-9)
        budget=max(metrics[k] for k in ['material_budget_error','lorenz_budget_error','total_budget_error']) if metrics else None
        positive=audit['minimum_swept_clearance']>0 and audit['minimum_swept_self_clearance']>0 and audit['minimum_swept_area_ratio']>(.01 if row['old_material'] else m.p.area_floor) and audit['minimum_swept_quality_ratio']>(0 if row['old_material'] else m.p.shape_floor)
        passed=complete and positive and probe>=0 and crossings==0 and budget is not None and budget<1e-3 and metrics['maximum_wall_overlap']==0 and np.all(np.diff(y[:,-3:-1],axis=0)>=-1e-9)
        audits.append(dict(name=row['name'],completed=complete,passed=bool(passed),geometry_passed=bool(positive and probe>=0 and crossings==0),minimum_probe_gap=probe,proper_edge_crossings=crossings,maximum_budget_error=budget,
            actual_initial_lorenz=y[0,-7:-4].tolist(),sampled_peak_abs_X=float(np.max(abs(y[:,-7]))),sampled_peak_actuation=float(m.p.alpha/10*np.max(abs(y[:,-7]))),instantaneous_feedback_coefficient=float((m.p.alpha/10)**2*np.sum(m.mass)/(m.p.viscosity*m.p.lorenz_capacity)),
            raw_sha256=hashlib.sha256((DATA/(row['name']+'.npz')).read_bytes()).hexdigest()))
    pairs=[];blocks=[]
    for name,r in lookup.items():
        if not name.endswith('-one') or name[:-4]+'-two' not in lookup:continue
        b=lookup[name[:-4]+'-two']
        limited=bool(r['failure'] or b['failure']);end=min(r['completed_end'],b['completed_end'])
        if limited:end=np.floor((end+1e-10)/r['record_step'])*r['record_step']
        if end<=0:continue
        start=np.ceil((end/2-1e-10)/r['record_step'])*r['record_step'];a_stats=motion(r,start,end);b_stats=motion(b,start,end)
        changes={k:100*(b_stats[k]/a_stats[k]-1) if a_stats[k]>1e-14 else None for k in a_stats}
        gaps=[]
        for item in [r,b]:
            raw=np.load(DATA/(item['name']+'.npz'));idx=int(np.argmin(abs(raw['time']-end)))
            assert abs(raw['time'][idx]-end)<1e-10
            gaps.append(mesh(item).diagnostics(raw['state'][idx])['surface_gap_rms'])
        pairs.append(dict(key=name[:-4],start=start,end=end,limited_to_completed_common_window=limited,one=a_stats,two=b_stats,changes_percent=changes,
            one_gap=gaps[0],two_gap=gaps[1],
            forward=b['audit']['forward_work_step_trapezoid'],returned=b['audit']['return_work_step_trapezoid']))
        if r['group']=='long':
            for start in [4,10,16]:
                if start+6>end+1e-10:continue
                a_stats=motion(r,start,start+6);b_stats=motion(b,start,start+6)
                blocks.append(dict(key=name[:-4],start=start,end=start+6,changes_percent={k:100*(b_stats[k]/a_stats[k]-1) for k in a_stats}))
    checks=[]
    for label in ['prior-extreme-new','harder64']:
        for method in ['radau','bdf']:
            for direction in ['one','two']:
                a=label+'-'+direction;b=label+'-'+method+'-fine-'+direction
                if a in lookup and b in lookup and not lookup[a]['failure'] and not lookup[b]['failure']:
                    error=state_error(a,b);checks.append(dict(name=b,reference=a,error=error,target=1e-4,passed=error<1e-4))
                elif a in lookup and b in lookup and lookup[b]['completed_end']>0:
                    error,end=prefix_error(a,b);checks.append(dict(name=b+'-recorded-prefix',reference=a,end=end,error=error,target=1e-4,passed=error<1e-4,limited_to_recorded_prefix=True))
        a=label+'-bdf-fine-one';b=label+'-radau-fine-one'
        if a in lookup and b in lookup and min(lookup[a]['completed_end'],lookup[b]['completed_end'])>0:
            # Only exactly shared recording times; no dense interpolation is treated as a certificate.
            ra=np.load(DATA/(a+'.npz'));rb=np.load(DATA/(b+'.npz'));common=np.intersect1d(np.round(ra['time'],10),np.round(rb['time'],10))
            ia=[int(np.argmin(abs(ra['time']-t))) for t in common];ib=[int(np.argmin(abs(rb['time']-t))) for t in common]
            x=ra['state'][ia,:-4];y=rb['state'][ib,:-4];error=float(np.max(abs(x-y)/np.maximum(1,np.max(abs(y),axis=0))))
            checks.append(dict(name=label+'-tight-method-agreement',end=float(common[-1]),error=error,target=1e-4,passed=error<1e-4,limited_to_recorded_common_times=True))
    for direction in ['one','two']:
        a='old-method-'+direction;b='prior-extreme-'+direction
        if a in lookup and not lookup[a]['failure']:
            error=state_error(a,b);checks.append(dict(name='old-law-method-control-'+direction,error=error,target=1e-4,passed=error<1e-4))
    if all(n in lookup for n in ['zero-port-one','zero-port-two']):
        error=state_error('zero-port-one','zero-port-two');checks.append(dict(name='zero-port-control',error=error,target=1e-12,passed=error<1e-12))
    if 'old-method-fine-one' in lookup and not lookup['old-method-fine-one']['failure']:
        error=state_error('old-method-fine-one','prior-extreme-one');checks.append(dict(name='old-law-method-control-fine-one',reference='historical RK4',error=error,target=1e-4,passed=error<1e-4))
        if error<1e-4:
            for check in checks:
                if check['name']=='old-law-method-control-one':check['superseded_by']='old-law-method-control-fine-one'
    if 'sparse-agreement' in lookup and not lookup['sparse-agreement']['failure']:
        independent=np.load(DATA/'sparse-agreement.npz');normal=np.load(DATA/'prior-extreme-new-one.npz')
        ids=[int(np.argmin(abs(normal['time']-t))) for t in independent['time']]
        assert np.max(abs(normal['time'][ids]-independent['time']))<1e-10
        x=normal['state'][ids,:-4];y=independent['state'][:,:-4];error=float(np.max(abs(x-y)/np.maximum(1,np.max(abs(y),axis=0))))
        checks.append(dict(name='dense-sparse-agreement',end=.2,error=error,target=1e-4,passed=error<1e-4))
    for label,reference in [('prior','prior-extreme-new-one'),('hardest','harder64-one')]:
        name='numeric-jac-'+label
        if name in lookup and reference in lookup and not lookup[name]['failure']:
            independent=np.load(DATA/(name+'.npz'));normal=np.load(DATA/(reference+'.npz'))
            ids=[int(np.argmin(abs(normal['time']-t))) for t in independent['time']]
            assert np.max(abs(normal['time'][ids]-independent['time']))<1e-10
            x=normal['state'][ids,:-4];y=independent['state'][:,:-4];error=float(np.max(abs(x-y)/np.maximum(1,np.max(abs(y),axis=0))))
            checks.append(dict(name=name,reference=reference,end=.08,error=error,target=1e-4,passed=error<1e-4))
    work_checks=[]
    for label in ['prior-extreme-new','harder64']:
        for direction in ['one','two']:
            main=lookup.get(label+'-'+direction);fine=lookup.get(label+'-bdf-fine-'+direction)
            if main and fine and not main['failure'] and not fine['failure']:
                differences={k:abs(main['audit'][k]-fine['audit'][k])/max(1,abs(fine['audit'][k])) for k in ['forward_work_step_trapezoid','return_work_step_trapezoid']}
                work_checks.append(dict(name=label+'-'+direction,scaled_differences=differences))
    prior=[]
    for direction in ['one','two']:
        raw=np.load(DATA/('prior-extreme-'+direction+'.npz'));oldrow=json.loads((DATA/('prior-extreme-'+direction+'.json')).read_text())
        m=OldMesh(Parameters(**oldrow['parameters']));values=[geometry(m,m.unpack(s)[0]) for s in raw['state']]
        prior.append(dict(direction=direction,min_quality_ratio=min(v['min_quality_ratio'] for v in values),min_altitude_ratio=min(v['min_altitude_ratio'] for v in values),minimum_saved_area_ratio=min(v['min_area_ratio'] for v in values),minimum_swept_area_ratio=oldrow['audit']['minimum_swept_area_ratio']))
    mesh_checks=[]
    for coarse,fine in [('steady-mesh24','steady-mesh48'),('steady-mesh48','steady-mesh72')]:
        if all(n in lookup and not lookup[n]['failure'] for n in [coarse,fine]):
            x,y=[lookup[n]['metrics']['final'] for n in [coarse,fine]]
            mesh_checks.append(dict(coarse=coarse,fine=fine,end=6,coarse_observables={k:x[k] for k in ['imprint','surface_gap_rms']},fine_observables={k:y[k] for k in ['imprint','surface_gap_rms']},relative_differences={k:abs(x[k]-y[k])/max(1e-14,abs(y[k])) for k in ['imprint','surface_gap_rms']}))
    for direction in ['one','two']:
        coarse='prior-extreme-new-'+direction;fine='mesh48-hard-'+direction
        if all(n in lookup for n in [coarse,fine]):
            raw=np.load(DATA/(coarse+'.npz'));ref=np.load(DATA/(fine+'.npz'));common=np.intersect1d(np.round(raw['time'],10),np.round(ref['time'],10))
            end=float(common[-1]);idx=int(np.argmin(abs(raw['time']-end)));jdx=int(np.argmin(abs(ref['time']-end)))
            x=mesh(lookup[coarse]).diagnostics(raw['state'][idx]);y=mesh(lookup[fine]).diagnostics(ref['state'][jdx])
            mesh_checks.append(dict(coarse=coarse,fine=fine,end=end,limited_to_recorded_common_prefix=bool(lookup[fine]['failure']),coarse_observables={k:x[k] for k in ['imprint','surface_gap_rms']},fine_observables={k:y[k] for k in ['imprint','surface_gap_rms']},relative_differences={k:abs(x[k]-y[k])/max(1e-14,abs(y[k])) for k in ['imprint','surface_gap_rms']}))
    if all(n in lookup for n in ['old-form','steady-mesh24']):
        old=np.load(DATA/'old-form.npz')['state'][-1,:-7];new=np.load(DATA/'steady-mesh24.npz')['state'][-1,:-7]
        error=float(np.max(abs(old-new)/np.maximum(1,abs(old))))
        checks.append(dict(name='gentle-forming-material-preservation',end=6,error=error,target=1e-4,passed=error<1e-4,excludes_uncoupled_lorenz=True))
    protocol=json.loads((DATA/'protocol.json').read_text());upgrade=json.loads((DATA/'solver-upgrade-protocol.json').read_text())
    expected={c['name'] for c in protocol['configurations']+upgrade['configurations']}|set(protocol['memory_cases'])|{'zero-port-one','zero-port-two','numeric-jac-prior','numeric-jac-hardest'}
    out=dict(analysis_source_sha256=analysis_source_sha256,configurations=len(allrows),declared_configurations=len(expected),missing_recordings=sorted(expected-set(lookup)),unexpected_recordings=sorted(set(lookup)-expected),all_declared_attempts_recorded=expected<=set(lookup),completed=sum(a['completed'] for a in audits),audits=audits,pairs=pairs,long_blocks=blocks,solver_checks=checks,work_refinement=work_checks,mesh_refinement=mesh_checks,prior_baseline=prior,
        failures=[dict(name=r['name'],failure=r['failure'],completed_end=r['completed_end'],requested_end=r['requested_end']) for r in allrows if r['failure']],
        maximum_completed_budget_error=max((a['maximum_budget_error'] for a in audits if a['completed']),default=None),
        all_completed_audits_passed=all(a['passed'] for a in audits if a['completed']),all_solver_checks_passed=bool(checks) and all(c['passed'] for c in checks),current_solver_checks_passed=bool(checks) and all(c['passed'] for c in checks if 'superseded_by' not in c),
        notes='Finite windows and three starts provide observed variation, not confidence intervals, global stability or experimental calibration.')
    assert analysis_source_sha256==hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'Audit source changed during execution; rerun it.'
    (DATA/'all-runs.json').write_text(json.dumps(allrows,indent=2)+'\n');(DATA/'analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['audits','pairs','long_blocks']},indent=2),flush=True)
if __name__=='__main__':main()
