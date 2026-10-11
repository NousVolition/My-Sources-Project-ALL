"""Audit recorded states and compare matched times without fitting trajectories."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'model'))
from mesh_hug import Mesh,Parameters
import importlib.util
spec=importlib.util.spec_from_file_location('production_hug',ROOT/'production/mesh_hug.py');production=importlib.util.module_from_spec(spec);sys.modules[spec.name]=production;spec.loader.exec_module(production)
def mesh_for(parameters):
    if 'self_neighborhood' in parameters:return production.Mesh(production.Parameters(**parameters))
    return Mesh(Parameters(**parameters))
def read(folder,name):
    r=json.loads((folder/(name+'.json')).read_text());a=np.load(folder/(name+'.npz'));return r,a
def compare(a,b):
    ta,ya=a['time'],a['state'];tb,yb=b['time'],b['state']
    pairs=[(i,int(np.argmin(abs(tb-t)))) for i,t in enumerate(ta) if np.min(abs(tb-t))<1e-10]
    assert len(pairs)>1
    ia,ib=np.array(pairs).T;aa=ya[ia,:-4];bb=yb[ib,:-4];scale=np.maximum(1,np.max(abs(bb),axis=0))
    e=float(np.max(abs(aa-bb)/scale));return dict(end=float(ta[ia[-1]]),scaled_state_error=e,passed=e<=1e-4)
def main():
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rows=[];raws={}
    for folder,prefix in [(ROOT/'data','pilot/'),(ROOT/'model/data',''),(ROOT/'refined/data','analytic/'),(ROOT/'production/data','production/'),(ROOT/'resolved/data','resolved/'),(ROOT/'solver/data','solver/')]:
        for f in sorted(folder.glob('*.json')):
            r=json.loads(f.read_text())
            if 'parameters' not in r or 'source_hashes' not in r:continue
            r,a=read(folder,f.stem);key=prefix+r['name'];m=mesh_for(r['parameters'])
            for name,sha in r['source_hashes'].items():assert hashlib.sha256((folder.parent/name).read_bytes()).hexdigest()==sha,(key,name)
            y=a['state'];t=a['time'];assert np.isfinite(y).all() and np.all(np.diff(t)>0)
            assert y.shape==(len(t),m.size) and len(r['diagnostics'])==len(t)
            assert np.array_equal(a['triangles'],m.tri) and np.array_equal(a['initial_positions'],m.initial)
            assert abs(t[-1]-r['completed_end'])<1e-12
            d=r['diagnostics'];mat=np.array([v['material_energy'] for v in d]);lor=np.array([v['lorenz_energy'] for v in d]);native,visc,plastic,port=y[:,-4:].T
            scale=lambda residual,E,W:float(np.max(abs(residual))/max(1,np.max(abs(E)),np.max(abs(W))))
            energy=dict(material=scale(mat-mat[0]+visc+plastic-port,mat,port),lorenz=scale(lor-lor[0]-native+(port if m.p.feedback else 0),lor,native),combined=r['energy_balance_scaled'])
            assert min(np.min(np.diff(visc)),np.min(np.diff(plastic)))>=-1e-9
            # All recorded meshes get a fresh crossing/clearance check.
            lo=1.;selflo=1.;walllo=1.;Jlo=1.
            for state in y:
                x=m.unpack(state)[0]
                lo=min(lo,float(np.min(m.clearance(x,False))))
                selflo=min(selflo,float(np.min(m.self_gap(x)[0]))) # crossing-domain guard included
                walllo=min(walllo,float(np.min(m.p.wall_half-abs(x))))
                Jlo=min(Jlo,float(np.min(m.areas(x)[0]/m.a0)))
            audit=r['audit'];safe=min(lo,selflo,walllo,audit['minimum_swept_clearance'],audit['minimum_swept_self_clearance'])>0 and audit['minimum_swept_area_ratio']>m.p.area_floor and audit['minimum_swept_quality_ratio']>m.p.shape_floor
            row=dict(key=key,name=r['name'],relative_path=f.relative_to(ROOT).as_posix(),parameters=r['parameters'],end=r['completed_end'],requested_end=r['requested_end'],failure=r['failure'],seconds=audit['seconds'],geometry_passed=bool(safe),energy_errors=energy,energy_passed=max(energy.values())<1e-3,minimum_area=audit['minimum_swept_area_ratio'],minimum_distortion=audit['minimum_swept_quality_ratio'],minimum_altitude=audit['minimum_swept_altitude_ratio'],final=d[-1],raw_sha256=hashlib.sha256(f.with_suffix('.npz').read_bytes()).hexdigest())
            rows.append(row);raws[key]=a
    by={r['key']:r for r in rows}
    mesh=[]
    for mode in ['one','two']:
        names=['pilot/revised24','mesh48-one','mesh72-one'] if mode=='one' else ['mesh24-two','mesh48-two','mesh72-two']
        for ca,fi in zip(names,names[1:]):
            if ca not in by or fi not in by:continue
            a,b=by[ca],by[fi]
            if abs(a['end']-.2)>1e-10 or abs(b['end']-.2)>1e-10:continue
            out=dict(mode=mode,coarse=ca,fine=fi,time=.2)
            for k in ['surface_gap_rms','gap_arclength_rms','imprint']:
                x,z=a['final'][k],b['final'][k];out[k]=dict(coarse=x,fine=z,absolute_difference=abs(x-z),relative_difference=abs(x-z)/max(1e-12,abs(z)))
            out['gap_passed']=out['gap_arclength_rms']['relative_difference']<=.05
            mesh.append(out)
    for mode in ['one','two']:
        for na,nb in [(24,48),(48,72),(72,96)]:
            ca=f'production/mesh{na}-{mode}';fi=f'production/mesh{nb}-{mode}'
            if ca not in by and f'resolved/mesh{na}-{mode}' in by:ca=f'resolved/mesh{na}-{mode}'
            if fi not in by and f'resolved/mesh{nb}-{mode}' in by:fi=f'resolved/mesh{nb}-{mode}'
            if ca not in by or fi not in by or min(by[ca]['end'],by[fi]['end'])<.2-1e-10:continue
            a,b=by[ca],by[fi];out=dict(mode='production '+mode,coarse=ca,fine=fi,time=.2)
            for k in ['surface_gap_rms','gap_arclength_rms','imprint']:
                x,z=a['final'][k],b['final'][k];out[k]=dict(coarse=x,fine=z,absolute_difference=abs(x-z),relative_difference=abs(x-z)/max(1e-12,abs(z)))
            out['gap_passed']=out['gap_arclength_rms']['relative_difference']<=.05;mesh.append(out)
    for group,end in [('mesh',.2),('hard',.1)]:
        for mode in ['one','two']:
            for na,nb in ([(24,48),(48,72)] if group=='mesh' else [(48,72)]):
                ca=f'solver/{group}{na}-{mode}';fi=f'solver/{group}{nb}-{mode}'
                if ca not in by or fi not in by or min(by[ca]['end'],by[fi]['end'])<end-1e-10:continue
                aa,bb=by[ca],by[fi];out=dict(mode='final '+group+' '+mode,coarse=ca,fine=fi,time=end)
                for k in ['surface_gap_rms','gap_arclength_rms','imprint']:
                    x,z=aa['final'][k],bb['final'][k];out[k]=dict(coarse=x,fine=z,absolute_difference=abs(x-z),relative_difference=abs(x-z)/max(1e-12,abs(z)))
                out['gap_passed']=out['gap_arclength_rms']['relative_difference']<=.05;mesh.append(out)
    comparisons=[]
    for a,b in [('pilot/revised48','mesh48-one'),('mesh48-one','analytic/analytic48-one'),('stress16-one','accurate16-one'),('stress16-two','accurate16-two'),('zero-one','zero-two'),('pilot/revised24','quad16-one'),('mesh24-two','quad16-two'),('quad16-one','analytic/quad32-one'),('resolved/stress16-one','resolved/accurate16-one'),('resolved/stress16-two','resolved/accurate16-two'),('production/mesh24-one','resolved/quad64-one'),('resolved/zero-one','resolved/zero-two'),('solver/stress16-one','solver/accurate16-one'),('solver/stress16-two','solver/accurate16-two'),('solver/mesh24-one','solver/quad64-one'),('solver/zero-one','solver/zero-two')]:
        if a in raws and b in raws:
            c=compare(raws[a],raws[b]);c.update(a=a,b=b,full_requested_interval=bool(min(by[a]['end'],by[b]['end'])>=min(by[a]['requested_end'],by[b]['requested_end'])-1e-10));comparisons.append(c)
    pairs=[]
    for name in [prefix+n for prefix in ['','resolved/','solver/'] for n in ['stress16','stress32','stress64','stress128','window0','window1','window2']]:
        one=name+'-one';two=name+'-two'
        if one not in by and 'analytic/'+one in by:one='analytic/'+one
        if one not in by or two not in by:continue
        common=min(by[one]['end'],by[two]['end']);start=2. if 'window' in name else 0.
        end=np.floor((common+1e-10)*100)/100
        if end<=start:continue
        metrics=[]
        for key in [one,two]:
            a=raws[key];tt=a['time'];y=a['state'];i=int(np.argmin(abs(tt-start)));j=int(np.argmin(abs(tt-end)));assert max(abs(tt[i]-start),abs(tt[j]-end))<1e-8
            t=tt[i:j+1];z=y[i:j+1,-5];mean=np.trapezoid(z,t)/(end-start);std=np.sqrt(np.trapezoid((z-mean)**2,t)/(end-start));m=mesh_for(by[key]['parameters']);p=m.p
            motion=np.sqrt((y[j,-3]-y[i,-3])/(p.viscosity*sum(m.mass)*(end-start)))
            metrics.append(dict(motion_rms=float(motion),z_std=float(std),net_port=float(y[j,-1]-y[i,-1])))
        pairs.append(dict(name=name,start=start,end=end,one=metrics[0],two=metrics[1],motion_change_percent=100*(metrics[1]['motion_rms']/metrics[0]['motion_rms']-1),z_std_change_percent=100*(metrics[1]['z_std']/metrics[0]['z_std']-1)))
    assert source_hash==hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'Analysis source changed while auditing'
    result=dict(runs=rows,mesh_comparisons=mesh,solver_comparisons=comparisons,pairs=pairs,completed=sum(r['failure'] is None for r in rows),attempted=len(rows),all_geometry_passed=all(r['geometry_passed'] for r in rows),all_energy_passed=all(r['energy_passed'] for r in rows),analysis_sha256=source_hash)
    (ROOT/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='runs'},indent=2));assert result['all_geometry_passed'] and result['all_energy_passed']
if __name__=='__main__':main()
