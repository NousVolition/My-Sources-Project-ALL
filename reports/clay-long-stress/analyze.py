"""Audit continuations, source identity, recorded geometry and matched intervals."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
from contact_hug import Mesh as Exhaustive
from adaptive_mesh import select

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def source_check(record,folder):
    for name,h in record['source_hashes'].items():assert sha(folder/name)==h,(name,'source mismatch')

def load(path):
    r=read(path);raw=np.load(path.with_suffix('.npz'));t=raw['time'];y=raw['state'];d=list(r['diagnostics'])
    source_check(r,ROOT/('legacy_source' if path.parent.name=='parents' else 'solver'))
    assert len(d)==len(t)==len(y) and np.isfinite(y).all() and np.all(np.diff(t)>0)
    assert abs(t[-1]-r['completed_end'])<1e-12
    if r.get('parent'):
        p=ROOT/r['parent']['metadata'];pr=read(p);z=np.load(p.with_suffix('.npz'));source_check(pr,ROOT/'legacy_source')
        assert sha(p)==r['parent']['metadata_sha256'] and sha(p.with_suffix('.npz'))==r['parent']['raw_sha256']
        assert r['parameters']==pr['parameters'] and r['initial']==pr['initial']
        assert np.array_equal(z['state'][-1],y[0]) and z['time'][-1]==t[0]==r['start_time']
        t=np.r_[z['time'],t[1:]];y=np.concatenate([z['state'],y[1:]]);d=pr['diagnostics']+d[1:]
    return r,raw,t,y,d

def matched(a,b,start=0.):
    ta,ya,da=a;tb,yb,db=b
    pairs=[(i,int(np.argmin(abs(tb-t)))) for i,t in enumerate(ta) if t>=start-1e-10 and np.min(abs(tb-t))<1e-10]
    ia,ib=np.array(pairs).T
    result=dict(start=float(ta[ia[0]]),end=float(ta[ia[-1]]),samples=len(ia),independent_interval=len(ia)>1)
    for key in ['gap_arclength_rms','surface_gap_rms','imprint']:
        x=np.array([da[i][key] for i in ia]);z=np.array([db[i][key] for i in ib]);scaled=abs(x-z)/np.maximum(.01,abs(z))
        result[key]=dict(maximum_scaled=float(max(scaled)),final_scaled=float(scaled[-1]),coarse_final=float(x[-1]),fine_final=float(z[-1]))
    if ya.shape[1]==yb.shape[1]:result['state_scaled']=float(np.max(abs(ya[ia,:-4]-yb[ib,:-4])/np.maximum(1,np.max(abs(yb[ib,:-4]),axis=0))))
    return result

def main():
    source=sha(__file__);rows=[];series={}
    for path in sorted((ROOT/'solver/data').glob('*.json')):
        if 'parameters' not in read(path):continue
        r,raw,t,y,d=load(path);m=Mesh(Parameters(**r['parameters']))
        assert np.array_equal(m.initial,raw['initial_positions']) and np.array_equal(m.tri,raw['triangles']) and y.shape[1]==m.size
        object_min=wall_min=self_min=1.;area_min=quality_min=1.;energy=[]
        for state,saved in zip(y,d):
            x,rest,l=m.unpack(state);e=m.energy_force(x,rest)[0]
            assert abs(e-saved['material_energy'])<1e-9*max(1,abs(e))
            energy.append(e);object_min=min(object_min,float(min(m.clearance(x,False))));self_min=min(self_min,float(min(Exhaustive.self_gap(m,x)[0])))
            wall_min=min(wall_min,float(np.min(m.p.wall_half-abs(x))))
            J=m.areas(x)[0]/m.a0;Q=2*J/np.sum(m.deformation(x)**2,axis=(1,2))
            area_min=min(area_min,float(min(J)));quality_min=min(quality_min,float(min(Q)))
        E=np.array(energy);L=np.array([v['lorenz_energy'] for v in d]);native,visc,plastic,port=y[:,-4:].T
        assert min(np.min(np.diff(visc)),np.min(np.diff(plastic)))>-1e-8
        scale=lambda z,e,w:float(max(abs(z))/max(1,max(abs(e)),max(abs(w))))
        errors=dict(material=scale(E-E[0]+visc-visc[0]+plastic-plastic[0]-port+port[0],E,port-port[0]),lorenz=scale(L-L[0]-native+native[0]+(port-port[0] if m.p.feedback else 0),L,native-native[0]))
        combined=E+L;work=native-visc-plastic+(0 if m.p.feedback else port)
        errors['combined']=scale(combined-combined[0]-(work-work[0]),combined,work-work[0])
        audit=r['audit'];geometry=min(object_min,wall_min,self_min,audit['minimum_swept_clearance'],audit['minimum_swept_self_clearance'])>0 and min(area_min,audit['minimum_swept_area_ratio'])>m.p.area_floor and min(quality_min,audit['minimum_swept_quality_ratio'])>m.p.shape_floor
        row=dict(name=r['name'],path=path.relative_to(ROOT).as_posix(),parameters=r['parameters'],start=r.get('start_time',0),end=r['completed_end'],requested_end=r['requested_end'],failure=r['failure'],seconds=audit['seconds'],steps=audit['accepted_steps'],minimum_object=object_min,minimum_wall=wall_min,minimum_self=self_min,minimum_area=min(area_min,audit['minimum_swept_area_ratio']),minimum_quality=min(quality_min,audit['minimum_swept_quality_ratio']),minimum_altitude=audit['minimum_swept_altitude_ratio'],geometry_passed=bool(geometry),energy_errors=errors,energy_passed=max(errors.values())<1e-3,raw_sha256=sha(path.with_suffix('.npz')),metadata_sha256=sha(path),final=d[-1])
        rows.append(row);series[r['name']]=(t,y,d)
    by={r['name']:r for r in rows}
    for path in (ROOT/'parents').glob('*.json'):
        r,z,t,y,d=load(path);series['parent/'+r['name']]=(t,y,d)
    comparisons=[]
    for a,b in [('fresh24-one','parent/mesh24-one'),('fresh24-two','parent/mesh24-two'),('long48-one','long72-one'),('long48-two','long72-two'),('long24-one','radau24-one'),('long24-two','radau24-two'),('long48-one','radau48-probe')]:
        if a not in series or b not in series:continue
        result=matched(series[a],series[b]);result.update(a=a,b=b)
        result['both_completed']=all(name.startswith('parent/') or by[name]['failure'] is None for name in [a,b])
        result['same_requested_interval']=a in by and b in by and by[a]['requested_end']==by[b]['requested_end']
        result['gap_target_met_on_shared_times']=all(result[k]['maximum_scaled']<=.05 for k in ['gap_arclength_rms','surface_gap_rms'])
        result['full_mesh_check_passed']=result['both_completed'] and result['same_requested_interval'] and result['gap_target_met_on_shared_times']
        comparisons.append(result)
    time_controls=[]
    for a,b in [('long24-one','radau24-one'),('long24-two','radau24-two'),('long48-one','radau48-probe')]:
        if a not in series or b not in series:continue
        start=max(by[a]['start'],by[b]['start']);c=matched(series[a],series[b],start)
        c.update(a=a,b=b,requested_control_end=by[b]['requested_end'])
        c['passed']=bool(c['independent_interval'] and by[b]['failure'] is None and c['end']>=by[b]['requested_end']-1e-10 and c['state_scaled']<=1e-4)
        time_controls.append(c)
    axes=[]
    for name in ['axis24x5','axis24x7','axis48x3']:
        if name in series:axes.append(dict(name=name,**matched(series[name],series['parent/mesh72-one'])))
    pairs=[]
    for n in [24,48,72]:
        aa=f'long{n}-one';bb=f'long{n}-two'
        if aa not in series or bb not in series:continue
        end=np.floor((min(by[aa]['end'],by[bb]['end'])+1e-10)*100)/100
        for start in [0.,.5 if n==24 else .1]:
            if end<=start:continue
            metrics=[]
            for name in [aa,bb]:
                t,y,d=series[name];i=int(np.argmin(abs(t-start)));j=int(np.argmin(abs(t-end)));assert abs(t[i]-start)<1e-9 and abs(t[j]-end)<1e-9
                ts=t[i:j+1];Z=y[i:j+1,-5];mean=np.trapezoid(Z,ts)/(end-start);std=np.sqrt(np.trapezoid((Z-mean)**2,ts)/(end-start))
                p=Parameters(**by[name]['parameters']);m=Mesh(p);motion=np.sqrt((y[j,-3]-y[i,-3])/(p.viscosity*sum(m.mass)*(end-start)))
                metrics.append(dict(motion_rms=float(motion),z_std=float(std),net_port=float(y[j,-1]-y[i,-1])))
            pairs.append(dict(mesh=n,start=start,end=end,one=metrics[0],two=metrics[1],motion_change_percent=100*(metrics[1]['motion_rms']/metrics[0]['motion_rms']-1),z_std_change_percent=100*(metrics[1]['z_std']/metrics[0]['z_std']-1)))
    decisions={mode:select([ROOT/'parents'/f'mesh{n}-{mode}.json' for n in [24,48,72]]) for mode in ['one','two']}
    # Run the same selection using the new independently executed coarse results.
    for mode in ['one','two']:
        if f'fresh24-{mode}' in by:decisions['fresh-'+mode]=select([ROOT/'solver/data'/f'fresh24-{mode}.json',*[ROOT/'parents'/f'mesh{n}-{mode}.json' for n in [48,72]]])
    result=dict(attempts=len(rows),completed=sum(r['failure'] is None for r in rows),runs=rows,comparisons=comparisons,time_controls=time_controls,axis_probes=axes,pairs=pairs,resolution=decisions,all_geometry_passed=all(r['geometry_passed'] for r in rows),all_energy_passed=all(r['energy_passed'] for r in rows),analysis_sha256=source,adaptive_sha256=sha(ROOT/'adaptive_mesh.py'))
    assert sha(__file__)==source
    (ROOT/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['runs','resolution']},indent=2));assert result['all_geometry_passed'] and result['all_energy_passed']
if __name__=='__main__':main()
