"""Independent raw-state, provenance, geometry, energy and convergence audit."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
from contact_hug import Mesh as Exhaustive
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def source_check(r):
 folders=[ROOT/'solver',ROOT/'legacy_source',ROOT/'checkpoint_source',*sorted((ROOT/'solver-versions').glob('*'))]
 assert any(all((f/n).exists() and sha(f/n)==h for n,h in r['source_hashes'].items()) for f in folders),'source mismatch'
def load(p,join=True):
 p=Path(p);r=read(p);source_check(r);z=np.load(p.with_suffix('.npz'));t=z['time'];y=z['state'];d=r['diagnostics']
 assert len(t)==len(y)==len(d) and np.isfinite(y).all() and np.all(np.diff(t)>0) and abs(t[-1]-r['completed_end'])<1e-12
 if r.get('parent'):
  pp=ROOT/r['parent']['metadata'];assert sha(pp)==r['parent']['metadata_sha256'] and sha(pp.with_suffix('.npz'))==r['parent']['raw_sha256']
  pr,pz,pt,py,pd=load(pp)
  assert r['parameters']==pr['parameters'] and r['initial']==pr['initial'] and np.array_equal(y[0],py[-1]) and t[0]==pt[-1]==r['start_time']
  if join:t=np.r_[pt,t[1:]];y=np.concatenate([py,y[1:]]);d=pd+d[1:]
 return r,z,t,y,d
def compare(a,b,start=0.,end=None):
 ra,_,ta,ya,da=a;rb,_,tb,yb,db=b
 pairs=[(i,int(np.argmin(abs(tb-t)))) for i,t in enumerate(ta) if t>=start-1e-10 and (end is None or t<=end+1e-10) and min(abs(tb-t))<1e-10]
 assert len(pairs)>1
 ia,ib=np.array(pairs).T;out=dict(a=ra['name'],b=rb['name'],start=float(ta[ia[0]]),end=float(ta[ia[-1]]),samples=len(pairs),coarse_mesh=[ra['parameters']['angles'],ra['parameters']['layers']],fine_mesh=[rb['parameters']['angles'],rb['parameters']['layers']])
 for key in ['gap_arclength_rms','surface_gap_rms','imprint']:
  x=np.array([da[i][key] for i in ia]);z=np.array([db[j][key] for j in ib]);err=abs(x-z)/np.maximum(.01,abs(z));idx=int(np.argmax(err))
  out[key]=dict(maximum_scaled=float(max(err)),maximum_absolute=float(max(abs(x-z))),time_of_maximum=float(ta[ia[idx]]),final_scaled=float(err[-1]),coarse_final=float(x[-1]),fine_final=float(z[-1]))
 if ya.shape[1]==yb.shape[1]:out['state_scaled']=float(np.max(abs(ya[ia,:-4]-yb[ib,:-4])/np.maximum(1,np.max(abs(yb[ib,:-4]),axis=0))))
 out['lorenz_scaled']=float(np.max(abs(ya[ia,-7:-4]-yb[ib,-7:-4])/np.maximum(1,np.max(abs(yb[ib,-7:-4]),axis=0))))
 out['drive_agreement_passed']=bool(ra['parameters']['feedback'] or out['lorenz_scaled']<=1e-4)
 same={k:v for k,v in ra['parameters'].items() if k not in ['angles','layers']}=={k:v for k,v in rb['parameters'].items() if k not in ['angles','layers']} and ra['initial']==rb['initial']
 out['same_physics']=same;out['gap_target_met']=all(out[k]['maximum_scaled']<=.05 for k in ['gap_arclength_rms','surface_gap_rms'])
 out['complete']=bool(not ra['failure'] and not rb['failure'] and abs(ra['completed_end']-ra['requested_end'])<1e-10 and abs(rb['completed_end']-rb['requested_end'])<1e-10)
 out['same_requested_end']=ra['requested_end']==rb['requested_end']
 out['refines_both_axes']=rb['parameters']['angles']>ra['parameters']['angles'] and rb['parameters']['layers']>ra['parameters']['layers']
 out['geometry_energy_passed']=all(r['energy_balance_scaled']<1e-3 and r['audit']['minimum_swept_area_ratio']>r['parameters']['area_floor'] and r['audit']['minimum_swept_quality_ratio']>r['parameters']['shape_floor'] and min(r['audit']['minimum_swept_clearance'],r['audit']['minimum_swept_self_clearance'])>0 for r in [ra,rb])
 out['full_mesh_check_passed']=bool(same and out['complete'] and out['refines_both_axes'] and out['geometry_energy_passed'] and out['drive_agreement_passed'] and out['same_requested_end'] and out['gap_target_met'] and out['start']<1e-10 and out['end']>=ra['requested_end']-1e-10)
 return out
def main():
 rows=[];series={};source=sha(__file__)
 paths=[*(ROOT/'solver/data').glob('*.json'),*(ROOT/'solver-versions').glob('*/data/*.json')]
 for path in sorted(paths):
  if 'parameters' not in read(path):continue
  r,z,t,y,d=load(path);m=Mesh(Parameters(**r['parameters']));assert np.array_equal(m.initial,z['initial_positions']) and np.array_equal(m.tri,z['triangles']) and y.shape[1]==m.size
  object_min=wall_min=self_min=area_min=quality_min=1.;E=[]
  for state,v in zip(y,d):
   x,rest,l=m.unpack(state);e=m.energy_force(x,rest)[0];assert abs(e-v['material_energy'])<1e-8*max(1,abs(e));E.append(e)
   object_min=min(object_min,float(min(m.clearance(x,False))));self_min=min(self_min,float(min(Exhaustive.self_gap(m,x)[0])));wall_min=min(wall_min,float(np.min(m.p.wall_half-abs(x))))
   J=m.areas(x)[0]/m.a0;Q=2*J/np.sum(m.deformation(x)**2,axis=(1,2));area_min=min(area_min,float(min(J)));quality_min=min(quality_min,float(min(Q)))
   assert abs(m.diagnostics(state)['gap_arclength_rms']-v['gap_arclength_rms'])<1e-12
  E=np.array(E);L=np.array([v['lorenz_energy'] for v in d]);native,visc,plastic,port=y[:,-4:].T
  assert min(np.min(np.diff(visc)),np.min(np.diff(plastic)))>-1e-6
  scale=lambda z,e,w:float(max(abs(z))/max(1,max(abs(e)),max(abs(w))))
  errors=dict(material=scale(E-E[0]+visc-visc[0]+plastic-plastic[0]-port+port[0],E,port-port[0]),lorenz=scale(L-L[0]-native+native[0]+(port-port[0] if m.p.feedback else 0),L,native-native[0]))
  energy=E+L;work=native-visc-plastic+(0 if m.p.feedback else port);errors['combined']=scale(energy-energy[0]-(work-work[0]),energy,work-work[0])
  audit=r['audit'];geom=min(object_min,wall_min,self_min,audit['minimum_swept_clearance'],audit['minimum_swept_self_clearance'])>0 and min(area_min,audit['minimum_swept_area_ratio'])>m.p.area_floor and min(quality_min,audit['minimum_swept_quality_ratio'])>m.p.shape_floor
  rows.append(dict(name=r['name'],path=path.relative_to(ROOT).as_posix(),parameters=r['parameters'],start=r['start_time'],end=r['completed_end'],requested_end=r['requested_end'],failure=r['failure'],seconds=audit['seconds'],steps=audit['accepted_steps'],minimum_object=object_min,minimum_wall=wall_min,minimum_self=self_min,minimum_area=min(area_min,audit['minimum_swept_area_ratio']),minimum_quality=min(quality_min,audit['minimum_swept_quality_ratio']),minimum_altitude=audit['minimum_swept_altitude_ratio'],energy_errors=errors,energy_passed=max(errors.values())<1e-3,geometry_passed=bool(geom),raw_sha256=sha(path.with_suffix('.npz')),metadata_sha256=sha(path),final=d[-1]))
  series[r['name']]=(r,z,t,y,d)
 for path in (ROOT/'parents').glob('*.json'):
  a=load(path);series['parent/'+a[0]['name']]=a
 comps=[]
 for a,b in [('parent/long48-one','parent/long72-one'),('extend72-one','fine96-one'),('parent/long72-two','fine96-two'),('fine96-one','fine120-one'),('fine96-two','fine120-two'),('extend72-two','extend120-two'),('extend120-two','extend144-two'),('extend144-two','shared168-long-two'),('fine120-one','fine144-one'),('indexed120-one','fine144-one'),('fine120-one','extend144-one'),('indexed120-one','extend144-one'),('fine120-two','fine144-two'),('shared144-one','shared168-one'),('shared144-one','extend168-one'),('fine144-two','shared168-two')]:
  if a in series and b in series:comps.append(compare(series[a],series[b]))
 controls=[]
 for a,b in [('extend72-one','control72-radau'),('fine96-one','control96-tight'),('fine96-two','control96-two-tight'),('fine120-one','indexed120-one'),('fine144-one','shared144-one'),('extend144-one','shared144-one'),('shared168-one','control168-shared-radau'),('shared144-one','resume-check-144-one'),('extend168-one','control168-resume-tight')]:
  if a in series and b in series:
   c=compare(series[a],series[b],start=series[b][0]['start_time'],end=series[b][0]['requested_end']);c['passed']=bool(not series[b][0]['failure'] and c['end']>=series[b][0]['requested_end']-1e-10 and c.get('state_scaled',1)<=1e-4);controls.append(c)
 pairs=[]
 for a,b in [('extend72-one','extend72-two'),('fine96-one','fine96-two'),('fine120-one','fine120-two'),('extend144-one' if 'extend144-one' in series else 'fine144-one','fine144-two'),('extend168-one' if 'extend168-one' in series else 'shared168-one','shared168-two')]:
  if a not in series or b not in series:continue
  A=series[a];B=series[b];end=np.floor(min(A[2][-1],B[2][-1])*200+1e-8)/200;metrics=[]
  for s in [A,B]:
   r,z,t,y,d=s;sel=(t<=end+1e-10);t=t[sel];y=y[sel];duration=t[-1]-t[0];Z=y[:,-5];mean=np.trapezoid(Z,t)/duration;std=np.sqrt(np.trapezoid((Z-mean)**2,t)/duration);m=Mesh(Parameters(**r['parameters']));motion=np.sqrt((y[-1,-3]-y[0,-3])/(m.p.viscosity*sum(m.mass)*duration));metrics.append(dict(motion_rms=float(motion),z_std=float(std),end=float(t[-1]),net_port=float(y[-1,-1]-y[0,-1])))
  assert abs(metrics[0]['end']-metrics[1]['end'])<1e-9
  pairs.append(dict(mesh=A[0]['parameters']['angles'],end=end,one=metrics[0],two=metrics[1],motion_change_percent=100*(metrics[1]['motion_rms']/metrics[0]['motion_rms']-1),z_std_change_percent=100*(metrics[1]['z_std']/metrics[0]['z_std']-1)))
 result=dict(attempts=len(rows),completed=sum(r['failure'] is None for r in rows),runs=rows,comparisons=comps,time_controls=controls,pairs=pairs,all_geometry_passed=all(r['geometry_passed'] for r in rows),all_energy_passed=all(r['energy_passed'] for r in rows),analysis_sha256=source)
 assert sha(__file__)==source
 (ROOT/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='runs'},indent=2));assert result['all_geometry_passed'] and result['all_energy_passed']
if __name__=='__main__':main()
