"""Recompute saved-state geometry, gap and separate energy balances."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def resolve(h):
 found=[p for p in R.rglob('*.json') if p.parent.name not in ['__pycache__'] and sha(p)==h]
 assert found,('Missing parent metadata',h)
 return found[0]
def source_check(r):
 assert any(all((f/n).exists() and sha(f/n)==h for n,h in r['source_hashes'].items()) for f in [R/v for v in ['baseline','baseline-fresh','numerical','solver']]),'Source mismatch in recording or parent'
def certificates(p):
 r=json.loads(p.read_text());source_check(r);out=[r['audit']]
 if r.get('parent'):out+=certificates(resolve(r['parent']['metadata_sha256']))
 return out
def load(p,join=True):
 r=json.loads(p.read_text());source_check(r);z=np.load(p.with_suffix('.npz'));t=z['time'];y=z['state'];d=r['diagnostics']
 assert y.shape[0]==len(t)==len(d) and np.isfinite(y).all() and np.all(np.diff(t)>0)
 assert abs(t[0]-r['start_time'])<1e-12 and abs(t[-1]-r['completed_end'])<1e-12
 if r.get('parent'):
  pp=resolve(r['parent']['metadata_sha256']);assert sha(pp.with_suffix('.npz'))==r['parent']['raw_sha256'];pr,pz,pt,py,pd=load(pp)
  assert r['parameters']==pr['parameters'] and r['initial']==pr['initial'] and t[0]==pt[-1] and np.array_equal(y[0],py[-1])
  if join:t=np.r_[pt,t[1:]];y=np.concatenate([py,y[1:]]);d=pd+d[1:]
 return r,z,t,y,d
def main():
 p=(R/sys.argv[1]).resolve();r,z,t,y,d=load(p);folders=[R/a for a in ['baseline','baseline-fresh','numerical','solver']]
 candidates=[f for f in folders if all((f/n).exists() and sha(f/n)==h for n,h in r['source_hashes'].items())];assert len(candidates)>=1,'Source mismatch'
 folder=candidates[0];sys.path.insert(0,str(folder));from mesh_hug import Mesh,Parameters
 from contact_hug import Mesh as Exhaustive
 m=Mesh(Parameters(**r['parameters']));assert np.array_equal(m.tri,z['triangles']) and np.array_equal(m.initial,z['initial_positions'])
 E=[];L=[];mins=dict(object=1.,wall=1.,self=1.,area=1.,quality=1.);max_gap_recompute=0
 for state,v in zip(y,d):
  x,rest,l=m.unpack(state);e=m.energy_force(x,rest)[0];E.append(e);L.append(.5*m.p.lorenz_capacity*np.dot(l,l))
  assert abs(e-v['material_energy'])<1e-8*max(1,abs(e))
  g=m.diagnostics(state);max_gap_recompute=max(max_gap_recompute,abs(g['gap_arclength_rms']-v['gap_arclength_rms']),abs(g['surface_gap_rms']-v['surface_gap_rms']))
  vals=dict(object=min(m.clearance(x,False)),wall=np.min(m.p.wall_half-abs(x)),self=min(Exhaustive.self_gap(m,x)[0]),area=min(m.areas(x)[0]/m.a0),quality=min(2*m.areas(x)[0]/m.a0/np.sum(m.deformation(x)**2,axis=(1,2))))
  for k,v in vals.items():mins[k]=min(mins[k],float(v))
 E=np.array(E);L=np.array(L);native,visc,plastic,port=y[:,-4:].T
 def error(z,e,w):return float(max(abs(z))/max(1,max(abs(e)),max(abs(w))))
 errs=dict(material=error(E-E[0]+visc-visc[0]+plastic-plastic[0]-port+port[0],E,port-port[0]),lorenz=error(L-L[0]-native+native[0]+(port-port[0] if m.p.feedback else 0),L,native-native[0]))
 work=native-visc-plastic+(0 if m.p.feedback else port);errs['combined']=error(E+L-E[0]-L[0]-(work-work[0]),E+L,work-work[0])
 a=r['audit'];chain=certificates(p);cert={k:min(v[k] for v in chain) for k in ['minimum_swept_clearance','minimum_swept_self_clearance','minimum_swept_area_ratio','minimum_swept_quality_ratio','minimum_swept_altitude_ratio']};geometry=mins['object']>0 and mins['wall']>0 and mins['self']>0 and min(mins['area'],cert['minimum_swept_area_ratio'])>m.p.area_floor and min(mins['quality'],cert['minimum_swept_quality_ratio'])>m.p.shape_floor and min(cert['minimum_swept_clearance'],cert['minimum_swept_self_clearance'])>0
 monotonic=min(np.min(np.diff(visc)),np.min(np.diff(plastic)))>=-1e-6
 out=dict(name=r['name'],path=p.relative_to(R).as_posix(),source=folder.name,parameters=r['parameters'],initial=r['initial'],start=r['start_time'],end=r['completed_end'],requested_end=r['requested_end'],failure=r['failure'],complete=r['failure'] is None and abs(r['completed_end']-r['requested_end'])<1e-12,samples=len(t),seconds=a['seconds'],steps=a['accepted_steps'],minimum=mins,minimum_swept_altitude_ratio=a['minimum_swept_altitude_ratio'],energy_errors=errs,geometry_passed=bool(geometry),energy_passed=bool(max(errs.values())<1e-3 and monotonic),gap_recompute_max=max_gap_recompute,raw_sha256=sha(p.with_suffix('.npz')),metadata_sha256=sha(p),final=d[-1])
 assert max_gap_recompute<1e-12
 out.update(certificate_minima=cert,certificate_sources=len(chain),minimum_swept_altitude_ratio=cert['minimum_swept_altitude_ratio'],certificate_scope='Current interval and all inherited accepted-path bounds; derived prefix bounds may conservatively cover the original full source')
 print(json.dumps(out,allow_nan=False))
if __name__=='__main__':main()
