"""Complete matched-interval comparisons; retain every attempted test."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,subprocess,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent
from audit_one import load
def compare(A,B,start=0,end=None):
 ra,_,ta,ya,da=A;rb,_,tb,yb,db=B;pairs=[(i,int(np.argmin(abs(tb-t)))) for i,t in enumerate(ta) if t>=start-1e-10 and (end is None or t<=end+1e-10) and min(abs(tb-t))<1e-10]
 if len(pairs)<2:return None
 ia,ib=np.array(pairs).T;out=dict(a=ra['name'],b=rb['name'],start=float(ta[ia[0]]),end=float(ta[ia[-1]]),samples=len(ia),coarse_mesh=[ra['parameters']['angles'],ra['parameters']['layers']],fine_mesh=[rb['parameters']['angles'],rb['parameters']['layers']])
 for key in ['gap_arclength_rms','surface_gap_rms','imprint']:
  a=np.array([da[i][key] for i in ia]);b=np.array([db[i][key] for i in ib]);err=abs(a-b)/np.maximum(.01,abs(b));j=int(np.argmax(err));out[key]=dict(maximum_scaled=float(max(err)),maximum_absolute=float(max(abs(a-b))),time_of_maximum=float(ta[ia[j]]),final_scaled=float(err[-1]))
 if ya.shape[1]==yb.shape[1]:out['state_scaled']=float(np.max(abs(ya[ia,:-4]-yb[ib,:-4])/np.maximum(1,np.max(abs(yb[ib,:-4]),axis=0))))
 out['lorenz_scaled']=float(np.max(abs(ya[ia,-7:-4]-yb[ib,-7:-4])/np.maximum(1,np.max(abs(yb[ib,-7:-4]),axis=0))))
 par=lambda r:{k:v for k,v in r['parameters'].items() if k not in ['angles','layers']}
 out['same_physics']=par(ra)==par(rb) and ra['initial']==rb['initial'];out['gap_target_met']=all(out[k]['maximum_scaled']<=.05 for k in ['gap_arclength_rms','surface_gap_rms']);out['same_requested_end']=ra['requested_end']==rb['requested_end']
 out['complete']=all(r['failure'] is None and abs(r['completed_end']-r['requested_end'])<1e-12 for r in [ra,rb]);out['refines_both_axes']=rb['parameters']['angles']>ra['parameters']['angles'] and rb['parameters']['layers']>ra['parameters']['layers']
 out['full_mesh_check_passed']=bool(out['same_physics'] and out['gap_target_met'] and out['complete'] and out['same_requested_end'] and out['refines_both_axes'] and out['start']<1e-10 and out['end']>=ra['requested_end']-1e-10 and (ra['parameters']['feedback'] or out['lorenz_scaled']<1e-4))
 return out
def main():
 paths=sorted([* (R/'baseline/data').glob('*.json'),*(R/'numerical/data').glob('*.json'),*(R/'solver/data').glob('*.json')]);rows=[];series={}
 cache=R/'data/audits';cache.mkdir(exist_ok=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();audit_sha=sha(R/'audit_one.py');source_sha=sha(Path(__file__))
 for p in paths:
  meta=json.loads(p.read_text())
  if 'parameters' not in meta:continue
  assert any(all((f/n).exists() and sha(f/n)==h for n,h in meta['source_hashes'].items()) for f in [R/v for v in ['baseline','baseline-fresh','numerical','solver']]),'Source mismatch'
  key=dict(audit_sha256=audit_sha,metadata_sha256=sha(p),raw_sha256=sha(p.with_suffix('.npz')));cp=cache/(meta['name']+'.json');cached=json.loads(cp.read_text()) if cp.exists() else None
  if cached and cached['key']==key:row=cached['result']
  else:
   v=subprocess.run([sys.executable,str(R/'audit_one.py'),p.relative_to(R).as_posix()],capture_output=True,text=True,check=True);row=json.loads(v.stdout);cp.write_text(json.dumps(dict(key=key,result=row),indent=2)+'\n')
  rows.append(row);series[row['name']]=load(p)
 for p in (R/'baseline-data').glob('*.json'):series[json.loads(p.read_text())['name']]=load(p)
 pairs=[('shared144-one','extend168-one'),('harder72-one','harder96-one'),('harder72-two','harder96-two'),('harder96-one','harder120-one'),('harder96-two','harder120-two')]
 for ell in ['050','100','200']:
  for mode in ['one','two']:
   for a,b in [(48,72),(72,96),(96,120),(120,144),(144,168),(168,192)]:pairs.append((f'length{ell}-{a}-{mode}',f'length{ell}-{b}-{mode}'))
 comparisons=[compare(series[a],series[b]) for a,b in pairs if a in series and b in series]
 controls=[]
 for n in [144,168]:
  for method in ['tight','radau']:
   a='extend168-one' if n==168 else 'shared144-one';b=f'peak{n}-{method}'
   if b in series:
    c=compare(series[a],series[b],start=series[b][0]['start_time'],end=series[b][0]['requested_end']);
    if c:c['passed']=bool(c['complete'] and c.get('state_scaled',1)<=1e-4);controls.append(c)
 if 'peak144-tight' in series and 'peak144-radau' in series:
  c=compare(series['peak144-tight'],series['peak144-radau'],start=.1,end=.14);c['passed']=bool(c['complete'] and c.get('state_scaled',1)<=1e-4);controls.append(c)
 if 'peak144-tight' in series and 'peak168-tight' in series:
  c=compare(series['peak144-tight'],series['peak168-tight'],start=.1,end=.14);c['purpose']='Spatial comparison after tightening the peak interval only';comparisons.append(c)
 for name,b in series.items():
  if name.endswith('-tight'):
   a=name[:-6]
   if a in series:
    c=compare(series[a],b,start=b[0]['start_time'],end=b[0]['requested_end']);c['passed']=bool(c['complete'] and c.get('state_scaled',1)<=1e-4);controls.append(c)
 if 'weak72-two-tight' in series and 'weak72-two-tighter' in series:
  c=compare(series['weak72-two-tight'],series['weak72-two-tighter']);c['passed']=bool(c['complete'] and c.get('state_scaled',1)<=1e-4);controls.append(c)
 for b in ['finest-weighted-radau','finest-original-radau']:
  if b in series:
   c=compare(series['length100-168-one'],series[b],start=series[b][0]['start_time'],end=series[b][0]['requested_end']);c['passed']=bool(c['complete'] and c.get('state_scaled',1)<=1e-4);controls.append(c)
 response=[]
 for n in [48,72,96,120,144,168,192]:
  a=f'length100-{n}-one';b=f'length100-{n}-two'
  if a not in series or b not in series:continue
  A,B=series[a],series[b];end=min(A[2][-1],B[2][-1]);metrics=[]
  for r,z,t,y,d in [A,B]:
   take=t<=end+1e-10;t=t[take];y=y[take];mean=np.trapezoid(y[:,-5],t)/end;std=np.sqrt(np.trapezoid((y[:,-5]-mean)**2,t)/end)
   metrics.append(dict(z_std=float(std),viscous_dissipation=float(y[-1,-3]),plastic_dissipation=float(y[-1,-2]),port_work=float(y[-1,-1]),final_X=float(y[-1,-7]),final_gap=float(d[len(t)-1]['gap_arclength_rms'])))
  response.append(dict(n=n,end=float(end),one=metrics[0],two=metrics[1],z_std_change_percent=100*(metrics[1]['z_std']/metrics[0]['z_std']-1)))
 sensitivities=[]
 for kind,x,y in [('material length','length050-72-one','length100-72-one'),('material length','length100-72-one','length200-72-one'),('initial perturbation 1e-8','weak72-two','weak72-perturb-two'),('changed constitutive model','original168-long-one','length100-168-one')]:
  if x in series and y in series:
   c=compare(series[x],series[y]);c['kind']=kind;sensitivities.append(c)
 for label,x,y in [('harder','harder72-one','harder72-two'),('intermediate capacity','middle72-one','middle72-two'),('larger capacity, loose','weak72-one','weak72-two'),('larger capacity, refined','weak72-one','weak72-two-tighter')]:
  if x not in series or y not in series:continue
  A,B=series[x],series[y];end=min(A[2][-1],B[2][-1]);metrics=[]
  for r,z,t,y,d in [A,B]:
   take=t<=end+1e-10;t=t[take];y=y[take];mean=np.trapezoid(y[:,-5],t)/end;std=np.sqrt(np.trapezoid((y[:,-5]-mean)**2,t)/end)
   metrics.append(dict(z_std=float(std),viscous_dissipation=float(y[-1,-3]),plastic_dissipation=float(y[-1,-2]),port_work=float(y[-1,-1]),final_X=float(y[-1,-7]),final_gap=float(d[len(t)-1]['gap_arclength_rms'])))
  response.append(dict(label=label,n=72,end=float(end),one=metrics[0],two=metrics[1],z_std_change_percent=100*(metrics[1]['z_std']/metrics[0]['z_std']-1)))
 assert source_sha==sha(Path(__file__)) and audit_sha==sha(R/'audit_one.py'),'Analysis source changed during audit'
 out=dict(attempts=len(rows),completed=sum(r['complete'] for r in rows),runs=rows,comparisons=comparisons,time_controls=controls,responses=response,sensitivities=sensitivities,all_geometry_passed=all(r['geometry_passed'] for r in rows),all_energy_passed=all(r['energy_passed'] for r in rows),analysis_sha256=source_sha,audit_one_sha256=audit_sha)
 for c in comparisons:
  participating=[r for r in rows if r['name'] in [c['a'],c['b']]]
  if participating:c['full_mesh_check_passed']=c['full_mesh_check_passed'] and all(r['geometry_passed'] and r['energy_passed'] for r in participating)
 (R/'analysis.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='runs'},indent=2))
if __name__=='__main__':main()
