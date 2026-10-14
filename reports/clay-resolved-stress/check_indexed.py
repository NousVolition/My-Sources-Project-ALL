"""Compare spatial-query forces against the original exhaustive implementation."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,hashlib,importlib.util,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
folder=ROOT/'solver-versions/indexed';sys.path.insert(0,str(folder))
spec=importlib.util.spec_from_file_location('indexed_mesh',folder/'mesh_hug.py');indexed=importlib.util.module_from_spec(spec);spec.loader.exec_module(indexed)
def main():
 checks=[];timing=[]
 for name,sub in [('long48-one','parents'),('long72-one','parents'),('extend72-one','solver/data'),('fine120-two','solver/data')]:
  path=ROOT/sub/(name+'.json');r=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));p=Parameters(**r['parameters']);a=Mesh(p);b=indexed.Mesh(p)
  ids=sorted(set([0,len(z['time'])-1,int(np.argmin([v['self_clearance'] for v in r['diagnostics']]))]))
  for i in ids:
   y=z['state'][i];old=a.rhs(0,y);new=b.rhs(0,y);err=float(np.max(abs(new-old))/max(1,np.max(abs(old))));assert err<2e-12
   checks.append(dict(record=name,index=i,error=err,passed=True))
  y=z['state'][-1];times=[]
  for m in [a,b]:
   t=time.perf_counter()
   for _ in range(20):m.rhs(0,y)
   times.append((time.perf_counter()-t)/20)
  timing.append(dict(record=name,old_rhs_seconds=times[0],indexed_rhs_seconds=times[1]))
 # Crossings whose endpoints are away from the contact range must still fail.
 m=indexed.Mesh(Parameters(angles=24,layers=3));x=m.initial.copy();x[0]=[-1.,-1.];x[1]=[1.,1.];x[4]=[-1.,1.];x[5]=[1.,-1.]
 try:m.self_contact_terms(x);crossing=False
 except FloatingPointError:crossing=True
 assert crossing;checks.append(dict(record='synthetic interior boundary crossing',passed=True))
 result=dict(all_passed=True,checks=checks,timing=timing,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(ROOT/'data/indexed-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
