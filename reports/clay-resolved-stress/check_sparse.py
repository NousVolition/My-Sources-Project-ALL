"""Independent dense-vs-sparse Jacobian checks on recorded hard states."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,hashlib,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def main():
 old_h=module('old_h',ROOT/'legacy_source/analytic_spatial.py');old_j=module('old_j',ROOT/'legacy_source/jacobian.py');old_j.hessian=old_h.hessian
 checks=[];timing=[];rng=np.random.default_rng(715)
 for name in ['long48-one','long48-two','long72-one','long72-two']:
  path=ROOT/'parents'/(name+'.json');r=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));m=Mesh(Parameters(**r['parameters']))
  for index in [0,len(z['time'])-1]:
   y=z['state'][index];t=time.perf_counter();A=coupled_jacobian(m,m.rhs,0,y);sparse_seconds=time.perf_counter()-t;t=time.perf_counter();B=old_j.coupled_jacobian(m,m.rhs,0,y);dense_seconds=time.perf_counter()-t
   error=float(np.max(abs((A-B).toarray()))/max(1,np.max(abs(B.data))));assert error<2e-12
   d=rng.normal(size=m.size);d/=np.linalg.norm(d);eps=1e-7
   fd=(m.rhs(0,y+eps*d)-m.rhs(0,y-eps*d))/(2*eps);derivative=float(np.max(abs(fd-A@d))/max(1,np.max(abs(fd))))
   assert derivative<2e-4
   checks.append(dict(record=name,index=index,assembly_scaled_error=error,directional_error=derivative,passed=True));timing.append(dict(record=name,index=index,sparse_seconds=sparse_seconds,dense_seconds=dense_seconds))
 result=dict(all_passed=True,checks=checks,timing=timing,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (ROOT/'data/sparse-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
