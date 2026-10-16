"""Independent full sparse solves of real and complex implicit stage matrices."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,time,hashlib,numpy as np
from pathlib import Path
from scipy.sparse import eye
from scipy.sparse.linalg import splu
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'numerical'))
from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian
from work_blocks import condensed_factor,condensed_solve
rows=[];rng=np.random.default_rng(43719)
for name in ['shared144-one','shared168-one','shared168-long-two']:
 r=json.loads((R/'baseline-data'/(name+'.json')).read_text());z=np.load(R/'baseline-data'/(name+'.npz'));m=Mesh(Parameters(**r['parameters']))
 for target in [.0,.1,.125,.2]:
  i=int(np.argmin(abs(z['time']-target)));y=z['state'][i];J=coupled_jacobian(m,m.rhs,z['time'][i],y)
  keep=np.arange(m.size) if m.p.feedback else np.r_[np.arange(m.size-7),np.arange(m.size-4,m.size)]
  J=J[keep][:,keep]
  for c in [1e-6,1e-3,1e-3+2e-3j]:
   A=(eye(len(keep),format='csc')-c*J).tocsc();b=rng.normal(size=len(keep)).astype(A.dtype)
   if np.iscomplexobj(A):b+=1j*rng.normal(size=len(keep))
   t=time.perf_counter();old=splu(A);a=old.solve(b);oldsec=time.perf_counter()-t
   t=time.perf_counter();f=condensed_factor(A,m.position_size,m.internal_size);v=condensed_solve(f,b);newsec=time.perf_counter()-t
   err=float(np.max(abs(a-v))/max(1,np.max(abs(a))));res=float(np.max(abs(A@v-b))/max(1,np.max(abs(A.data))*np.max(abs(v))))
   assert err<1e-8 and res<1e-12,(name,target,c,err,res)
   rows.append(dict(record=name,time=float(z['time'][i]),coefficient=str(c),scaled_difference=err,scaled_residual=res,full_seconds=oldsec,condensed_seconds=newsec))
result=dict(passed=True,checks=len(rows),maximum_difference=max(r['scaled_difference'] for r in rows),maximum_residual=max(r['scaled_residual'] for r in rows),rows=rows,source_sha256=hashlib.sha256((R/'numerical/work_blocks.py').read_bytes()).hexdigest(),check_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out=R/'data/factor-check.json'
if '--write' in sys.argv or not out.exists():out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
