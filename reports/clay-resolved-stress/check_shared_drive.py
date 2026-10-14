"""Independent Lorenz-reference and passive-work linear-solve controls."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver-versions/shared'))
from reference_drive import reference
from work_blocks import factor,solve
from scipy.sparse import bmat,csc_matrix,eye
from scipy.sparse.linalg import splu
from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian
def main():
 a=reference(.5,[-8,8,27],500);b=reference(.5,[-8,8,27],500,'Radau');t=np.linspace(0,.5,1001);u=a.sol(t);v=b.sol(t);error=float(np.max(abs(u-v)/np.maximum(1,np.max(abs(v),axis=1)[:,None])));assert error<1e-8
 rng=np.random.default_rng(72013);n=30;A=rng.normal(size=(n,n));A=A.T@A+np.eye(n);C=rng.normal(size=(4,n))*1e6;M=bmat([[csc_matrix(A),None],[csc_matrix(C),eye(4)]],format='csc');known=rng.normal(size=n+4);rhs=M@known;x=solve(factor(M),rhs);err=float(max(abs(x-known))/max(1,max(abs(known))));assert err<1e-7
 p=ROOT/'solver/data/fine96-one.json';r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));m=Mesh(Parameters(**r['parameters']));y=z['state'][int(np.argmin([d['min_area_ratio'] for d in r['diagnostics']]))];J=coupled_jacobian(m,m.rhs,0,y);keep=np.r_[np.arange(m.size-7),np.arange(m.size-4,m.size)];J=J[keep][:,keep];real=[]
 for step in [1e-6,1e-4,1e-3]:
  A=eye(len(keep),format='csc')-step*J;bvec=rng.normal(size=len(keep));new=solve(factor(A),bvec);old=splu(A).solve(bvec);difference=float(max(abs(new-old))/max(1,max(abs(old))));residual=float(max(abs(A@new-bvec))/max(1,max(abs(bvec))));assert difference<1e-7 and residual<1e-7
  real.append(dict(step=step,solution_difference=difference,residual=residual))
 result=dict(all_passed=True,lorenz_reference_scaled_error=error,lorenz_samples=len(t),passive_work_solution_error=err,real_stiff_matrix_checks=real,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(ROOT/'data/shared-drive-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
