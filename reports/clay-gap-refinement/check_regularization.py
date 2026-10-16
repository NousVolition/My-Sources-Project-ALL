"""Affine patch, rigid-motion, gradient and full RHS Jacobian checks."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'solver'))
from mesh_hug import Mesh,Parameters
from regularization import energy_gradient
from jacobian import coupled_jacobian
from reference_drive import reference
from scipy.integrate import solve_ivp
rows=[];rng=np.random.default_rng(88134)
for n in [24,48,72,144,168,192]:
 m=Mesh(Parameters(angles=n,layers=n//12+1));x=m.initial;F=np.array([[1.1,.15],[-.05,.9]])
 affine=x@F.T+np.array([.3,-.7]);e,g=energy_gradient(m,affine)
 assert abs(e)<1e-9 and np.max(abs(g))<1e-7,(n,e,np.max(abs(g)))
 d=rng.normal(size=x.shape);v=rng.normal(size=x.shape);d*=.001/np.max(abs(d));v/=np.linalg.norm(v);y=x+d;eps=1e-6
 e,g=energy_gradient(m,y);fd=(energy_gradient(m,y+eps*v)[0]-energy_gradient(m,y-eps*v)[0])/(2*eps);exact=float(np.sum(g*v));err=abs(fd-exact)/max(1,abs(exact));assert err<1e-7 and e>=-1e-12
 rows.append(dict(n=n,affine_energy=float(energy_gradient(m,affine)[0]),gradient_error=err))
for two in [False,True]:
 m=Mesh(Parameters(angles=48,layers=5,feedback=two,alpha=128,rho=500,viscosity=.01875));y=m.pack((-8,8,27));y[:m.position_size]+=.0001*rng.normal(size=m.position_size);y[m.position_size:m.position_size+m.internal_size]*=1.001
 J=coupled_jacobian(m,m.rhs,0,y)
 for i in range(8):
  d=rng.normal(size=m.size);d/=np.linalg.norm(d);eps=2e-7;fd=(m.rhs(0,y+eps*d)-m.rhs(0,y-eps*d))/(2*eps);exact=J@d;err=float(np.max(abs(fd-exact))/max(1,np.max(abs(fd))));assert err<2e-5,err
  rows.append(dict(feedback=two,directional=i,error=err))
a=reference(.5,(-8,8,27),500);b=reference(1.5,(-8,8,27),500);ts=np.linspace(0,.5,501);assert np.array_equal(a.sol(ts),b.sol(ts))
c=reference(1.5,(-8,8,27),500,method='Radau');err=float(np.max(abs(c.sol(np.linspace(0,1.5,1501))-b.sol(np.linspace(0,1.5,1501)))/np.maximum(1,abs(c.sol(np.linspace(0,1.5,1501))))));assert err<1e-8
result=dict(passed=True,rows=rows,reference_extension_bitwise_preserved=True,independent_reference_scaled=err,source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'solver').glob('*.py')},check_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out=R/'data/regularization-check.json'
if '--write' in sys.argv or not out.exists():out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(passed=True,checks=len(rows),reference_scaled=err)))
