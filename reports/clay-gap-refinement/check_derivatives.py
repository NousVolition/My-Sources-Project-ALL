"""Independent energy differences and RHS derivatives at actual stress states."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,tempfile,subprocess,hashlib,numpy as np
R=Path(__file__).resolve().parent
if len(sys.argv)>1 and sys.argv[1]=='--zero-worker':
 source=sys.argv[2];sys.path.insert(0,str(R/source));from mesh_hug import Mesh,Parameters
 from jacobian import coupled_jacobian
 r=json.loads((R/'baseline-data/shared168-one.json').read_text());z=np.load(R/'baseline-data/shared168-one.npz');i=int(np.argmin(abs(z['time']-.125)));p=dict(r['parameters']);
 if source=='solver':p['regularization_length']=0
 m=Mesh(Parameters(**p));y=z['state'][i];x,rest,l=m.unpack(y);e,g,*_=m.energy_force(x,rest);d=np.random.default_rng(317).normal(size=m.size)
 np.savez(sys.argv[3],energy=e,gradient=g,rhs=m.rhs(.125,y),jac_direction=coupled_jacobian(m,m.rhs,.125,y)@d);raise SystemExit
sys.path.insert(0,str(R/'solver'));from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian
rows=[];rng=np.random.default_rng(713551)
for name in ['length100-72-one','length100-96-one','harder72-one','weak72-two','length100-168-two','length100-144-one','length100-168-one','harder120-one']:
 p=R/'solver/data'/(name+'.json')
 if not p.exists():continue
 r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));m=Mesh(Parameters(**r['parameters']))
 for target in [.125,.2,float(z['time'][-1])]:
  i=int(np.argmin(abs(z['time']-target)));t=float(z['time'][i]);y=z['state'][i];x,rest,l=m.unpack(y);J=coupled_jacobian(m,m.rhs,t,y);g=m.energy_force(x,rest)[1]
  for j in range(3):
   d=rng.normal(size=m.size);d/=np.linalg.norm(d);eps=1e-7;fd=(m.rhs(t,y+eps*d)-m.rhs(t,y-eps*d))/(2*eps);exact=J@d;err=float(np.max(abs(fd-exact))/max(1,np.max(abs(fd))))
   dx=d[:m.position_size].reshape(-1,2);ef=(m.energy_force(x+eps*dx,rest)[0]-m.energy_force(x-eps*dx,rest)[0])/(2*eps);eg=float(np.sum(g*dx));energy_err=abs(ef-eg)/max(1,abs(eg))
   rows.append(dict(record=name,time=t,direction=j,rhs_error=err,energy_gradient_error=energy_err));assert err<2e-5 and energy_err<2e-5,rows[-1]
with tempfile.TemporaryDirectory(dir=R/'data',prefix='derivative-') as folder:
 assert Path(folder).resolve().is_relative_to((R/'data').resolve())
 files=[Path(folder)/f'{s}.npz' for s in ['numerical','solver']]
 for source,out in zip(['numerical','solver'],files):subprocess.run([sys.executable,str(Path(__file__)), '--zero-worker',source,str(out)],check=True)
 a,b=[np.load(p) for p in files];zero={k:float(np.max(abs(a[k]-b[k]))/max(1,np.max(abs(a[k])))) for k in a.files};a.close();b.close();assert max(zero.values())<1e-12,zero
result=dict(passed=True,probes=len(rows),rows=rows,zero_length_unchanged_model=zero,source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'solver').glob('*.py')},check_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out=R/'data/derivative-check.json'
if '--write' in sys.argv or not out.exists():out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(passed=True,probes=len(rows),zero_length_max=max(zero.values()),maximum_rhs_error=max(r['rhs_error'] for r in rows),maximum_energy_gradient_error=max(r['energy_gradient_error'] for r in rows))))
