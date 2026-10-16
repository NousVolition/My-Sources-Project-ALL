"""Independent 32/64/128 point contact and observable quadrature checks."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'solver'))
from mesh_hug import Mesh,Parameters
from clay_hug import object_gap
from contact_hug import barrier
def integrate(m,x,q):
 t,w=np.polynomial.legendre.leggauss(q);t=(t+1)/2;w=w/2;i,j=m.inner_i,m.inner_j;pts=(1-t[None,:,None])*x[i,None]+t[None,:,None]*x[j,None];g,n=object_gap(pts.reshape(-1,2),m.p);e,d=barrier(g,m.p.barrier_width,m.p.contact);weights=(m.surface_length[:,None]*w).ravel();cg=((weights*d)[:,None]*n).reshape(len(i),len(t),2)
 grad=m.scatter(i,((1-t[None,:,None])*cg).sum(axis=1))+m.scatter(j,(t[None,:,None]*cg).sum(axis=1));length=np.linalg.norm(x[j]-x[i],axis=1);rms=np.sqrt(np.sum(length[:,None]*w*g.reshape(len(i),q)**2)/sum(length))
 return float(weights@e),grad,float(rms)
rows=[]
for name in ['length100-72-one','length100-96-one','length100-144-one','length100-168-one','harder72-one','harder96-one','harder120-one','middle72-two','weak72-two-tighter']:
 p=R/'solver/data'/(name+'.json')
 if not p.exists():continue
 r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));m=Mesh(Parameters(**r['parameters']))
 for target in [.125,.2,float(z['time'][-1])]:
  i=int(np.argmin(abs(z['time']-target)));x=m.unpack(z['state'][i])[0];ref=integrate(m,x,128)
  for q in [32,64]:
   v=integrate(m,x,q);e=abs(v[0]-ref[0])/max(1,abs(ref[0]));g=float(np.max(abs(v[1]-ref[1]))/max(1,np.max(abs(ref[1]))));gap=abs(v[2]-ref[2])/max(.01,abs(ref[2]));rows.append(dict(record=name,time=float(z['time'][i]),quadrature=q,energy_scaled=e,force_scaled=g,gap_scaled=gap));assert e<1e-5 and g<1e-4 and gap<1e-7,rows[-1]
out=dict(passed=True,checks=len(rows),reference_points=128,rows=rows,source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'solver').glob('*.py')},check_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());path=R/'data/quadrature-check.json'
if '--write' in sys.argv or not path.exists():path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(passed=True,checks=len(rows),max_force=max(r['force_scaled'] for r in rows),max_gap=max(r['gap_scaled'] for r in rows))))
