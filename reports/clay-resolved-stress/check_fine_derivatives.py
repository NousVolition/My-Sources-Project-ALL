"""Force and Jacobian probes at actually recorded fine-mesh high-load states."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian
def main():
 checks=[];rng=np.random.default_rng(1731)
 for name in ['fine96-one','fine120-one','fine120-two','indexed120-one','fine144-one','extend144-one','fine144-two','shared144-one','shared168-one','extend168-one','shared168-two','shared168-long-two']:
  paths=[*(ROOT/'solver/data').glob(name+'.json'),*(ROOT/'solver-versions').glob('*/data/'+name+'.json')]
  if not paths:continue
  assert len(paths)==1;p=paths[0]
  r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));m=Mesh(Parameters(**r['parameters']));ids=sorted(set([int(np.argmin([d['min_area_ratio'] for d in r['diagnostics']])),int(np.argmin([d['self_clearance'] for d in r['diagnostics']])),len(z['time'])-1]))
  for i in ids:
   y=z['state'][i];x,rest,l=m.unpack(y);e,g,*_=m.energy_force(x,rest);J=coupled_jacobian(m,m.rhs,z['time'][i],y)
   for j in range(2):
    direction=rng.normal(size=m.size);direction/=np.linalg.norm(direction);field=direction[:m.position_size].reshape(m.n,2);errors=[]
    for eps in [3e-7,1e-7,3e-8]:
     try:
      fd=(m.rhs(0,y+eps*direction)-m.rhs(0,y-eps*direction))/(2*eps);derivative=float(np.max(abs(fd-J@direction))/max(1,np.max(abs(fd))))
      energy=(m.energy_force(x+eps*field,rest)[0]-m.energy_force(x-eps*field,rest)[0])/(2*eps);gradient=float(abs(energy-np.sum(g*field))/max(1,abs(energy)));errors.append(dict(epsilon=eps,jacobian=derivative,force=gradient))
     except FloatingPointError:continue
    assert errors
    best_j=min(v['jacobian'] for v in errors);best_f=min(v['force'] for v in errors);passed=best_j<2e-4 and best_f<2e-4;assert passed,(name,i,best_j,best_f)
    checks.append(dict(record=name,index=i,time=float(z['time'][i]),direction=j,probes=errors,passed=passed))
 assert checks
 result=dict(all_passed=True,checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(ROOT/'data/fine-derivative-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,checks=len(checks))))
if __name__=='__main__':main()
