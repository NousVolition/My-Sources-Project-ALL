"""Independent contact quadrature-order checks at actually computed states."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,numpy as np
from dataclasses import replace
from analyze import load,Mesh,Parameters
ROOT=Path(__file__).resolve().parent
def main():
 checks=[]
 for name in ['extend72-one','fine96-one','fine120-one','indexed120-one','fine144-one','extend144-one','shared144-one','shared168-one','extend168-one']:
  paths=[*(ROOT/'solver/data').glob(name+'.json'),*(ROOT/'solver-versions').glob('*/data/'+name+'.json')]
  if not paths:continue
  assert len(paths)==1;path=paths[0]
  r,z,t,y,d=load(path);p=Parameters(**r['parameters']);meshes={q:Mesh(replace(p,quadrature=q)) for q in [16,32,64,128]}
  ids=sorted(set([int(np.argmin([v['certified_clearance'] for v in d])),int(np.argmin([v['min_area_ratio'] for v in d])),int(np.argmin([v['gap_arclength_rms'] for v in d])),len(t)-1]))
  for i in ids:
   x,rest,l=meshes[32].unpack(y[i]);values={q:m.energy_force(x,rest)[:2] for q,m in meshes.items()};ref_e,ref_g=values[128];errors=[]
   for q in [16,32,64]:
    e,g=values[q];errors.append(dict(order=q,energy_scaled=float(abs(e-ref_e)/max(1,abs(ref_e))),force_scaled=float(np.max(abs(g-ref_g))/max(1,np.max(abs(ref_g))))))
   passed=all(max(v['energy_scaled'],v['force_scaled'])<1e-3 for v in errors if v['order'] in [32,64]);checks.append(dict(record=name,time=float(t[i]),errors=errors,passed=passed))
 result=dict(all_passed=all(v['passed'] for v in checks),checks=checks,target=.001,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(ROOT/'data/quadrature-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
