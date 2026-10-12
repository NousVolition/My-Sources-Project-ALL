"""Check the Jacobian at the exact difficult continuation starting states."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,hashlib,time,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian

def main():
    rng=np.random.default_rng(823);checks=[]
    for name in ['hard48-one','hard72-one']:
        r=json.loads((ROOT/'parents'/(name+'.json')).read_text());z=np.load(ROOT/'parents'/(name+'.npz'));y=z['state'][-1];m=Mesh(Parameters(**r['parameters']));f=m.rhs(0,y);J=coupled_jacobian(m,m.rhs,0,y)
        for trial in range(3):
            direction=rng.normal(size=m.size);direction[-4:]=0;direction/=np.linalg.norm(direction)
            errors=[]
            for h in [1e-5,1e-6,1e-7,1e-8]:
                fd=(m.rhs(0,y+h*direction)-m.rhs(0,y-h*direction))/(2*h)
                errors.append(float(np.max(abs(fd-J@direction))/max(1,np.max(abs(fd)))))
            checks.append(dict(name=name,trial=trial,steps=[1e-5,1e-6,1e-7,1e-8],errors=errors,passed=min(errors)<5e-5,maximum_speed=float(np.max(np.linalg.norm(f[:m.position_size].reshape(-1,2),axis=1)))))
    result=dict(all_passed=all(c['passed'] for c in checks),checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'data/stiff-jacobian-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['all_passed']
if __name__=='__main__':main()
