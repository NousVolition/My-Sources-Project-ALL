import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import importlib.util,json
from pathlib import Path
import numpy as np
from dataclasses import replace
from mesh_hug import Mesh,Parameters
from jacobian import coupled_jacobian,spatial_groups
root=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('broad_jacobian',root.parent/'jacobian.py');broad=importlib.util.module_from_spec(spec);spec.loader.exec_module(broad)
rows=[]
for label,p,scale in [('initial',Parameters(),1),('active-object',Parameters(),.86),('active-self',Parameters(self_width=.5),1),('fine-mesh',Parameters(angles=48,layers=5),1)]:
    m=Mesh(p);y=m.pack((2,3,4));y[:m.position_size]*=scale
    a=coupled_jacobian(m,m.rhs,0,y).toarray();b=broad.coupled_jacobian(m,m.rhs,0,y).toarray()
    e=float(np.max(abs(a-b))/max(1,np.max(abs(b))))
    rows.append(dict(name=label,error=e,passed=e<1e-12,probe_groups=len(spatial_groups(m,m.unpack(y)[0])[1])))
result=dict(checks=rows,all_passed=all(r['passed'] for r in rows))
(root/'data/jacobian-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['all_passed']
