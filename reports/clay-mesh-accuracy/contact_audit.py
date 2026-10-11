import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,numpy as np
from pathlib import Path
root=Path(__file__).resolve().parent;sys.path.insert(0,str(root/'model'))
from mesh_hug import Mesh,Parameters
from contact_hug import barrier
rows=[]
for path in [root/'data/revised24',root/'model/data/mesh48-one']:
    a=np.load(path.with_suffix('.npz'));r=json.loads(path.with_suffix('.json').read_text());m=Mesh(Parameters(**r['parameters']));n=m.p.angles
    for t,y in zip(a['time'],a['state']):
        x=m.unpack(y)[0];g=m.self_gap(x,False);active=g<m.p.self_width
        if not active.any():continue
        v=m.si[active];s=m.sa[active];b=m.sb[active]
        same=(v//n==s//n)
        distance=lambda a,b:np.minimum(abs(a%n-b%n),n-abs(a%n-b%n))
        refarc=np.where(same,np.minimum(distance(v,s),distance(v,b))*2*np.pi/n*np.linalg.norm(m.initial[v],axis=1),np.inf)
        e=barrier(g[active],m.p.self_width,m.p.self_strength)[0]*m.sw[active]
        rows.append(dict(run=r['name'],time=float(t),active_pairs=int(active.sum()),same_boundary_pairs=int(same.sum()),nearby_reference_pairs=int((refarc<.25).sum()),total_self_energy=float(e.sum()),nearby_reference_energy=float(e[refarc<.25].sum()),minimum_gap=float(g.min())))
result=dict(frames=rows)
(root/'data/contact-neighborhood-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(rows,indent=2))
