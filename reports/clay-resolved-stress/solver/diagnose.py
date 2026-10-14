"""Controlled spatial checks; no trajectory fitting or error-target adjustment."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,numpy as np
from dataclasses import replace
from mesh_hug import Mesh,Parameters
from cell_hug import Mesh as Previous
from clay_hug import Mesh as SpringMesh
from contact_hug import barrier
def main():
    out=[]
    F=np.array([[1.8,.2],[0,.65]])
    for n,l in [(24,3),(48,5),(72,7),(96,9)]:
        p=Parameters(angles=n,layers=l);m=Mesh(p);old=Previous(p)
        x=m.initial@F.T;oldbase=SpringMesh.material(m,x,m.l0)[0]
        Q=2*np.linalg.det(F)/np.sum(F*F)
        energy=(m.material(x,m.l0)[0]-oldbase)/sum(m.a0)
        expected=barrier(np.array([np.linalg.det(F)-p.area_floor]),p.area_activation-p.area_floor,p.area_strength)[0][0]+barrier(np.array([Q-p.shape_floor]),p.shape_activation-p.shape_floor,p.shape_strength)[0][0]
        # Exact same affine deformation, area-normalized. New invariant protection
        # must reproduce its closed formula regardless of triangulation.
        out.append(dict(angles=n,layers=l,new_affine_barrier_density=energy,exact_density=expected,absolute_error=abs(energy-expected),old_affine_barrier_density=(old.material(x,m.l0)[0]-oldbase)/sum(m.a0)))
    result=dict(affine_checks=out,all_passed=all(v['absolute_error']<1e-11 for v in out))
    (Path(__file__).parent/'data/spatial-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['all_passed']
if __name__=='__main__':main()
