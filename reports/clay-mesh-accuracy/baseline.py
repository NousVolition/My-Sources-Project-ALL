"""Preserve the previous mesh failure and recompute an independent gap metric."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,shutil,hashlib,numpy as np
from mesh_hug import Mesh,Parameters
root=Path(__file__).resolve().parent
def main():
    previous=root.parent/'clay-cell-stress/data';out=root/'historical';out.mkdir(exist_ok=True)
    rows=[]
    for label,name in [('coarse-one','prior-extreme-new-bdf-fine-one'),('fine-one','mesh48-hard-one'),('coarse-two','prior-extreme-new-bdf-fine-two'),('fine-two','mesh48-hard-two')]:
        a=np.load(previous/(name+'.npz'));r=json.loads((previous/(name+'.json')).read_text())
        i=int(np.argmin(abs(a['time']-.2)));assert abs(a['time'][i]-.2)<1e-10
        m=Mesh(Parameters(**r['parameters']));x=m.unpack(a['state'][i])[0]
        # Geometry-only diagnostic: do not evaluate the revised energy on old states.
        from clay_hug import object_gap
        t,w=np.polynomial.legendre.leggauss(64);t=(t+1)/2;w=w/2
        ai=m.inner_i;bi=m.inner_j
        points=(1-t[None,:,None])*x[ai,None]+t[None,:,None]*x[bi,None]
        g=object_gap(points.reshape(-1,2),m.p)[0].reshape(len(ai),-1)
        lengths=np.linalg.norm(x[bi]-x[ai],axis=1)
        arc=float(np.sqrt(np.sum(lengths[:,None]*w*g*g)/lengths.sum()))
        samples=np.linspace(0,1,17);s=(1-samples[None,:,None])*x[ai,None]+samples[None,:,None]*x[bi,None]
        gap=float(np.sqrt(np.mean(np.maximum(object_gap(s.reshape(-1,2),m.p)[0],0)**2)))
        np.savez_compressed(out/(label+'.npz'),time=a['time'][i:i+1],positions=x,triangles=m.tri,outline=a['outline'])
        rows.append(dict(label=label,source_name=name,time=.2,angles=m.p.angles,layers=m.p.layers,gap_legacy=gap,gap_arclength=arc,raw_sha256=hashlib.sha256((previous/(name+'.npz')).read_bytes()).hexdigest(),source_commit='893a778806bff49e400dbab6ce3eea4bd08f4dbc',integration_tolerance=r['audit']['rtol']))
    (out/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
