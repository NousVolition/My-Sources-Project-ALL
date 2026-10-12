"""Separate radial and angular under-resolution without changing material laws."""
import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from experiment import run,Parameters
cases=[(24,5),(24,7),(48,3)]
if __name__=='__main__':
    (ROOT/'data/refinement-probe-protocol.json').write_text(json.dumps(dict(meshes=cases,end=.2,alpha=16,rho=160,viscosity=.15,feedback=False,scope='Axis-isolation diagnostic; not tuned material constants'),indent=2)+'\n')
    for n,l in cases:
        run(f'axis{n}x{l}',Parameters(angles=n,layers=l,alpha=16,rho=160,viscosity=.15,preload=.35,feedback=False),end=.2,max_seconds=900)
