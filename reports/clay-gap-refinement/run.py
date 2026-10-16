"""Run the declared material experiment; never overwrite a recording."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import argparse,sys,json,hashlib
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--source',choices=['baseline','numerical','solver'],default='solver')
ap.add_argument('--parent')
ap.add_argument('--n',type=int,default=168)
ap.add_argument('--layers',type=int,default=15)
ap.add_argument('--end',type=float,default=.5)
ap.add_argument('--seconds',type=float,default=7200)
ap.add_argument('--rtol',type=float,default=2e-7)
ap.add_argument('--method',choices=['BDF','Radau'],default='BDF')
ap.add_argument('--two',action='store_true')
ap.add_argument('--alpha',type=float,default=128)
ap.add_argument('--rho',type=float,default=500)
ap.add_argument('--viscosity',type=float,default=.01875)
ap.add_argument('--length',type=float,default=.1)
ap.add_argument('--capacity',type=float,default=.01)
ap.add_argument('--perturb',type=float,default=0.)
a=ap.parse_args()
sys.path.insert(0,str(R/a.source));from experiment import run,Parameters
kw=dict(initial=(-8+a.perturb,8,27));extra=dict(regularization_length=a.length) if a.source=='solver' else {};p=Parameters(angles=a.n,layers=a.layers,feedback=a.two,alpha=a.alpha,rho=a.rho,viscosity=a.viscosity,preload=.35,lorenz_capacity=a.capacity,**extra)
if a.parent:
 path=(R/a.parent).resolve();assert path.is_relative_to(R);r=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));p=Parameters(**r['parameters']);kw=dict(state=z['state'][-1],start_time=float(z['time'][-1]),initial=r['initial'],parent=dict(metadata=path.relative_to(R).as_posix(),metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),raw_sha256=hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest()))
run(a.name,p,end=a.end,rtol=a.rtol,atol=a.rtol/100,method=a.method,max_seconds=a.seconds,**kw)
