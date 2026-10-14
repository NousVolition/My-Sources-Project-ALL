"""Fresh high-resolution integrations with conservative spatial self-contact."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,argparse,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver-versions/indexed'))
from experiment import run,Parameters
def main():
 ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--n',type=int,default=144);ap.add_argument('--layers',type=int,default=13);ap.add_argument('--two',action='store_true');ap.add_argument('--end',type=float,default=.5);ap.add_argument('--seconds',type=float,default=10800);ap.add_argument('--parent');a=ap.parse_args();kw=dict(end=a.end,max_seconds=a.seconds)
 if a.parent:
  path=ROOT/'parents'/(a.parent+'.json');r=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));assert z['time'][-1]==r['completed_end'] and np.isfinite(z['state']).all()
  p=Parameters(**r['parameters']);kw.update(state=z['state'][-1],start_time=float(z['time'][-1]),initial=r['initial'],parent=dict(metadata=path.relative_to(ROOT).as_posix(),metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),raw_sha256=hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest(),checkpoint_origin='local-recording',source_record='solver-versions/indexed/data/'+a.parent+'.json'))
 else:p=Parameters(angles=a.n,layers=a.layers,feedback=a.two,alpha=128,rho=500,viscosity=.01875,preload=.35)
 run(a.name,p,**kw)
if __name__=='__main__':main()
