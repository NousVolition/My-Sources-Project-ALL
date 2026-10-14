"""Run a declared fresh case or an exact, hash-recorded checkpoint continuation."""
from pathlib import Path
import sys,json,hashlib,argparse
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
import numpy as np
from experiment import run,Parameters,hashes
def main():
 ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--n',type=int,default=96);ap.add_argument('--layers',type=int,default=9);ap.add_argument('--two',action='store_true');ap.add_argument('--end',type=float,default=.5);ap.add_argument('--seconds',type=float,default=7200);ap.add_argument('--method',default='BDF');ap.add_argument('--rtol',type=float,default=2e-7);ap.add_argument('--record-step',type=float,default=.005);ap.add_argument('--max-step',type=float,default=.005);ap.add_argument('--parent');ap.add_argument('--alpha',type=float,default=128);ap.add_argument('--rho',type=float,default=500);ap.add_argument('--viscosity',type=float,default=.01875);a=ap.parse_args()
 kw=dict(end=a.end,max_seconds=a.seconds,method=a.method,rtol=a.rtol,atol=a.rtol/100,record_step=a.record_step,max_step=a.max_step)
 if a.parent:
  path=ROOT/'parents'/(a.parent+'.json');r=json.loads(path.read_text());raw=path.with_suffix('.npz');z=np.load(raw)
  assert z['time'][-1]==r['completed_end'] and np.isfinite(z['state']).all()
  local=r['source_hashes']==hashes()
  origin=dict(checkpoint_origin='local-recording',source_record='solver/data/'+a.parent+'.json') if local else dict(checkpoint_origin='published-recording',upstream_commit='fbf4ab15bc1e0875929955e5f1720748464768b6')
  kw.update(initial=r['initial'],state=z['state'][-1],start_time=float(z['time'][-1]),parent=dict(metadata=path.relative_to(ROOT).as_posix(),metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),**origin))
  p=Parameters(**r['parameters'])
 else:p=Parameters(angles=a.n,layers=a.layers,feedback=a.two,alpha=a.alpha,rho=a.rho,viscosity=a.viscosity,preload=.35)
 run(a.name,p,**kw)
if __name__=='__main__':main()
