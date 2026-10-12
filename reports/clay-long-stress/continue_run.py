"""Continue a verified recorded state with its clock and work history intact."""
from pathlib import Path
import json,hashlib,sys,argparse
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from experiment import run,Parameters
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('parent');ap.add_argument('name');ap.add_argument('--end',type=float,required=True);ap.add_argument('--seconds',type=float,default=2400);ap.add_argument('--method',default='BDF');ap.add_argument('--rtol',type=float,default=2e-9);a=ap.parse_args()
    path=(ROOT/'parents'/(a.parent+'.json'));r=json.loads(path.read_text());raw=path.with_suffix('.npz');z=np.load(raw)
    assert z['time'][-1]==r['completed_end'] and np.isfinite(z['state']).all()
    provenance=dict(metadata=path.relative_to(ROOT).as_posix(),metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),upstream_commit='fb5e18e33fb48cd1928c556c45903735bbc86b3b')
    run(a.name,Parameters(**r['parameters']),end=a.end,max_seconds=a.seconds,method=a.method,rtol=a.rtol,atol=a.rtol/100,initial=r['initial'],state=z['state'][-1],start_time=float(z['time'][-1]),parent=provenance)
if __name__=='__main__':main()
