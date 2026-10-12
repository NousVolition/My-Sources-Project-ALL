"""Replay actual attempt settings; optionally test fresh starts instead of resumes."""
from pathlib import Path
import sys,json,argparse
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from experiment import run,Parameters
import numpy as np
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--only',nargs='+');ap.add_argument('--from-zero',action='store_true');ap.add_argument('--part',type=int,default=0);ap.add_argument('--parts',type=int,default=1);a=ap.parse_args()
    if not 0<=a.part<a.parts:raise ValueError('Require 0 <= part < parts')
    cases=json.loads((ROOT/'data/executed-protocol.json').read_text())['attempts']
    if a.only:
        if not set(a.only)<=set(c['name'] for c in cases):raise ValueError('Unknown attempt')
        cases=[c for c in cases if c['name'] in a.only]
    for i,c in enumerate(cases):
        if i%a.parts!=a.part:continue
        state=None;parent=None;start=0.;name=c['name']
        if c['parent'] and not a.from_zero:
            parent=c['parent'];raw=np.load((ROOT/parent['metadata']).with_suffix('.npz'));state=raw['state'][-1];start=float(raw['time'][-1])
        elif a.from_zero:name+='-from-zero'
        s=c['audit'];run(name,Parameters(**c['parameters']),end=c['requested_end'],rtol=s['rtol'],atol=s['atol'],method=s['method'],max_step=s['max_step'],max_seconds=s['max_seconds'],initial=c['initial'],state=state,start_time=start,parent=parent)
if __name__=='__main__':main()
