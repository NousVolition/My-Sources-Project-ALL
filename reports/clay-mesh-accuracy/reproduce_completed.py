"""Replay the released-model protocol in a fresh copy prepared by prepare_rerun."""
import sys,json,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from experiment import run,Parameters
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--part',type=int,default=0);ap.add_argument('--parts',type=int,default=1);ap.add_argument('--only',nargs='*');a=ap.parse_args()
    if not 0<=a.part<a.parts:raise ValueError('Require 0 <= part < parts')
    cases=json.loads((ROOT/'data/executed-protocol.json').read_text())['attempts']
    if a.only:
        assert set(a.only)<=set(c['name'] for c in cases);cases=[c for c in cases if c['name'] in a.only]
    for i,c in enumerate(cases):
        if i%a.parts!=a.part:continue
        s=c['audit'];assert s['max_rhs']==2000000 and s['max_steps']==50000
        run(c['name'],Parameters(**c['parameters']),end=c['requested_end'],initial=c['initial'],method=s['method'],rtol=s['rtol'],atol=s['atol'],max_step=s['max_step'],max_seconds=s['max_seconds'])
if __name__=='__main__':main()
