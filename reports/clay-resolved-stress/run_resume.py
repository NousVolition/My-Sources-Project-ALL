"""Continue an accepted common-input state without resetting work or memory."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,argparse,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver-versions/shared-resume'))
from experiment import run,Parameters

def main():
 ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--parent',required=True);ap.add_argument('--end',type=float,default=.5);ap.add_argument('--seconds',type=float,default=3600);ap.add_argument('--rtol',type=float,default=2e-7);ap.add_argument('--next-two-records',action='store_true');a=ap.parse_args()
 path=ROOT/'parents'/(a.parent+'.json');r=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));assert not r['parameters']['feedback'] and r['audit']['split_one_way'] and z['time'][-1]==r['completed_end']
 start=float(z['time'][-1]);end=min(.5,(np.floor(start/.005+1e-9)+2)*.005) if a.next_two_records else a.end
 originals=[*(ROOT/'solver/data').glob(a.parent+'.json'),*(ROOT/'solver-versions').glob('*/data/'+a.parent+'.json')]
 source=r.get('checkpoint_extract',{}).get('source_record')
 if not source:assert len(originals)==1;source=originals[0].relative_to(ROOT).as_posix()
 parent=dict(metadata=path.relative_to(ROOT).as_posix(),metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),raw_sha256=hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest(),checkpoint_origin='local-recording',source_record=source)
 run(a.name,Parameters(**r['parameters']),state=z['state'][-1],start_time=start,initial=r['initial'],parent=parent,end=float(end),max_seconds=a.seconds,rtol=a.rtol,atol=a.rtol/100)
if __name__=='__main__':main()
