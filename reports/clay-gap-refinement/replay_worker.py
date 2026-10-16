import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent;r=json.loads((R/'replay-request.json').read_text())['record']
folders=[R/k for k in ['baseline','baseline-fresh','numerical','solver']];choices=[f for f in folders if all((f/n).exists() and hashlib.sha256((f/n).read_bytes()).hexdigest()==h for n,h in r['source_hashes'].items())];assert len(choices)==1
sys.path.insert(0,str(choices[0]));from experiment import run,Parameters
kw={}
if r.get('parent'):
 paths=[p for p in R.rglob('*.json') if hashlib.sha256(p.read_bytes()).hexdigest()==r['parent']['metadata_sha256']];assert paths;p=paths[0];z=np.load(p.with_suffix('.npz'));kw=dict(state=z['state'][-1],parent=r['parent'])
a=r['audit'];run(r['name'],Parameters(**r['parameters']),end=r['requested_end'],rtol=a['rtol'],atol=a['atol'],method=a['method'],max_step=a['max_step'],max_seconds=a['max_seconds'],initial=r['initial'],start_time=r['start_time'],record_step=r['record_step'],**kw)
