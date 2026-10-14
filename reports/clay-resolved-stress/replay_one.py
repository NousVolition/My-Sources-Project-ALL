"""Replay one immutable protocol; NumPy imports before the SciPy thread setup."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
name=sys.argv[1];protocol=json.loads((ROOT/'data/executed-protocol.json').read_text())['attempts'];r=next(v for v in protocol if v['name']==name)
folders=[ROOT/'solver',*sorted((ROOT/'solver-versions').glob('*'))];matches=[f for f in folders if all((f/n).exists() and hashlib.sha256((f/n).read_bytes()).hexdigest()==h for n,h in r['source_hashes'].items())];assert len(matches)==1
folder=matches[0];sys.path.insert(0,str(folder))
import numpy as np
from experiment import run,Parameters
kw=dict(end=r['requested_end'],rtol=r['audit']['rtol'],atol=r['audit']['atol'],method=r['audit']['method'],max_step=r['audit']['max_step'],max_seconds=r['audit']['max_seconds'],initial=r['initial'],record_step=r['record_step'])
if r['parent']:
 z=np.load((ROOT/r['parent']['metadata']).with_suffix('.npz'));kw.update(state=z['state'][-1],start_time=r['start_time'],parent=r['parent'])
run(r['name'],Parameters(**r['parameters']),**kw)
