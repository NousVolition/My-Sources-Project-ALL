"""Extract an exact saved prefix for a local independent integration check."""
import argparse,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('name');ap.add_argument('--time',type=float,required=True);a=ap.parse_args();p=(R/a.source).resolve();assert p.is_relative_to(R);r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));ids=np.flatnonzero(abs(z['time']-a.time)<1e-10);assert len(ids)==1;i=int(ids[0]);assert i>0
out=R/'parents'/(a.name+'.json');assert not out.exists() and not out.with_suffix('.npz').exists()
r.update(name=a.name,completed_end=float(z['time'][i]),diagnostics=r['diagnostics'][:i+1],failure=dict(time=float(z['time'][i]),message='Exact derived saved-state prefix; original source continued'))
r['checkpoint_extract']=dict(source=p.relative_to(R).as_posix(),source_metadata_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_raw_sha256=hashlib.sha256(p.with_suffix('.npz').read_bytes()).hexdigest(),index=i,audit_scope='Original full integration, not prefix step counts')
out.write_text(json.dumps(r,indent=2)+'\n');np.savez_compressed(out.with_suffix('.npz'),**{k:z[k][:i+1] if k in ['state','time'] else z[k] for k in z.files});print(out)
