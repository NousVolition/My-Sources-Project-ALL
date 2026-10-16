"""Create a result-free copy and replay one recorded experiment exactly."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import argparse,sys,json,shutil,subprocess,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('record');ap.add_argument('destination',type=Path);a=ap.parse_args();src=(R/a.record).resolve();r=json.loads(src.read_text());dest=a.destination.resolve();dest.mkdir(exist_ok=False)
for folder in ['baseline','baseline-fresh','numerical','solver']:
 (dest/folder).mkdir()
 for p in (R/folder).glob('*.py'):shutil.copy2(p,dest/folder/p.name)
for folder in ['parents','baseline-data']:shutil.copytree(R/folder,dest/folder)
for name in ['replay_worker.py','audit_one.py']:shutil.copy2(R/name,dest/name)
payload=dest/'replay-request.json';payload.write_text(json.dumps(dict(record=r,original_path=a.record)))
subprocess.run([sys.executable,str(dest/'replay_worker.py')],check=True,cwd=dest)
out=dest/a.record;z=np.load(src.with_suffix('.npz'));v=np.load(out.with_suffix('.npz'));equal={k:bool(np.array_equal(z[k],v[k])) for k in z.files};assert all(equal.values()),equal
proof=dict(record=a.record,destination=str(dest),all_arrays_bitwise_equal=True,arrays=equal,source_metadata_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),source_raw_sha256=hashlib.sha256(src.with_suffix('.npz').read_bytes()).hexdigest())
(R/'data'/('replay-'+r['name']+'.json')).write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
