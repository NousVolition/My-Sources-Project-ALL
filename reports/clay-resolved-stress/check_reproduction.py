"""Check a fresh replay against saved raw states, sources and physical setup."""
from pathlib import Path
import argparse,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser();ap.add_argument('copy',type=Path);ap.add_argument('--additional-copy',type=Path,nargs='*');ap.add_argument('--only',nargs='+',required=True);a=ap.parse_args();checks=[]
 for name in a.only:
  paths=[*(ROOT/'solver/data').glob(name+'.json'),*(ROOT/'solver-versions').glob('*/data/'+name+'.json')];assert len(paths)==1;p=paths[0]
  candidates=[a.copy]+(a.additional_copy or []);copies=[(folder,folder/p.relative_to(ROOT)) for folder in candidates if (folder/p.relative_to(ROOT)).exists()];assert len(copies)==1
  folder,q=copies[0];r=json.loads(p.read_text());s=json.loads(q.read_text());u=np.load(p.with_suffix('.npz'));v=np.load(q.with_suffix('.npz'))
  assert r['source_hashes']==s['source_hashes'] and r['parameters']==s['parameters'] and r['initial']==s['initial'] and not r['failure'] and not s['failure']
  equal=all(np.array_equal(u[k],v[k]) for k in u.files);assert equal
  checks.append(dict(name=name,replay_copy_name=folder.name,all_raw_arrays_bitwise_equal=equal,passed=True))
 result=dict(all_passed=True,checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(ROOT/'data/reproduction-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
