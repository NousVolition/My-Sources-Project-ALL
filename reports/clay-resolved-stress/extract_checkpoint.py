"""Losslessly extract a saved prefix; retain the completed source and its audit."""
import argparse,json,hashlib,numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('name');ap.add_argument('--time',type=float,required=True);a=ap.parse_args();p=a.source.resolve();r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));ids=np.flatnonzero(abs(z['time']-a.time)<1e-10);assert len(ids)==1;i=int(ids[0]);assert i>0
 out=ROOT/'parents'/(a.name+'.json');assert not out.exists() and not out.with_suffix('.npz').exists();r.update(name=a.name,completed_end=float(z['time'][i]),failure=dict(time=float(z['time'][i]),message='Derived saved-state checkpoint; source integration continued'),diagnostics=r['diagnostics'][:i+1])
 r['checkpoint_extract']=dict(source_record=p.relative_to(ROOT).as_posix(),source_metadata_sha256=sha(p),source_raw_sha256=sha(p.with_suffix('.npz')),index=i,time=float(z['time'][i]),audit_scope='Original full source integration; not reinterpreted as prefix step counts')
 arrays={k:z[k][:i+1] if k in ['state','time'] else z[k] for k in z.files};np.savez_compressed(out.with_suffix('.npz'),**arrays);out.write_text(json.dumps(r,indent=2)+'\n');print(out)
if __name__=='__main__':main()
