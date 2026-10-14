"""Accept an actual finer recording only after full matched-history checks."""
from pathlib import Path
import argparse,json,hashlib
from analyze import load,compare
ROOT=Path(__file__).resolve().parent
def select(paths):
 paths=[Path(p).resolve() for p in paths]
 checks=[]
 for coarse,fine in zip(paths,paths[1:]):
  c=compare(load(coarse),load(fine));c['inputs']=[dict(metadata=hashlib.sha256(Path(p).read_bytes()).hexdigest(),raw=hashlib.sha256(Path(p).with_suffix('.npz').read_bytes()).hexdigest()) for p in [coarse,fine]];checks.append(c)
  if c['full_mesh_check_passed']:
   r=json.loads(Path(fine).read_text());return dict(accepted=True,selected=Path(fine).relative_to(ROOT).as_posix(),mesh=[r['parameters']['angles'],r['parameters']['layers']],checks=checks)
 return dict(accepted=False,selected=None,mesh=None,checks=checks)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('paths',nargs='+',type=Path);ap.add_argument('--output',type=Path);a=ap.parse_args();result=select(a.paths)
 if a.output:assert not a.output.exists();a.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2));raise SystemExit(0 if result['accepted'] else 1)
if __name__=='__main__':main()
