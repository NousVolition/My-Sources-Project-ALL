"""Compute real common-input meshes and refuse unresolved accuracy at the cap."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
from select_mesh import select
ROOT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--evidence',nargs='+',type=Path);ap.add_argument('--start',type=int,default=144);ap.add_argument('--max-mesh',type=int,default=168);ap.add_argument('--end',type=float,default=.5);ap.add_argument('--seconds',type=float,default=10800);a=ap.parse_args()
 out=ROOT/'data'/(a.name+'-matched-decision.json');assert not out.exists()
 if a.evidence:paths=a.evidence;result=select(paths)
 else:
  assert a.start>=24 and a.start%24==0 and a.max_mesh>=a.start+24 and a.max_mesh%24==0
  paths=[];result=dict(accepted=False,selected=None,mesh=None,checks=[])
  for n in range(a.start,a.max_mesh+1,24):
   name=f'{a.name}-{n}x{n//12+1}'
   subprocess.run([sys.executable,str(ROOT/'run_shared.py'),name,'--n',str(n),'--layers',str(n//12+1),'--end',str(a.end),'--seconds',str(a.seconds)],check=True)
   paths.append(ROOT/'solver-versions/shared/data'/(name+'.json'))
   if len(paths)>1:result=select(paths)
   if result['accepted']:break
 result.update(driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),record_paths=[p.resolve().relative_to(ROOT).as_posix() for p in paths])
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));raise SystemExit(0 if result['accepted'] else 1)
if __name__=='__main__':main()
