"""Compute and select a real finer mesh; report failure at the declared cap."""
from pathlib import Path
import sys,argparse,json,hashlib
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from experiment import run,Parameters
from select_mesh import select
def main():
 ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--evidence',nargs='+',type=Path);ap.add_argument('--start',type=int,default=72);ap.add_argument('--max-mesh',type=int,default=120);ap.add_argument('--two',action='store_true');ap.add_argument('--alpha',type=float,default=128);ap.add_argument('--rho',type=float,default=500);ap.add_argument('--viscosity',type=float,default=.01875);ap.add_argument('--end',type=float,default=.5);ap.add_argument('--seconds',type=float,default=10800);args=ap.parse_args()
 out=ROOT/'data'/(args.name+'-resolution.json');assert not out.exists()
 if args.evidence:result=select(args.evidence)
 else:
  assert args.start>=24 and args.start%24==0 and args.max_mesh>=args.start+24 and args.max_mesh%24==0
  paths=[];result=dict(accepted=False,selected=None,mesh=None,checks=[])
  for n in range(args.start,args.max_mesh+1,24):
   name=f'{args.name}-{n}x{n//12+1}';run(name,Parameters(angles=n,layers=n//12+1,feedback=args.two,alpha=args.alpha,rho=args.rho,viscosity=args.viscosity,preload=.35),end=args.end,max_seconds=args.seconds)
   paths.append(ROOT/'solver/data'/(name+'.json'))
   if len(paths)>1:result=select(paths)
   if result['accepted']:break
 result['driver_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));raise SystemExit(0 if result['accepted'] else 1)
if __name__=='__main__':main()
