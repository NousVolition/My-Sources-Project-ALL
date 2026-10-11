"""Run both directions on two meshes and check time integration separately.

The default starts at 48 x 5, not the known-underresolved 24 x 3 preview.
An agreement flag is never issued for incomplete or unchecked comparisons.
"""
import sys,argparse,json
from pathlib import Path
from dataclasses import replace
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver'))
from experiment import run,Parameters
from accuracy_gate import assess
def main():
    ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--end',type=float,default=.2);ap.add_argument('--alpha',type=float,default=16);ap.add_argument('--rho',type=float,default=160);ap.add_argument('--viscosity',type=float,default=.15);ap.add_argument('--budget',type=float,default=1200);a=ap.parse_args()
    base=Parameters(alpha=a.alpha,rho=a.rho,viscosity=a.viscosity,preload=.35,quadrature=32)
    checks=[]
    for mode,fb in [('one',False),('two',True)]:
        names=[]
        for n,l in [(48,5),(72,7)]:
            name=f'{a.name}-{mode}-{n}';run(name,replace(base,angles=n,layers=l,feedback=fb),end=a.end,max_seconds=a.budget);names.append(ROOT/'solver/data'/(name+'.json'))
        checks.append(assess(*names,'mesh'))
        name=f'{a.name}-{mode}-radau';run(name,replace(base,angles=48,layers=5,feedback=fb),end=a.end,method='Radau',rtol=2e-10,atol=2e-12,max_seconds=a.budget)
        checks.append(assess(names[0],ROOT/'solver/data'/(name+'.json'),'time'))
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,scope='Agreement over saved times; not a bound against an exact solution or a physical clay calibration.')
    (ROOT/'solver/data'/(a.name+'-verification.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
if __name__=='__main__':main()
