import argparse,json
from pathlib import Path
from run_molecular import run

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);ap.add_argument('--stage',choices=['main','controls','temperature'],default='main');ap.add_argument('--platform',default='OpenCL');a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    cases=[]
    if a.stage=='main':
        for seed in (1,2,3):
            for isotope in ('H2O','D2O'):
                for beads in (1,32): cases.append(dict(n=64,isotope=isotope,beads=beads,T=300.,seed=seed,dt=.0005,eq=5.,production=5.,sample=.1,mobility=5.))
    elif a.stage=='controls':
        cases=[dict(n=n,isotope='H2O',beads=b,T=300.,seed=1,dt=dt,eq=5.,production=5.,sample=.1,mobility=5. if dt<.0005 else 0.) for n,b,dt in [(64,16,.0005),(64,64,.0005),(64,32,.00025),(216,32,.0005)]]
    else:
        cases=[dict(n=64,isotope=i,beads=b,T=T,seed=1,dt=.0005,eq=3.,production=3.,sample=.1,mobility=0.) for T in (280.,320.) for i in ('H2O','D2O') for b in (1,32)]
    (a.out/f'{a.stage}-plan.json').write_text(json.dumps(cases,indent=2))
    for c in cases:
        tag=f"n{c['n']}-{c['isotope']}-P{c['beads']}-T{c['T']}-s{c['seed']}-dt{c['dt']}"
        out=a.out/tag
        # No silent cache reuse or overwriting: invoke a new output root to rerun.
        print('START '+tag,flush=True);run(out,platform=a.platform,**c)

if __name__=='__main__': main()
