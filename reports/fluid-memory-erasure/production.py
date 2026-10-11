"""Production tasks after the documented numerical pilot; resume-safe."""
import argparse,json
from pathlib import Path
from run import ROOT,P,run_one,digest,dump

AMENDMENT=json.loads((ROOT/'production-amendment.json').read_text())
N=AMENDMENT['production_n']


def tasks(stage):
    if stage=='primary':return [dict(seed=s,n=N) for s in P['train_seeds']+P['validation_seeds']+P['test_seeds']]
    if stage=='refinement':
        return ([dict(seed=s,n=n) for s in P['test_seeds'] for n in P['grids'] if n!=N]+
                [dict(seed=s,n=N,dt=P['dt']/2) for s in P['test_seeds']]+
                [dict(seed=s,n=max(P['grids']),dt=P['dt']/2) for s in P['test_seeds'][:2]])
    if stage=='sensitivity':
        return ([dict(seed=s,n=N,end=max(P['long_horizons'])) for s in P['long_horizon_seeds']]+
                [dict(seed=s,n=N,amplitude=P['probe_velocity_rms']/2) for s in P['amplitude_control_seeds']]+
                [dict(seed=s,n=N,nu=nu) for s in P['viscosity_control_seeds'] for nu in P['viscosity_controls']])
    raise ValueError(stage)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['primary','refinement','sensitivity','all'],default='all')
    ap.add_argument('--backend',choices=['cpu','gpu'],default='gpu');ap.add_argument('--out',type=Path,default=ROOT/'data');args=ap.parse_args()
    audit=json.loads((ROOT/'pilot-results/numerical-audit.json').read_text())
    if not audit['all_passed']:raise RuntimeError('Pilot gates must pass before production')
    design={p:digest(ROOT/p) for p in ['protocol.json','production-amendment.json','production.py','solver.py','run.py']}
    plan=[]
    for stage in ['primary','refinement','sensitivity'] if args.stage=='all' else [args.stage]:
        planned=tasks(stage);complete=[]
        for task in planned:
            m=run_one(**task,backend=args.backend,out=args.out);complete.append(m['config']);plan.append(m['config'])
            dump(args.out/('production-execution-'+stage+'.json'),dict(stage=stage,planned=len(planned),completed=len(complete),configs=complete,design_sha256=design))
    dump(args.out/'production-execution.json',dict(stages=args.stage,completed_this_invocation=len(plan),configs=plan,design_sha256=design))

if __name__=='__main__':main()
