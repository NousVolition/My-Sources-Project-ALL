"""New independent fields using the unchanged, previously validated fluid solver."""
import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'fluid'))
from recovery import run
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int);ap.add_argument('--label',default='confirmation');ap.add_argument('--n',type=int);ap.add_argument('--dt',type=float);args=ap.parse_args()
    if args.seed is not None:
        run(args.seed,n=args.n,dt=args.dt,label=args.label);raise SystemExit(0)
    p=json.loads((ROOT/'protocol.json').read_text())
    for seed in p['fluid_replication']['new_test_seeds']:run(seed,label='replication')
    # Two new-field numerical controls, declared before seeing forecast scores.
    for seed in p['fluid_replication']['new_test_seeds'][:2]:
        run(seed,dt=.01,label='half');run(seed,n=26,dt=.01,label='fine')
