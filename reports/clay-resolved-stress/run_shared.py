"""Matched high-precision external drive; ordinary reciprocal coupling control."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,argparse
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'solver-versions/shared'))
from experiment import run,Parameters
ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--n',type=int,default=168);ap.add_argument('--layers',type=int,default=15);ap.add_argument('--two',action='store_true');ap.add_argument('--end',type=float,default=.5);ap.add_argument('--seconds',type=float,default=10800);ap.add_argument('--rtol',type=float,default=2e-7);ap.add_argument('--method',default='BDF');a=ap.parse_args()
run(a.name,Parameters(angles=a.n,layers=a.layers,feedback=a.two,alpha=128,rho=500,viscosity=.01875,preload=.35),end=a.end,max_seconds=a.seconds,rtol=a.rtol,atol=a.rtol/100,method=a.method)
