"""Additional zero-port and unmodified SciPy finite-Jacobian controls."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib
from pathlib import Path
from dataclasses import replace
import numpy as np
from run_experiments import run,Parameters,hashes
import implicit_solver
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def main():
    protocol=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model_sources=hashes(),zero_port_end=2,independent_end=.08,independent_method='Radau',rtol=2e-9,atol=2e-11,max_step=.005,spatial_jacobian='SciPy numerical Jacobian, independently from custom blocks')
    (DATA/'independent-protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for fb in [False,True]:run('zero-port-'+('two' if fb else 'one'),Parameters(preload=.35,alpha=0,feedback=fb),end=2,group='zero-port')
    # Strip the supplied Jacobian before calling the independent SciPy class.
    # Continuous forces, domain guards and output targets remain identical.
    original=implicit_solver.Radau
    class NumericRadau(original):
        def __init__(self,*args,**kw):
            kw['jac']=None;super().__init__(*args,**kw)
    implicit_solver.Radau=NumericRadau
    try:
        for label,p in [('prior',Parameters(preload=.35,alpha=16,rho=160,viscosity=.15,feedback=False)),('hardest',Parameters(preload=.35,alpha=64,rho=300,viscosity=.0375,feedback=False))]:
            run('numeric-jac-'+label,p,end=.08,method='Radau',max_step=.005,rtol=2e-9,atol=2e-11,group='independent')
    finally:implicit_solver.Radau=original
if __name__=='__main__':main()
