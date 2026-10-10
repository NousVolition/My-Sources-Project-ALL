"""Check the strongest vortex initial velocity at every swept temperature."""
import argparse,json
from pathlib import Path
import numpy as np
from run_fluid import properties,vortex,dump
from vendor.ns_solver import Solver,difference

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',required=True,type=Path);a=ap.parse_args();checks=[]
    for T in [283.15,298.15,313.15]:
        for isotope,p in properties(T).items():
            h16=np.load(a.data/f'vortex-T{T}-{isotope}-U0.02-n16.npz')['final_hat']
            r24,h24,i24,s24=vortex(p,U=.02,n=24)
            r32,h32,i32,s32=vortex(p,U=.02,n=32)
            rhalf,half,_,_=vortex(p,U=.02,n=24,dt=.005)
            np.savez_compressed(a.data/f'extended-T{T}-{isotope}.npz',n24=h24,n32=h32,n24_half_dt=half,initial24=i24,initial32=i32)
            checks.append(dict(T_K=T,isotope=isotope,n16_n24=difference(h16,h24,Solver(16,p['nu']/(.02*.001))),n24_n32=difference(h24,h32,s24),half_dt=np.sqrt(s24.inner(half-h24,half-h24)/s24.inner(half,half)),fine_result=r32))
    dump(a.data/'extended_convergence.json',checks)
    print(checks,flush=True)

if __name__=='__main__':main()
