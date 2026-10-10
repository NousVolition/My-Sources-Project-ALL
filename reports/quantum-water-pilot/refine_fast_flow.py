"""Bounded follow-up to the failed 100 mm/s spatial convergence screen.

Preserve the original 24/32/48 results. Add 64 and 96 points per axis and
a half-step at 64. A small grid difference is evidence, not an exact-error bound.
"""
import argparse,time
from pathlib import Path
import numpy as np
from run_fluid import properties,dump
from vendor.ns_solver import Solver,difference


def run(root):
    out=root/'fluid-refinement';out.mkdir(exist_ok=False)
    plan=dict(temperature_K=298.15,U_ms=.1,end_seconds=.05,
              cases=[dict(n=64,dt=.01),dict(n=96,dt=.01),dict(n=64,dt=.005)],
              spatial_screen=.01,time_screen=.001,
              reason='Follow up failed 32-to-48 grid comparison, without replacing it.',
              limit='Stop at 96; report an unresolved result if this still fails. No exact error bound or experimental validation is inferred.')
    dump(out/'protocol.json',plan);rows=[];checks=[]
    for iso,p in properties(plan['temperature_K']).items():
        fields={48:np.load(root/'fluid-stress'/f'{iso}-U0.1-n48.npz')['final_hat']}
        for case in plan['cases']:
            n=case['n'];dt=case['dt'];s=Solver(n,p['nu']/(.1*.001),workers=1)
            h=s.initial('taylor_green');h0=h.copy();e0=.5*s.inner(h,h)
            loss=0.;cfl=0.;steps=round(5/dt);tick=time.perf_counter()
            for j in range(steps):
                h,dl,c=s.step(h,5/steps);loss+=dl;cfl=max(cfl,c)
                if (j+1)%100==0:print(iso,n,dt,j+1,'/',steps,round(time.perf_counter()-tick,1),'s',flush=True)
            obs,_,_=s.observe(h,5.,loss,e0)
            rows.append(dict(isotope=iso,n=n,dt=dt,cfl=cfl,seconds=time.perf_counter()-tick,**obs))
            np.savez_compressed(out/f'{iso}-n{n}-dt{dt}.npz',final_hat=h,initial_hat=h0)
            if dt==.01:
                coarse=48 if n==64 else 64
                change=difference(fields[coarse],h,Solver(coarse,s.nu,workers=1))
                checks.append(dict(isotope=iso,kind='grid',coarse_n=coarse,fine_n=n,relative_field_change=change,passed=change<plan['spatial_screen']))
                fields[n]=h
            else:
                # Direct subtraction avoids cancellation when temporal error is tiny.
                change=np.sqrt(s.inner(h-fields[n],h-fields[n])/s.inner(h,h))
                checks.append(dict(isotope=iso,kind='time',n=n,coarse_dt=.01,fine_dt=dt,relative_field_change=change,passed=change<plan['time_screen']))
            dump(out/'results.json',rows);dump(out/'convergence.json',checks)
            print('FINISHED',rows[-1],checks[-1],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True)
    run(ap.parse_args().data)
