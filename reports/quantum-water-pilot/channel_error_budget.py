"""Post-primary diagnostic: separate spatial and temporal channel errors.

Use exact eigenmodes of the finite-difference Laplacian. These are analytical
evaluations of the recorded scheme, not additional physical simulation runs.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.fft import dst,idst
from run_fluid import properties,exact_startup,dump


def run(root):
    saved=json.loads((root/'fluid-data/channel_convergence.json').read_text())
    rows=[];n=129;H=.001;G=10.;end=.1;y=np.linspace(0,H,n);dy=y[1]-y[0]
    for iso,p in properties(298.15).items():
        k=np.arange(1,n-1)
        lam=4*p['nu']/dy**2*np.sin(k*np.pi/(2*(n-1)))**2
        steady=G*y*(H-y)/(2*p['mu']);a=dst(steady[1:-1],type=1,norm='ortho')
        discrete=steady.copy();discrete[1:-1]-=idst(a*np.exp(-lam*end),type=1,norm='ortho')
        continuum=exact_startup(y,end,p['rho'],p['mu'],G,H)
        norm=np.linalg.norm(continuum);space=discrete-continuum
        for dt in [.001,.0005,.00025]:
            amplification=(1-lam*dt/2)/(1+lam*dt/2)
            cn=steady.copy();cn[1:-1]-=idst(a*amplification**round(end/dt),type=1,norm='ortho')
            temporal=cn-discrete;total=cn-continuum
            measured=next(r['relative_L2_error'] for r in saved if r['material']==iso and r['n']==n and r['dt_s']==dt)
            relative=np.linalg.norm(total)/norm
            assert abs(relative-measured)<1e-11
            rows.append(dict(isotope=iso,n=n,dt_s=dt,
                spatial_error_relative_to_continuum=np.linalg.norm(space)/norm,
                temporal_error_relative_to_continuum=np.linalg.norm(temporal)/norm,
                total_error_relative_to_continuum=relative,
                error_direction_cosine=np.dot(space,temporal)/(np.linalg.norm(space)*np.linalg.norm(temporal)),
                reconstructed_vs_recorded_error_difference=abs(relative-measured)))
    result=dict(purpose='Explain error cancellation in the existing channel time-step sweep.',
                method='Exact discrete-Laplacian sine modes separate finite-grid and Crank-Nicolson errors at 0.1 s. Both components use the continuum norm as denominator.',
                comparisons=rows)
    dump(root/'channel-error-budget.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True)
    run(ap.parse_args().data)
