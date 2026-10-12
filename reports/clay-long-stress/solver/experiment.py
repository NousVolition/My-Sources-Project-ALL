"""Record all attempts, actual states, parameters, source hashes and work audits."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,hashlib,time
from pathlib import Path
from dataclasses import asdict
import numpy as np
from mesh_hug import Parameters,Mesh,object_outline
from implicit_solver import integrate
ROOT=Path(__file__).resolve().parent
def hashes():return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py') if p.name in ['clay_hug.py','contact_hug.py','cell_hug.py','mesh_hug.py','implicit_solver.py','jacobian.py','analytic_spatial.py','numerical_jacobian.py','experiment.py','self_contact.py']}
def run(name,p,end=.2,rtol=2e-9,atol=2e-11,method='BDF',max_step=.01,max_seconds=600,initial=(-8,8,27),state=None,start_time=0.,parent=None):
    source=hashes();path=ROOT/'data'/name;path.parent.mkdir(exist_ok=True)
    if path.with_suffix('.json').exists() or path.with_suffix('.npz').exists():raise ValueError('Never overwrite an attempt: '+name)
    if state is not None and parent is None:raise ValueError('Continued states require source provenance')
    m,t,y,f,a=integrate(p,end=end,rtol=rtol,atol=atol,method=method,max_step=max_step,max_seconds=max_seconds,initial=initial,state=state,record_step=.01,start_time=start_time)
    assert source==hashes(),'Sources changed during integration'
    d=[m.diagnostics(v) for v in y]
    E=np.array([v['material_energy']+v['lorenz_energy'] for v in d]);W=y[:,-4]-y[:,-3]-y[:,-2]+(0 if p.feedback else y[:,-1])
    balance=float(np.max(abs(E-E[0]-(W-W[0])))/max(1,np.max(abs(E)),np.max(abs(W-W[0]))))
    result=dict(name=name,parameters=asdict(p),requested_end=end,completed_end=float(t[-1]),initial=list(initial),start_time=float(t[0]),parent=parent,failure=f,audit=a,source_hashes=source,energy_balance_scaled=balance,diagnostics=d)
    np.savez_compressed(path.with_suffix('.npz'),time=t,state=y,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p))
    path.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(name=name,end=float(t[-1]),seconds=a['seconds'],failure=f,gap=d[-1]['surface_gap_rms'],arcgap=d[-1]['gap_arclength_rms'],imprint=d[-1]['imprint'],energy=balance)),flush=True)
    return result,y[-1]
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--n',type=int,default=48);ap.add_argument('--layers',type=int,default=5);ap.add_argument('--two',action='store_true');ap.add_argument('--end',type=float,default=.2);ap.add_argument('--alpha',type=float,default=16);ap.add_argument('--rho',type=float,default=160);ap.add_argument('--viscosity',type=float,default=.15);ap.add_argument('--old-contact',action='store_true');ap.add_argument('--old-shape',action='store_true');ap.add_argument('--rtol',type=float,default=2e-9);ap.add_argument('--method',default='BDF');ap.add_argument('--seconds',type=float,default=600);ap.add_argument('--quadrature',type=int,default=32)
    a=ap.parse_args();p=Parameters(angles=a.n,layers=a.layers,alpha=a.alpha,rho=a.rho,viscosity=a.viscosity,preload=.35,feedback=a.two,surface_contact=not a.old_contact,continuum_shape=not a.old_shape,quadrature=a.quadrature)
    run(a.name,p,end=a.end,rtol=a.rtol,atol=a.rtol/100,method=a.method,max_seconds=a.seconds)
