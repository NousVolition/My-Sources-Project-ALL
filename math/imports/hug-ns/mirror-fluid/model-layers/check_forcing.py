"""Analytic forced-shear check of the separate driven Heun implementation."""
import hashlib,json
from pathlib import Path
import numpy as np
from run_fluid import ROOT,solver,driven_step
m=solver();n=17;x,y,z,*ops0=m.grids(n);ops=tuple(ops0);k=2*np.pi/m.L
shape=np.array([np.zeros_like(x),np.sin(k*x),np.zeros_like(x)])
a0=.3;B=.2;P=.2;T=.4;lam=m.NU*k*k;omega=2*np.pi/P
exact=a0*np.exp(-lam*T)+B*(lam*np.sin(omega*T)-omega*np.cos(omega*T)+omega*np.exp(-lam*T))/(lam*lam+omega*omega)
errors=[]
for dt in (.01,.005,.0025):
 u=a0*shape
 for step in range(round(T/dt)):u,_=driven_step(m,u,dt,step*dt,ops,B*shape,P)
 amp=np.mean(np.sum(u*shape,axis=0))/np.mean(np.sum(shape*shape,axis=0));error=abs(amp-exact);errors.append(dict(dt=dt,error=float(error)))
assert errors[0]['error']>3.5*errors[1]['error'] and errors[1]['error']>3.5*errors[2]['error']
v=dict(status='passed',exact_amplitude=float(exact),errors=errors,driver_sha256=hashlib.sha256((ROOT/'run_fluid.py').read_bytes()).hexdigest())
(ROOT/'forcing-verification.json').write_bytes(json.dumps(v,indent=2).encode());print(json.dumps(v))
