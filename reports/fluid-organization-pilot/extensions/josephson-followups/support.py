"""Shared integration and run-group inference; no dependence on prior folders."""
from pathlib import Path
import json,csv
import numpy as np
from scipy.integrate import solve_ivp
from circuit import stable_rhs
ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'protocol.json').read_text(encoding='utf-8'))
def save(name,obj):
    path=ROOT/name;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(name,rows):
    p=ROOT/'results'/name;p.parent.mkdir(exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def interval(x,seed=43105):
    x=np.asarray(x,float);r=np.random.default_rng(seed)
    b=x[r.integers(0,len(x),(4000,len(x)))].mean(1)
    return {'mean':float(x.mean()),'ci95':np.quantile(b,[.025,.975]).tolist(),'n':len(x)}
def evolve(s,bias,alpha,duration,times=None,max_step=.5,rtol=1e-10,atol=1e-12):
    s=np.asarray(s,float);shape=s.shape
    def fun(t,y):return stable_rhs(y.reshape(shape),bias(t) if callable(bias) else bias,alpha).ravel()
    sol=solve_ivp(fun,(0,duration),s.ravel(),method='DOP853',rtol=rtol,atol=atol,max_step=max_step,t_eval=times,dense_output=times is None)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y.T.reshape((-1,)+shape) if times is not None else sol
