"""The uncoupled Lorenz input is independent of every material discretization."""
import numpy as np
from scipy.integrate import solve_ivp
from types import SimpleNamespace
def rhs(t,y,rho):
 X,Y,Z=y;return [10*(Y-X),X*(rho-Z)-Y,X*Y-(8/3)*Z]
def reference(end,initial,rho,method='DOP853'):
 # Fixed half-unit segments preserve all prior input values when extending.
 # In particular [0,.5] is exactly the published reference calculation.
 horizon=max(.5,np.ceil(end/.5)*.5);parts=[];y=np.asarray(initial,dtype=float)
 for a in np.arange(0,horizon,.5):
  z=solve_ivp(lambda t,y:rhs(t,y,rho),(a,a+.5),y,method=method,rtol=1e-12,atol=1e-14,max_step=.002,dense_output=True)
  if not z.success:raise RuntimeError(z.message)
  parts.append(z);y=z.y[:,-1]
 def evaluate(t):
  ts=np.asarray(t);flat=ts.reshape(-1);ids=np.maximum(0,np.ceil(flat/.5).astype(int)-1);ids=np.minimum(ids,len(parts)-1)
  if np.any(flat<0) or np.any(flat>horizon):raise ValueError('Outside reference horizon')
  out=np.empty((3,len(flat)))
  for i in np.unique(ids):out[:,ids==i]=parts[i].sol(flat[ids==i])
  return out[:,0] if ts.ndim==0 else out.reshape((3,)+ts.shape)
 return SimpleNamespace(sol=evaluate,nfev=sum(z.nfev for z in parts),success=True,horizon=horizon)
