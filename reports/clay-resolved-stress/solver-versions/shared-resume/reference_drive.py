"""The uncoupled Lorenz input is independent of every material discretization."""
import numpy as np
from scipy.integrate import solve_ivp
def rhs(t,y,rho):
 X,Y,Z=y;return [10*(Y-X),X*(rho-Z)-Y,X*Y-(8/3)*Z]
def reference(end,initial,rho,method='DOP853'):
 sol=solve_ivp(lambda t,y:rhs(t,y,rho),(0,end),np.asarray(initial,dtype=float),method=method,rtol=1e-12,atol=1e-14,max_step=.002,dense_output=True)
 if not sol.success:raise RuntimeError(sol.message)
 return sol
