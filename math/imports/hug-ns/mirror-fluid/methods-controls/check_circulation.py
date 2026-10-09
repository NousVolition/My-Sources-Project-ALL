import numpy as np
from audit_saved import circulation,source_audit
from study import ROOT,load_solver,save
m=load_solver();n=33;x,y,z,*ops=m.grids(n);k=2*np.pi/m.L
u=np.stack([np.sin(k*y),np.zeros_like(x),np.zeros_like(x)])
actual=circulation(u,m,tuple(ops),.9)
expected=.6*(np.sin(k*(-.3))-np.sin(k*.3))
assert abs(actual-expected)<1e-12,(actual,expected)
source_audit(m)
save(ROOT/'circulation-verification.json',dict(status='passed',analytic_shear_loop=expected,measured=actual,
 correction='Fixed missing N factor in FFT frequencies before publication. Recomputed source initialization audit; no fluid field or trajectory changed.'))
print(actual,expected)
