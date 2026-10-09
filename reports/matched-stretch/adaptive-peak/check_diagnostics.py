import numpy as np
from audit import ROOT,Flow,detailed,gate,save,L
f=Flow(32,workers=1);axis=np.arange(32)*f.dx-L/2
x,y,z=np.meshgrid(axis,axis,axis,indexing='ij');k=2*np.pi/L
u=np.stack([np.sin(k*y),np.zeros_like(x),np.zeros_like(x)]);h=f.project(f.hat(u))
r=detailed(f,h,0.)
assert abs(r['Wmax']-k)<1e-12
assert abs(r['energy']-L**3/4)<1e-10
assert abs(r['enstrophy']-L**3*k*k/4)<1e-10
assert abs(sum(r['spectra']['energy_sum'])-r['energy'])<1e-10
assert abs(sum(r['spectra']['enstrophy_sum'])-r['enstrophy'])<1e-10
assert abs(r['peak_budget']['transport_LHS']['along_vorticity'])<1e-12
assert abs(r['peak_budget']['stretching_RHS']['along_vorticity'])<1e-12
assert abs(r['peak_budget']['viscosity_RHS']['along_vorticity']+f.nu*k**3)<1e-12
assert abs(r['alignment_degrees']-90)<1e-10
assert r['peak_budget']['closure_relative_l2']<1e-10
assert gate(r)['pass']
save(ROOT/'diagnostics-verification.json',dict(status='passed',analytic_test='Single sine shear: exact W, energy, enstrophy, zero stretching/transport, viscous decay, 90-degree strain alignment, Parseval and a passing resolved-width gate.',budget_closure=r['peak_budget']['closure_relative_l2']))
print('Analytic shear diagnostics passed')
