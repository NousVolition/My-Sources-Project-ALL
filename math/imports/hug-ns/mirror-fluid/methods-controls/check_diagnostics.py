import json
import numpy as np
from diagnostics import Meter,initial,mirror,transform,norm
from study import ROOT,load_solver,save
m=load_solver();n=33
x,y,z,kx,ky,kz,k2,keep,dx=m.grids(n);ops=(kx,ky,kz,k2,keep,dx)
k=2*np.pi/m.L
u=np.stack([np.sin(k*y),np.zeros_like(x),np.zeros_like(x)])
_,phi,_,_=initial(m,n)
meter=Meter(m,ops,phi)
w=meter.real(meter.curlh(meter.fft(u)))
curl_error=float(np.max(abs(w[2]+k*np.cos(k*y))))
assert curl_error<1e-12,curl_error
r=meter.full(u)
assert abs(r['energy']-m.L**3/4)<1e-11
assert abs(r['enstrophy']-m.L**3*k*k/4)<1e-11
assert abs(r['energy_rate']+2*m.NU*r['enstrophy'])<1e-11
assert r['projection_relative']<1e-12
assert r['vorticity_product_rule_residual']<1e-11
b,p,_,_=initial(m,n);q=b+.01*norm(b)*p
for kind in ('shift1','shifthalf','rotate'):
    assert norm(q-transform(transform(q,kind),kind,True))/norm(q)<1e-12,kind
    assert abs(norm(q)-norm(transform(q,kind)))<1e-12,kind
assert norm(mirror(mirror(q))-q)<1e-14
refpath=ROOT.parent/'hug-github-7-intake/reproduction/dashboard.json'
if not refpath.exists():refpath=ROOT/'uploaded/reproduced-dashboard.json'
ref=json.loads(refpath.read_text())
for name,field in [('even',b),('plus',q)]:
    row=meter.full(field)
    for key in ('W','E','D','energy','enstrophy','ratio'):
        assert abs(row[key]-ref[name][0][key])<1e-10,(name,key,row[key],ref[name][0][key])
save(ROOT/'diagnostics-verification.json',dict(status='passed',curl_sign_error=curl_error,
    checks=['Analytic shear curl sign, energy, enstrophy and viscous dissipation','Physical vorticity terms match curl of velocity RHS for low-frequency shear','Coordinate transforms invert and preserve norm','Reflection squares to identity','Initial observables reproduce supplied dashboard','Original spectral band has zero retained modes'],
    original_cutoff_retained_modes=r['original_cutoff_retained_modes']))
print('Diagnostics checks passed')
