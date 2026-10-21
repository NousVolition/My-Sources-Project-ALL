from pathlib import Path
import json,time
import numpy as np
from model import Model,cases
R=Path(__file__).resolve().parent
p=json.loads((R/'protocol.json').read_text());cs=cases(p);checks=[]
def check(name,value,limit):
    checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
    if value>limit:raise AssertionError(checks[-1])
m=Model(24,p,cs);h=m.start();f=m.real(h)
# Nonuniform guards exercise all variable-viscosity derivative terms.
f[:,4]+=.06*np.cos(m.x)*np.sin(m.y);f[:,5]+=.04*np.sin(2*m.x+m.y)
h=m.project(m.fft(f));rhs,aux=m.rhs(.7,h);other=m.independent_rhs(.7,h)
check('independent_curl_and_conservative_rhs',np.max(np.abs(rhs-other)),1e-12)
fr=m.real(rhs);diag=m.diagnostics(h)
check('semidiscrete_kinetic_energy_identity',np.max(np.abs(np.mean(f[:,0]*fr[:,0]+f[:,1]*fr[:,1],axis=(-2,-1))+diag['dissipation'])),1e-12)
check('momentum_derivative',np.max(np.abs(fr[:,:2].mean(axis=(-2,-1)))),1e-12)
fine=Model(48,p,cs)
check('identical_initial_fourier_fields',np.max(np.abs(m.real(m.start())-fine.real(fine.start())[:,:,::2,::2])),1e-12)
check('identical_disturbance_profile',np.max(np.abs(m.profile-fine.profile[::2,::2])),1e-12)
# Analytic variable stress check: u=sin(y), v=0, gAB=.5+.1cos(x).
one=[dict(id=0,mode='directional',weight='equal',k=[1,1],pulse='none')];q=Model(24,p,one)
f=q.real(q.start());f[:,0]=np.sin(q.y);f[:,1]=0;f[:,4]=.5+.1*np.cos(q.x);f[:,5]=.5
hh=q.project(q.fft(f));got=q.rhs(0,hh)[0]
nu=p['mu']+p['a']*(1+.1*np.cos(q.x));nux=-.1*p['a']*np.sin(q.x)
expect=np.zeros_like(f);expect[:,0]=-nu*np.sin(q.y);expect[:,1]=nux*np.cos(q.y)
check('full_variable_stress_gradient_terms',np.max(np.abs(got[:,:2]-q.project(q.fft(expect))[:,:2])),1e-12)
# Exact decaying shear for constant viscosity, while SIMS and guards can evolve.
one[0]['mode']='no_feedback';q=Model(24,p,one);f=q.real(q.start());f[:,0]=np.sin(q.y);f[:,1]=0;hh=q.project(q.fft(f));dt=.01
for j in range(100):
    t=j*dt;a=q.rhs(t,hh)[0];b=q.rhs(t+dt/2,hh+dt*a/2)[0];c=q.rhs(t+dt/2,hh+dt*b/2)[0];d=q.rhs(t+dt,hh+dt*c)[0]
    hh=q.project(hh+dt*(a+2*b+2*c+d)/6)
check('analytic_shear_decay',np.max(np.abs(q.real(hh)[0,0]-np.sin(q.y)*np.exp(-(p['mu']+p['a'])))),1e-10)
(R/'preflight.json').write_text(json.dumps(dict(passed=all(c['passed'] for c in checks),checks=checks),indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
