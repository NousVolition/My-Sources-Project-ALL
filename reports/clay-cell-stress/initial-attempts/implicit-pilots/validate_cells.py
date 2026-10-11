"""Independent energy, gradient, geometry and operational-limit checks."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from cell_hug import Parameters,Mesh,geometry,swept_shape
from contact_hug import Mesh as OldMesh
from implicit_solver import integrate

def main():
    checks=[];rng=np.random.default_rng(51)
    def check(name,value,limit):
        checks.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(value<=limit)))
    p=Parameters();m=Mesh(p);old=OldMesh(p)
    # Active compression and active distortion at two smooth affine configurations.
    for label,F in [('compression',np.diag([.67,.67])),('shape',np.diag([1.5,.7]))]:
        x=m.initial@F.T;rest=m.l0.copy();e,g,*_=m.material(x,rest)
        for j in range(5):
            direction=rng.normal(size=x.shape);direction/=np.linalg.norm(direction);h=1e-6
            fd=(m.material(x+h*direction,rest)[0]-m.material(x-h*direction,rest)[0])/(2*h)
            analytic=float(np.sum(g*direction))
            check(label+'-gradient-'+str(j),abs(fd-analytic)/max(1,abs(fd)),2e-7)
        added=g-old.material(x,rest)[1]
        check(label+'-net-force',np.linalg.norm(added.sum(axis=0)),1e-10)
        check(label+'-net-torque',abs(np.sum(x[:,0]*added[:,1]-x[:,1]*added[:,0])),1e-10)
    check('inactive-energy',abs(m.material(m.initial,m.l0)[0]-old.material(m.initial,m.l0)[0]),1e-14)
    for label,F in [('area-domain',np.diag([.4,.4])),('shape-domain',np.diag([2.3,.2]))]:
        caught=False
        try:m.material(m.initial@F.T,m.l0)
        except FloatingPointError:caught=True
        check(label,0 if caught else 1,0)
    for j in range(4):
        a=m.pack();b=a.copy();b[:m.position_size]+=.025*rng.normal(size=m.position_size)
        bounds=swept_shape(m,a,b)
        dense=np.min([[geometry(m,m.unpack((1-t)*a+t*b)[0])[k] for k in ['min_area_ratio','min_quality_ratio','min_altitude_ratio']] for t in np.linspace(0,1,101)],axis=0)
        check('swept-lower-bound-'+str(j),max(0,float(np.max(np.array(bounds)-dense))),1e-12)
    # Energy derivatives and cancellation of the reciprocal port with barriers active.
    p=replace(p,area_activation=.95,shape_activation=.95);m=Mesh(p)
    y=m.pack();y[:m.position_size]*=.9;y[-7:-4]=[2,3,4]
    f=m.rhs(0,y);h=1e-7
    def energy(v):
        x,r,l=m.unpack(v);return m.energy_force(x,r)[0]+.5*p.lorenz_capacity*np.dot(l,l)
    fd=(energy(y+h*f)-energy(y-h*f))/(2*h)
    expected=f[-4]-f[-3]-f[-2]
    check('total-work-identity',abs(fd-expected)/max(1,abs(expected)),2e-7)
    result=integrate(p,end=.01,max_rhs=2)
    check('budget-failure-preserved',0 if result[3] and result[1][-1]==0 else 1,0)
    out=dict(all_passed=all(c['passed'] for c in checks),checks=checks)
    path=Path(__file__).resolve().parent/'data/validation.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2));assert out['all_passed']
if __name__=='__main__':main()
