import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np,json
from pathlib import Path
from dataclasses import replace
from mesh_hug import Mesh,Parameters,swept_shape
from jacobian import coupled_jacobian
def main():
    checks=[];rng=np.random.default_rng(42)
    def check(name,error,tol):checks.append(dict(name=name,error=float(error),tolerance=tol,passed=bool(error<=tol)))
    p=Parameters();m=Mesh(p)
    check('identity-deformation',np.max(abs(m.deformation(m.initial)-np.eye(2))),1e-13)
    for label,F in [('compressed',np.diag([.67,.67])),('distorted',np.diag([1.5,.7]))]:
        x=m.initial@F.T;e,g,*_=m.material(x,m.l0)
        for j in range(4):
            d=rng.normal(size=x.shape);d/=np.linalg.norm(d);h=1e-6
            fd=(m.material(x+h*d,m.l0)[0]-m.material(x-h*d,m.l0)[0])/(2*h)
            check(label+str(j),abs(fd-np.sum(g*d))/max(1,abs(fd)),3e-7)
    # Contact is active near the tips of the ellipse.
    y=m.pack((2,3,4));y[:m.position_size]*=.86
    e,g,*_=m.energy_force(*m.unpack(y)[:2])
    for j in range(4):
        d=rng.normal(size=(m.n,2));d/=np.linalg.norm(d);h=1e-7;x,r,_=m.unpack(y)
        fd=(m.energy_force(x+h*d,r)[0]-m.energy_force(x-h*d,r)[0])/(2*h)
        check('contact-gradient'+str(j),abs(fd-np.sum(g*d))/max(1,abs(fd)),2e-6)
    for fb in [False,True]:
        m=Mesh(replace(p,feedback=fb));f=m.rhs(0,y);h=1e-8
        def E(v):
            x,r,l=m.unpack(v);return m.energy_force(x,r)[0]+.5*m.p.lorenz_capacity*np.dot(l,l)
        fd=(E(y+h*f)-E(y-h*f))/(2*h);expected=f[-4]-f[-3]-f[-2]+(0 if fb else f[-1])
        check('work-identity-'+str(fb),abs(fd-expected)/max(1,abs(expected)),3e-6)
        J=coupled_jacobian(m,m.rhs,0,y)
        for j in range(3):
            d=rng.normal(size=y.shape);d/=np.linalg.norm(d);h=1e-7
            fd=(m.rhs(0,y+h*d)-m.rhs(0,y-h*d))/(2*h)
            check('jacobian-'+str(fb)+'-'+str(j),np.max(abs(fd-J@d))/max(1,np.max(abs(fd))),5e-5)
    out=dict(all_passed=all(c['passed'] for c in checks),checks=checks)
    (Path(__file__).parent/'data/validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));assert out['all_passed']
if __name__=='__main__':main()
