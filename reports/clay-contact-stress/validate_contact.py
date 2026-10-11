from pathlib import Path
from dataclasses import replace
import json
import numpy as np
from contact_hug import *
ROOT=Path(__file__).resolve().parent

def main():
    rng=np.random.default_rng(20261011);rows=[]
    def check(name,error,target):rows.append(dict(name=name,error=float(error),target=target,passed=bool(error<target)))
    for shape in ['circle','ellipse','lobed']:
        p=Parameters(object_shape=shape,preload=.35);m=Mesh(p)
        # Bring the annulus close enough to activate a barrier while staying admissible.
        x=m.initial.copy();scale=1.
        while True:
            trial=x*.998
            try:g=m.clearance(trial,False)
            except FloatingPointError:break
            if min(g)<.012:break
            x=trial
        rest=m.l0*.98;direction=rng.normal(size=x.shape);direction/=np.linalg.norm(direction);h=1e-7
        e,grad,*_=m.energy_force(x,rest)
        fd=(m.energy_force(x+h*direction,rest)[0]-m.energy_force(x-h*direction,rest)[0])/(2*h)
        check('barrier-energy-gradient-'+shape,abs(fd-np.sum(grad*direction))/max(1,abs(fd)),2e-6)
        y=m.pack();y[:m.position_size]=x.ravel();y[m.position_size:m.position_size+m.internal_size]=rest;y[-7:-4]=[3,4,20]
        rhs=m.rhs(0,y)
        energy=lambda v:m.energy_force(*m.unpack(v)[:2])[0]+p.lorenz_capacity*np.dot(v[-7:-4],v[-7:-4])/2
        hw=1e-6/max(1,float(np.max(abs(rhs)))) # Resolve the much stiffer barrier direction.
        fd=(energy(y+hw*rhs)-energy(y-hw*rhs))/(2*hw)
        check('reciprocal-work-'+shape,abs(fd-(rhs[-4]-rhs[-3]-rhs[-2]))/max(1,abs(fd)),2e-6)
        # Dense independent probes of every material edge must exceed its lower bound.
        f=np.linspace(0,1,1001)
        points=((1-f[None,:,None])*x[m.i,None]+f[None,:,None]*x[m.j,None]).reshape(-1,2)
        gap=object_gap(points,p)[0].reshape(len(m.edges),-1)
        lower=m.clearance(x,False)
        if shape=='lobed':lower=lower.reshape(len(m.edges),-1).min(axis=1)
        check('edge-certificate-vs-dense-probes-'+shape,max(0,float(max(lower-gap.min(axis=1)))),1e-10)
    p=Parameters();m=Mesh(p);x=m.initial.copy();x[:,0]*=1.31
    d=rng.normal(size=x.shape);d/=np.linalg.norm(d);rest=m.l0
    e,g,*_=m.energy_force(x,rest);h=1e-7
    fd=(m.energy_force(x+h*d,rest)[0]-m.energy_force(x-h*d,rest)[0])/(2*h)
    check('wall-barrier-gradient',abs(fd-np.sum(g*d))/max(1,abs(fd)),2e-6)
    # Deliberately invalid inputs are refused instead of clipped into a valid image.
    for name,x in [('object',m.initial*.5),('wall',m.initial*1.5)]:
        rejected=False
        try:m.energy_force(x,rest)
        except FloatingPointError:rejected=True
        check('reject-'+name,0 if rejected else 1,.5)
    # Compress radial thickness without inversion to activate self-contact.
    x=m.initial.copy();r=np.linalg.norm(x,axis=1)
    x*=np.repeat([1.25,1.256,1.262],p.angles)[:,None]/r[:,None]
    direction=rng.normal(size=x.shape);direction/=np.linalg.norm(direction);h=1e-8
    e,g,*_=m.energy_force(x,m.l0)
    fd=(m.energy_force(x+h*direction,m.l0)[0]-m.energy_force(x-h*direction,m.l0)[0])/(2*h)
    check('active-self-contact-energy-gradient',abs(fd-np.sum(g*direction))/max(1,abs(fd)),2e-6)
    gap,normal,t=m.self_gap(x);_,d=barrier(gap,p.self_width,p.self_strength)
    force=(m.sw*d)[:,None]*normal
    net=m.scatter(m.si,force)-m.scatter(m.sa,(1-t[:,None])*force)-m.scatter(m.sb,t[:,None]*force)
    check('self-contact-action-reaction',np.linalg.norm(net.sum(axis=0)),1e-8)
    raw=np.load(ROOT/'obstacle-only/data/amplitude4-one.npz');index=np.argmin(abs(raw['time']-4.6));crossed=raw['state'][index]
    rejected=False
    try:m.energy_force(*m.unpack(crossed)[:2])
    except FloatingPointError as error:rejected='self-crossing' in str(error)
    check('reject-recorded-self-crossing',0 if rejected else 1,.5)
    _,_,certificate=m.swept_clearance(raw['state'][index-1],crossed)
    check('swept-self-certificate-rejects-crossing',0 if certificate<=0 else 1,.5)
    result=dict(checks=rows,all_passed=all(r['passed'] for r in rows))
    (ROOT/'data/validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if not result['all_passed']:raise SystemExit(1)
if __name__=='__main__':main()
