"""Assemble material/contact energy Hessians; intended for independent checking."""
import numpy as np
from contact_hug import barrier
from clay_hug import object_gap
def bd2(g,w,k):
    d=np.minimum(g,w);z=d-w
    return -k*(2*np.log(d/w)+4*z/d-z*z/(d*d))
def hessian(m,x,rest,X):
    p=m.p;k=m.position_size;H=np.zeros((k,k));I=np.eye(2)
    def add(nodes,local):
        ids=(2*nodes[:,:,None]+np.arange(2)).reshape(len(nodes),-1)
        np.add.at(H,(ids[:,:,None],ids[:,None,:]),local)
    delta=x[m.i]-x[m.j];length=np.linalg.norm(delta,axis=1);u=delta/length[:,None]
    B=m.k[:,None,None]*((1-rest/length)[:,None,None]*I+(rest/length)[:,None,None]*u[:,:,None]*u[:,None,:])
    add(np.column_stack([m.i,m.j]),np.concatenate([np.concatenate([B,-B],axis=2),np.concatenate([-B,B],axis=2)],axis=1))
    area,ag=m.areas(x);J=area/m.a0;gA=ag.reshape(-1,6)
    HA=np.zeros((6,6))
    for a,b,s in [(0,3,1),(1,2,-1),(0,5,-1),(1,4,1),(2,5,1),(3,4,-1)]:HA[a,b]=HA[b,a]=s*.5
    local=(p.bulk/(J*J*m.a0))[:,None,None]*gA[:,:,None]*gA[:,None,:]+(p.bulk*(1-1/J))[:,None,None]*HA
    _,d=barrier(J-p.area_floor,p.area_activation-p.area_floor,p.area_strength);dd=bd2(J-p.area_floor,p.area_activation-p.area_floor,p.area_strength)
    local+=(dd/m.a0)[:,None,None]*gA[:,:,None]*gA[:,None,:]+d[:,None,None]*HA
    F=m.deformation(x);S=np.sum(F*F,axis=(1,2));Q=2*J/S
    gJ=gA/m.a0[:,None];HJ=HA[None,:,:]/m.a0[:,None,None]
    gS=(2*np.einsum('tij,taj->tai',F,m.G)).reshape(-1,6)
    HS=2*np.einsum('tai,tbi,cd->tacbd',m.G,m.G,I).reshape(-1,6,6)
    gQ=2*(gJ/S[:,None]-J[:,None]*gS/S[:,None]**2)
    HQ=2*(HJ/S[:,None,None]-(gJ[:,:,None]*gS[:,None,:]+gS[:,:,None]*gJ[:,None,:])/S[:,None,None]**2-J[:,None,None]*HS/S[:,None,None]**2+2*J[:,None,None]*gS[:,:,None]*gS[:,None,:]/S[:,None,None]**3)
    _,d=barrier(Q-p.shape_floor,p.shape_activation-p.shape_floor,p.shape_strength);dd=bd2(Q-p.shape_floor,p.shape_activation-p.shape_floor,p.shape_strength)
    local+=m.a0[:,None,None]*(dd[:,None,None]*gQ[:,:,None]*gQ[:,None,:]+d[:,None,None]*HQ)
    add(m.tri,local)
    radius=np.linalg.norm(x,axis=1);unit=x/radius[:,None]
    diag=(m.mass*(p.preload+p.alpha/10*X)/radius)[:,None,None]*(I-unit[:,:,None]*unit[:,None,:])
    for sign in [-1,1]:
        dd=bd2(p.wall_half-sign*x,p.barrier_width,p.wall)
        diag+=m.mass[:,None,None]*dd[:,:,None]*I
    add(np.arange(m.n)[:,None],diag)
    i,j=m.inner_i,m.inner_j;t=m.surface_t
    pts=((1-t[None,:,None])*x[i,None]+t[None,:,None]*x[j,None]).reshape(-1,2)
    g,n=object_gap(pts,p);_,d=barrier(g,p.barrier_width,p.contact);dd=bd2(g,p.barrier_width,p.contact)
    # Only a two-dimensional normal derivative is approximated, independently
    # at quadrature points. No full-system finite-difference probes are needed.
    hn=np.empty((len(pts),2,2));eps=1e-6
    for axis in [0,1]:
        shift=np.zeros_like(pts);shift[:,axis]=eps
        hn[:,:,axis]=(object_gap(pts+shift,p)[1]-object_gap(pts-shift,p)[1])/(2*eps)
    hn=(hn+hn.transpose(0,2,1))/2
    weights=(m.surface_length[:,None]*m.surface_w).ravel()
    hc=weights[:,None,None]*(dd[:,None,None]*n[:,:,None]*n[:,None,:]+d[:,None,None]*hn)
    blend=np.tile(np.column_stack([1-t,t]),(len(i),1))
    block=np.einsum('qa,qb,qij->qaibj',blend,blend,hc).reshape(-1,4,4)
    add(np.repeat(np.column_stack([i,j]),len(t),axis=0),block)
    # Self-contact is uncommon in these runs. Isolate that force for finite
    # differences when active, so the derivative remains complete.
    g=m.self_gap(x,False)
    if np.any(g<p.self_width+1e-5):
        def grad(v):
            gap,n,t=m.self_gap(v);_,d=barrier(gap,p.self_width,p.self_strength)
            f=(m.sw*d)[:,None]*n
            return (m.scatter(m.si,f)-m.scatter(m.sa,(1-t[:,None])*f)-m.scatter(m.sb,t[:,None]*f)).ravel()
        for col in range(k):
            step=1e-7*max(1,abs(x.ravel()[col]));a=x.copy();b=x.copy();a.ravel()[col]+=step;b.ravel()[col]-=step
            H[:,col]+=(grad(a)-grad(b))/(2*step)
    return H
