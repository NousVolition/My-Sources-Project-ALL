"""Assemble material/contact energy Hessians; intended for independent checking."""
import numpy as np
from scipy.sparse import coo_matrix
from contact_hug import barrier
from clay_hug import object_gap
def bd2(g,w,k):
    d=np.minimum(g,w);z=d-w
    return -k*(2*np.log(d/w)+4*z/d-z*z/(d*d))
def hessian(m,x,rest,X):
    p=m.p;k=m.position_size;rows=[];cols=[];values=[];I=np.eye(2)
    def add(nodes,local):
        ids=(2*nodes[:,:,None]+np.arange(2)).reshape(len(nodes),-1)
        rows.append(np.broadcast_to(ids[:,:,None],local.shape).ravel())
        cols.append(np.broadcast_to(ids[:,None,:],local.shape).ravel())
        values.append(local.ravel())
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
    # Differentiate each active three-vertex contact locally, in a single batch.
    # This replaces full-system probes without changing the contact energy.
    g=m.self_gap(x,False)
    active=g<p.self_width+1e-5
    if np.any(active):
        ids=np.column_stack([m.si[active],m.sa[active],m.sb[active]])
        local=x[ids];w=m.sw[active]
        def pairgrad(z):
            v,a,b=z[:,0],z[:,1],z[:,2];edge=b-a
            tau=np.clip(np.sum((v-a)*edge,axis=1)/np.sum(edge*edge,axis=1),0,1)
            diff=v-a-tau[:,None]*edge;gap=np.linalg.norm(diff,axis=1)
            _,d=barrier(gap,p.self_width,p.self_strength)
            f=(w*d/gap)[:,None]*diff
            return np.stack([f,-(1-tau[:,None])*f,-tau[:,None]*f],axis=1).reshape(-1,6)
        hs=np.empty((len(ids),6,6))
        for col in range(6):
            step=1e-7*np.maximum(1,abs(local.reshape(-1,6)[:,col]))
            aa=local.copy();bb=local.copy()
            aa.reshape(-1,6)[:,col]+=step;bb.reshape(-1,6)[:,col]-=step
            hs[:,:,col]=(pairgrad(aa)-pairgrad(bb))/(2*step[:,None])
        add(ids,(hs+hs.transpose(0,2,1))/2)
    H=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(k,k)).tocsc()
    H.eliminate_zeros()
    return H
