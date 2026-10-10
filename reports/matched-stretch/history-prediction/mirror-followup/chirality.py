"""Label-invariant signed geometry; uses only supplied past/current positions."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from features import neighbors, wrap


def signed_features(past, k=8, reverse=False):
    past=np.asarray(past)
    if past.ndim!=3 or past.shape[-1]!=3 or len(past)<3 or len(past)%2!=1:
        raise ValueError("Need an odd number >=3 of past snapshots with shape (time,markers,3)")
    if not np.isfinite(past).all():
        raise ValueError("Nonfinite positions")
    idx=neighbors(past[-1],k)
    r=wrap(past[:,idx]-past[:,:,None,:])
    d=np.linalg.norm(r,axis=-1)
    v=r/np.maximum(d[...,None],1e-12)
    relative_d=d[-1]/np.maximum(d[-1].mean(axis=1,keepdims=True),1e-12)
    moments=np.stack([(relative_d[...,None]**j*v[-1]).mean(axis=1) for j in range(3)],axis=1)
    den=np.linalg.norm(moments,axis=-1).prod(axis=-1)
    c=np.divide(np.linalg.det(moments),den,out=np.zeros_like(den),where=den>1e-12)
    if reverse:
        v=v[::-1]
    q=np.linalg.det(np.stack([v[0],v[len(v)//2],v[-1]],axis=-2))
    p=q.mean(axis=1)
    return c[:,None],np.column_stack([p,np.abs(q).mean(axis=1),p*c])
