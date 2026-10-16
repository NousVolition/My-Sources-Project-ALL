"""Exact triangular solve: work counters never influence physical variables."""
import numpy as np
from scipy.sparse.linalg import splu
from scipy.sparse import diags
def factor(A):
 k=A.shape[0]-4;D=A[k:,k:].toarray();diag=np.diag(D)
 assert not A[:k,k:].nnz and np.array_equal(D,np.diag(diag)) and np.all(diag!=0)
 return splu(A[:k,:k].tocsc()),A[k:,:k].tocsr(),diag
def solve(f,b):
 lu,cross,diag=f;x=lu.solve(b[:-4]);w=(b[-4:]-cross@x)/diag;return np.r_[x,w]
def condensed_factor(A,position_size,internal_size):
 """Exact elimination of diagonal plastic rows and the small Lorenz border."""
 n=position_size;e=internal_size;end=A.shape[0]-4;q=end-n-e
 assert q in [0,3] and A[:end,end:].nnz==0
 D=A[end:,end:].toarray();wd=np.diag(D)
 assert np.array_equal(D,np.diag(wd)) and np.all(wd!=0)
 rr=A[n:n+e,n:n+e];rd=rr.diagonal()
 assert rr.nnz==np.count_nonzero(rd) and np.all(rd!=0) and A[n:n+e,n+e:end].nnz==0
 B=A[:n,n:n+e].tocsc();C=A[n:n+e,:n].tocsc()
 BC=B@diags(1/rd);S=(A[:n,:n]-BC@C).tocsc();lu=splu(S)
 if q:
  border=A[:n,n+e:end].toarray();lower=A[n+e:end,:n].toarray();lr=A[n+e:end,n:n+e].tocsc()
  lower=lower-(lr@diags(1/rd)@C).toarray();u=lu.solve(border)
  small=A[n+e:end,n+e:end].toarray()-lower@u
 else:lower=u=small=lr=None
 return n,e,lu,BC,C,rd,lower,u,small,lr,A[end:,:end].tocsr(),wd

def condensed_solve(f,b):
 n,e,lu,BC,C,rd,lower,u,small,lr,cross,wd=f;end=n+e+(0 if small is None else 3)
 right=b[:n]-BC@b[n:n+e];x=lu.solve(right)
 if small is not None:
  z=np.linalg.solve(small,b[n+e:end]-lr@(b[n:n+e]/rd)-lower@x);x-=u@z
 else:z=np.empty(0,dtype=b.dtype)
 r=(b[n:n+e]-C@x)/rd;physical=np.r_[x,r,z]
 return np.r_[physical,(b[end:]-cross@physical)/wd]

def install(solver,mesh):
 def lu(A):solver.nlu+=1;return condensed_factor(A,mesh.position_size,mesh.internal_size)
 solver.lu=lu;solver.solve_lu=condensed_solve
