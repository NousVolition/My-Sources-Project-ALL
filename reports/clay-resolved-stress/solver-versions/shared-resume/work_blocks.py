"""Exact triangular solve: work counters never influence physical variables."""
import numpy as np
from scipy.sparse.linalg import splu
def factor(A):
 k=A.shape[0]-4;D=A[k:,k:].toarray();diag=np.diag(D)
 assert not A[:k,k:].nnz and np.array_equal(D,np.diag(diag)) and np.all(diag!=0)
 return splu(A[:k,:k].tocsc()),A[k:,:k].tocsr(),diag
def solve(f,b):
 lu,cross,diag=f;x=lu.solve(b[:-4]);w=(b[-4:]-cross@x)/diag;return np.r_[x,w]
def install(solver):
 def lu(A):solver.nlu+=1;return factor(A)
 solver.lu=lu;solver.solve_lu=solve
