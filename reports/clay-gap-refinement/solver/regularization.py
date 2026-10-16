"""Proposed reference-gradient jump energy, with a declared material length.

E = mu*ell^2/2 sum_e (L_e/d_e) ||F_left-F_right||_F^2.
The reference geometry and ell stay fixed under mesh refinement. This is a
constitutive revision, not an error estimator or an experimentally fitted law.
It is positive semidefinite and zero on every affine deformation.
"""
import numpy as np
from scipy.sparse import coo_matrix,kron,eye

def build(m):
 edges={};rows=[];cols=[];values=[];weights=[]
 for t,tri in enumerate(m.tri):
  for a,b in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]:edges.setdefault(tuple(sorted((int(a),int(b)))),[]).append(t)
 centers=m.initial[m.tri].mean(axis=1)
 for (a,b),ts in sorted(edges.items()):
  if len(ts)!=2:continue
  t,u=ts;edge=m.initial[b]-m.initial[a];length=np.linalg.norm(edge);normal=np.array([-edge[1],edge[0]])/length
  distance=abs(np.dot(centers[t]-centers[u],normal));assert distance>0
  weight=length/distance
  for j in range(2):
   row=len(weights);weights.append(weight)
   for sign,v in [(1,t),(-1,u)]:
    rows.extend([row]*3);cols.extend(m.tri[v]);values.extend(sign*m.G[v,:,j])
 C=coo_matrix((values,(rows,cols)),shape=(len(weights),m.n)).tocsc();C.sum_duplicates();C.eliminate_zeros()
 w=np.array(weights);K=(m.p.modulus*m.p.regularization_length**2)*(C.T@C.multiply(w[:,None])).tocsc();K.sum_duplicates();K.eliminate_zeros()
 m.jump_operator=C;m.jump_weights=w;m.regularizer=K;m.regularizer_hessian=kron(K,eye(2),format='csc')

def energy_gradient(m,x):
 # A reference-centered evaluation avoids cancellation for rigid translations.
 displacement=x-m.initial;g=m.regularizer@displacement
 return .5*float(np.sum(displacement*g)),g
