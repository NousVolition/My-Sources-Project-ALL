"""Conservative spatial queries for the unchanged self-contact forces.

If point-to-segment distance is <= w, distance to its midpoint is <= w+L/2.
If two segments intersect, their midpoint distance is <= (L1+L2)/2.
Queries using global maximum L/2 therefore cannot omit either event. All
selected distances, gradients and proper crossing tests remain exact. The
original exhaustive conservative bounds are retained for swept certificates.
"""
import numpy as np
from scipy.spatial import cKDTree
from contact_hug import barrier

class IndexedContact:
 def __init__(self,p):
  super().__init__(p);n=len(self.self_edges)
  self.lookup=np.full((n,n),-1,dtype=int)
  vertex_rank=np.full(self.n,-1,dtype=int);vertex_rank[self.self_edges[:,0]]=np.arange(n)
  self.lookup[vertex_rank[self.si],vertex_rank[self.sa]]=np.arange(len(self.si))

 def self_contact_terms(self,x):
  if not np.isfinite(x).all():raise FloatingPointError('Nonfinite contact positions')
  edges=x[self.self_edges];mid=edges.mean(axis=1);radius=np.linalg.norm(edges[:,1]-edges[:,0],axis=1)/2;largest=float(max(radius));tree=cKDTree(mid)
  pairs=tree.query_pairs(2*largest+1e-12,output_type='ndarray')
  if len(pairs):
   ai,bi=pairs.T;labels=self.self_edges
   nonincident=(labels[ai,0]!=labels[bi,0])&(labels[ai,0]!=labels[bi,1])&(labels[ai,1]!=labels[bi,0])&(labels[ai,1]!=labels[bi,1])
   ai=ai[nonincident];bi=bi[nonincident];u,v=edges[ai,0],edges[ai,1];w,z=edges[bi,0],edges[bi,1]
   def cross(a,b):return a[:,0]*b[:,1]-a[:,1]*b[:,0]
   if np.any((cross(v-u,w-u)*cross(v-u,z-u)<0)&(cross(z-w,u-w)*cross(z-w,v-w)<0)):raise FloatingPointError('Material boundary self-crossing trial')
  groups=tree.query_ball_point(x[self.self_edges[:,0]],largest+self.p.self_width+1e-12)
  counts=np.fromiter((len(g) for g in groups),dtype=int,count=len(groups));point=np.repeat(np.arange(len(groups)),counts)
  segment=np.concatenate(groups).astype(int) if sum(counts) else np.empty(0,dtype=int)
  ids=self.lookup[point,segment];ids=np.sort(ids[ids>=0]);gradient=np.zeros_like(x)
  if not len(ids):return 0.,gradient
  a=x[self.sa[ids]];b=x[self.sb[ids]];v=x[self.si[ids]];edge=b-a
  tau=np.clip(np.sum((v-a)*edge,axis=1)/np.maximum(np.sum(edge*edge,axis=1),1e-24),0,1);diff=v-a-tau[:,None]*edge;gap=np.linalg.norm(diff,axis=1)
  keep=gap<self.p.self_width;ids=ids[keep];tau=tau[keep];diff=diff[keep];gap=gap[keep]
  if not len(ids):return 0.,gradient
  e,d=barrier(gap,self.p.self_width,self.p.self_strength);force=(self.sw[ids]*d/gap)[:,None]*diff
  gradient+=self.scatter(self.si[ids],force)-self.scatter(self.sa[ids],(1-tau[:,None])*force)-self.scatter(self.sb[ids],tau[:,None]*force)
  return float(np.dot(self.sw[ids],e)),gradient
