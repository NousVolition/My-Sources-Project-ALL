"""Conservative broad phase; omitted distance work has exactly zero barrier force.

Point-to-box distance is a lower bound on point-to-segment distance. Far pairs
return that lower bound for conservative path certificates. All potentially
active pairs retain the original exact distance and gradient. Segment crossings
are tested exactly whenever their bounding boxes overlap.
"""
import numpy as np

class FastSelfContact:
    def self_gap(self, x, gradient=True):
        a=x[self.sa];b=x[self.sb];v=x[self.si]
        delta=np.maximum(np.maximum(np.minimum(a,b)-v,v-np.maximum(a,b)),0)
        gap=np.linalg.norm(delta,axis=1)
        near=gap < self.p.self_width + .001
        edge=b[near]-a[near];rel=v[near]-a[near]
        tau=np.clip(np.sum(rel*edge,axis=1)/np.maximum(np.sum(edge*edge,axis=1),1e-24),0,1)
        diff=rel-tau[:,None]*edge;distance=np.linalg.norm(diff,axis=1)
        gap[near]=distance
        if not gradient:return gap
        edges=x[self.self_edges]
        low=edges.min(axis=1);high=edges.max(axis=1)
        overlap=np.all(np.maximum(low[self.pair_a],low[self.pair_b]) <= np.minimum(high[self.pair_a],high[self.pair_b]),axis=1)
        pa=self.pair_a[overlap];pb=self.pair_b[overlap]
        u,v1=edges[pa,0],edges[pa,1];w,z=edges[pb,0],edges[pb,1]
        def cross(a,b):return a[:,0]*b[:,1]-a[:,1]*b[:,0]
        if np.any((cross(v1-u,w-u)*cross(v1-u,z-u)<0)&(cross(z-w,u-w)*cross(z-w,v1-w)<0)):
            raise FloatingPointError('Material boundary self-crossing trial')
        normal=np.zeros((len(gap),2));normal[near]=diff/np.maximum(distance[:,None],1e-24)
        t=np.zeros(len(gap));t[near]=tau
        return gap,normal,t
