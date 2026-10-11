"""Periodic four-point tensor Lagrange velocity and its exact cellwise Jacobian.

This changes tracer observation/advection accuracy, not the Eulerian NS solver.
The interpolant is continuous, but its derivative can still jump at cell faces.
It is not exactly divergence-free; convergence and volume diagnostics are required.
"""
import numpy as np

L=6.0


def sample(u,positions):
    n=u.shape[1]
    q=(positions+L/2)%L*n/L
    a=np.floor(q).astype(int);b=q-a
    w=np.stack([-b*(b-1)*(b-2)/6,(b+1)*(b-1)*(b-2)/2,
                -(b+1)*b*(b-2)/2,(b+1)*b*(b-1)/6],axis=-1)
    dw=np.stack([-(3*b*b-6*b+2)/6,(3*b*b-4*b-1)/2,
                 (-3*b*b+2*b+2)/2,(3*b*b-1)/6],axis=-1)*n/L
    velocity=np.zeros((len(q),3));gradient=np.zeros((len(q),3,3))
    for indices in np.ndindex(4,4,4):
        ix=(a+np.asarray(indices)-1)%n
        value=u[:,ix[:,0],ix[:,1],ix[:,2]].T
        weight=[w[:,j,indices[j]] for j in range(3)]
        velocity+=value*(weight[0]*weight[1]*weight[2])[:,None]
        for j in range(3):
            others=[k for k in range(3) if k!=j]
            derivative=dw[:,j,indices[j]]*weight[others[0]]*weight[others[1]]
            gradient[:,:,j]+=value*derivative[:,None]
    return velocity,gradient
