"""Smooth periodic vorticity construction for three explicitly new experiments."""
import numpy as np
from scipy import fft, special

L = 6.0
CORE_RADIUS = 0.70
AXIAL_RADIUS = 1.0
OUTER_RADIUS = 0.70
CORE_PEAK = 80.0
DISTANCE = 1.5
EXODUS_HALF_SEPARATION = 0.9


def geometry(n):
    axis = np.arange(n) * L / n - L / 2
    x, y, z = np.meshgrid(axis, axis, axis, indexing='ij', sparse=True)
    k = 2 * np.pi * fft.fftfreq(n, d=L/n)
    kr = 2 * np.pi * fft.rfftfreq(n, d=L/n)
    ks = (k[:,None,None], k[None,:,None], kr[None,None,:])
    idx = fft.fftfreq(n) * n
    ir = fft.rfftfreq(n) * n
    # Strict 2/3 cutoff: exclude equality to avoid retained endpoint aliases.
    keep = (np.abs(idx[:,None,None]) < n/3) & (np.abs(idx[None,:,None]) < n/3) & (ir[None,None,:] < n/3)
    k2 = sum(kj*kj for kj in ks)
    return (x,y,z), ks, k2, keep


def bump(a, b, ca, cb, radius):
    k = 2*np.pi/L
    return np.exp((np.cos(k*(a-ca)) + np.cos(k*(b-cb)) - 2)/(k*radius)**2)


def initial_parts(n, case, gamma=1.0, workers=4):
    (x,y,z),ks,k2,keep = geometry(n)
    k = 2*np.pi/L
    # A localized azimuthal velocity makes a central finite tube with smooth
    # return vorticity. Written as curl(0,0,psi), it is divergence-free exactly.
    # Axial localization prevents reversing the background from merely
    # translating an axially uniform central tube half a periodic box.
    envelope=bump(x,y,0,0,CORE_RADIUS)*np.exp((np.cos(k*z)-1)/(k*AXIAL_RADIUS)**2)
    core_u=np.stack(np.broadcast_arrays(-CORE_PEAK/2*np.sin(k*y)/k*envelope,
                       CORE_PEAK/2*np.sin(k*x)/k*envelope,np.zeros((n,n,n))))
    core_h=fft.rfftn(core_u,axes=(1,2,3),workers=workers)*keep
    omega = np.zeros((3,n,n,n))
    positions=[]
    for sy in (-1,1):
        for sz in (-1,1):
            if case == 'aligned':
                cz, sign = sz*DISTANCE, -sy*sz
            elif case == 'compressive':
                cz, sign = sz*DISTANCE, sy*sz
            elif case == 'exodus':
                cz, sign = sz*EXODUS_HALF_SEPARATION, sy*sz
            else:
                raise ValueError(case)
            cy=sy*DISTANCE
            omega[0] += gamma*sign*bump(y,z,cy,cz,OUTER_RADIUS)
            positions.append((cy,cz,sign))
    bg_h=fft.rfftn(omega, axes=(1,2,3),workers=workers)*keep
    core_h[:,0,0,0]=0
    bg_h[:,0,0,0]=0
    def velocity(wh):
        kx,ky,kz=ks
        crossed=np.stack((ky*wh[2]-kz*wh[1],kz*wh[0]-kx*wh[2],kx*wh[1]-ky*wh[0]))
        inv=np.zeros_like(k2);np.divide(1,k2,out=inv,where=k2>0)
        uh=1j*crossed*inv
        uh[:,0,0,0]=0
        return uh
    return core_h,velocity(bg_h),positions


def inspect(n, case, gamma=1.0, workers=4):
    core,bg,positions=initial_parts(n,case,gamma,workers)
    _,ks,k2,keep=geometry(n)
    uh=core+bg
    real=lambda h:fft.irfftn(h,s=(n,n,n),axes=(-3,-2,-1),workers=workers)
    u=real(uh)
    center=(n//2,)*3
    jac=np.empty((3,3))
    for i in range(3):
        for j in range(3):jac[i,j]=real(1j*ks[j]*bg[i])[center]
    strain=(jac+jac.T)/2
    vals,vecs=np.linalg.eigh(strain)
    kx,ky,kz=ks
    wh=1j*np.stack((ky*uh[2]-kz*uh[1],kz*uh[0]-kx*uh[2],kx*uh[1]-ky*uh[0]))
    w=real(wh)
    wc=w[(slice(None),)+center]
    align=abs(float(wc@vecs[:,-1]))/np.linalg.norm(wc)
    outward=[]
    # Sample each initial outer centreline with periodic trilinear interpolation.
    from scipy.ndimage import map_coordinates
    sample_x=np.linspace(-L/2,L/2,32,endpoint=False)
    for cy,cz,sign in positions:
        coords=np.array([(sample_x+L/2)/(L/n),np.full(32,(cy+L/2)/(L/n)),np.full(32,(cz+L/2)/(L/n))])
        vy=map_coordinates(u[1],coords,order=1,mode='grid-wrap')
        outward.append(float(np.mean(np.sign(cy)*vy)))
    speed=float(np.sqrt(np.sum(u*u,axis=0)).max())
    return {'n':n,'case':case,'gamma':gamma,'center_spin':wc.tolist(),
            'max_spin':float(np.sqrt(np.sum(w*w,axis=0)).max()),
            'center_background_strain':strain.tolist(),'eigenvalues':vals.tolist(),
            'alignment_to_largest_eigenvector':align,
            'parallel_stretch':float(wc@strain@wc/(wc@wc)),
            'outer_mean_outward_speeds':outward,'max_speed':speed,
            'mean_velocity':[float(c.mean()) for c in u],
            'core_radius_cells':CORE_RADIUS/(L/n),'outer_radius_cells':OUTER_RADIUS/(L/n)}


if __name__=='__main__':
    import json,time
    started=time.perf_counter()
    unit=inspect(80,'aligned')
    gamma=4.0/unit['parallel_stretch']
    print(json.dumps({'unit':unit,'gamma_for_stretch_4':gamma},indent=2),flush=True)
    for case in ('aligned','compressive','exodus'):
        print(json.dumps(inspect(48,case,gamma),indent=2),flush=True)
    print('elapsed',time.perf_counter()-started)
