"""Read-only diagnostics; no measured quantity feeds back into the solver."""
import numpy as np

L = 6.0
NU = .01

def dot(a,b):
    return float(np.mean(np.sum(np.asarray(a)*b,axis=0)))

def norm(a):
    return np.sqrt(dot(a,a))

def mirror(u):
    n=u.shape[1]; ref=(-np.arange(n))%n
    out=u[:,ref].copy();out[0]*=-1
    return out

def transform(u,kind,inverse=False):
    n=u.shape[1]
    if kind.startswith('shift'):
        cells={'shift1':1.,'shifthalf':.5}[kind]
        if inverse: cells=-cells
        phase=np.exp(-2j*np.pi*np.fft.fftfreq(n)*cells)
        return np.fft.ifftn(np.fft.fftn(u,axes=(1,2,3))*phase[None,:,None,None],axes=(1,2,3)).real
    if kind=='rotate':
        # Proper 90 degree rotation about z: new u(x,y,z)=R u(y,-x,z).
        if inverse:
            for _ in range(3):u=transform(u,kind)
            return u
        ref=(-np.arange(n))%n
        a=u.transpose(0,2,1,3)[:,ref]
        return np.stack([-a[1],a[0],a[2]])
    return u.copy()

class Meter:
    def __init__(self,m,ops,phi,kind='none'):
        self.m=m;self.ops=ops;self.phi=phi;self.kind=kind
        self.k=np.asarray(ops[:3]);self.k2=ops[3];self.keep=ops[4];self.dx=ops[5]
        self.n=phi.shape[1];self.vol=m.L**3
        self.kidx=self.k/(2*np.pi/m.L)
        self.edge=self.keep & (np.max(np.abs(self.kidx),axis=0)>=.8*(self.n//3))
        self.shell=np.rint(np.sqrt(self.k2)/(2*np.pi/m.L)).astype(int)
    def fft(self,u):return np.fft.fftn(u,axes=(-3,-2,-1))*self.keep
    def real(self,uh):return np.fft.ifftn(uh,axes=(-3,-2,-1)).real
    def curlh(self,h):
        k=self.k
        return 1j*np.stack([k[1]*h[2]-k[2]*h[1],k[2]*h[0]-k[0]*h[2],k[0]*h[1]-k[1]*h[0]])
    def projecth(self,h):
        div=np.sum(self.k*h,axis=0);fac=np.zeros_like(self.k2)
        np.divide(1.,self.k2,out=fac,where=self.k2>0)
        return (h-self.k*div*fac)*self.keep
    def basic(self,u):
        h=self.fft(u); wh=self.curlh(h);w=self.real(wh);mag=np.sqrt(np.sum(w*w,axis=0))
        return dict(W=float(mag.max()),energy=.5*self.vol*dot(u,u),enstrophy=.5*self.vol*dot(w,w))
    def full(self,u,previous_loc=None):
        m=self.m;h=self.fft(u);w=self.real(self.curlh(h));mag=np.sqrt(np.sum(w*w,axis=0))
        loc=np.unravel_index(int(mag.argmax()),mag.shape);W=float(mag[loc]);idx=(slice(None),)+loc
        # All velocity and vorticity terms use one shared spatial location.
        gu=np.asarray([self.real(1j*self.k[j]*h) for j in range(3)]) # derivative,component,...
        gw=np.asarray([self.real(1j*self.k[j]*self.curlh(h)) for j in range(3)])
        advu=-sum(u[j]*gu[j] for j in range(3))
        advuh=self.fft(advu); advu=self.real(advuh)
        pforce=self.real(self.projecth(advuh)-advuh)
        vu=self.real(-m.NU*self.k2*h)
        advw=self.real(self.fft(-sum(u[j]*gw[j] for j in range(3))))
        stretch=self.real(self.fft(sum(w[j]*gu[j] for j in range(3))))
        vw=self.real(-m.NU*self.k2*self.curlh(h))
        omega_rhs=self.real(self.curlh(self.projecth(advuh)-m.NU*self.k2*h))
        speed=np.linalg.norm(u[idx]);wh=w[idx]/max(W,1e-300);uv=u[idx]/max(speed,1e-300)
        fields={'advection':advw,'stretching':stretch,'viscosity':vw}
        vterms={'advection':advu,'pressure':pforce,'viscosity':vu}
        vt={key:{'vector':v[idx].tolist(),'along_vorticity':float(np.dot(v[idx],wh))} for key,v in fields.items()}
        vel={key:{'vector':v[idx].tolist(),'along_velocity':float(np.dot(v[idx],uv))} for key,v in vterms.items()}
        canonical=transform(u,self.kind,True);a=canonical-mirror(canonical)
        U=norm(u);A=norm(a);D=dot(a,self.phi)
        coeff=.5*self.vol/self.n**6
        modalE=coeff*np.sum(abs(h)**2,axis=0);modalZ=self.k2*modalE
        E=float(modalE.sum());Z=float(modalZ.sum())
        div=self.real(1j*np.sum(self.k*h,axis=0))
        ph=self.projecth(h);projection=norm(self.real(ph)-u)/max(U,1e-300)
        volumes={str(c):float(np.mean(mag>=c)*self.vol) for c in (50,100,200,400)}
        half=float(np.mean(mag>=W*.5)*self.vol)
        equivdiam=2*(3*half/(4*np.pi))**(1/3)
        delta=None if previous_loc is None else np.minimum(abs(np.array(loc)-previous_loc),self.n-abs(np.array(loc)-previous_loc))*self.dx
        energy_rate=self.vol*dot(u,advu+pforce+vu)
        zrate=self.vol*dot(w,omega_rhs)
        zphysical=self.vol*dot(w,advw+stretch+vw)
        oldband=self.k2>.6*self.k2.max()
        return dict(E=A/U,D=D,U=U,A=A,B=norm(a-D*self.phi),W=W,energy=E,enstrophy=Z,
            ratio=W/max(float(mag.mean()),1e-300),peak_index=list(map(int,loc)),
            peak_position=(-m.L/2+np.array(loc)*self.dx).tolist(),
            peak_index_changed=previous_loc is not None and list(loc)!=list(previous_loc),
            peak_displacement=None if delta is None else float(np.linalg.norm(delta)),
            peak_policy='Fresh global maximum at each output; displacement is not a material trajectory.',
            fixed_threshold_volumes=volumes,half_peak_volume=half,
            half_peak_equivalent_sphere_diameter=equivdiam,
            half_peak_equivalent_sphere_diameter_cells=equivdiam/self.dx,
            vorticity_terms=vt,velocity_terms=vel,
            vorticity_product_rule_residual=norm(omega_rhs-advw-stretch-vw)/max(norm(omega_rhs),1e-300),
            vorticity_actual_rhs_along_peak=float(np.dot(omega_rhs[idx],wh)),
            energy_rate=energy_rate,viscous_energy_rate=-2*m.NU*Z,
            nonlinear_energy_rate=self.vol*dot(u,advu+pforce),
            enstrophy_rate=zrate,enstrophy_physical_terms_rate=zphysical,
            divergence_max=float(abs(div).max()),divergence_rms=float(np.sqrt(np.mean(div*div))),
            projection_relative=projection,high_band_energy_fraction=float(modalE[self.edge].sum()/max(E,1e-300)),
            high_band_enstrophy_fraction=float(modalZ[self.edge].sum()/max(Z,1e-300)),
            original_cutoff_retained_modes=int(np.sum(oldband & self.keep)),
            spectra=dict(shell_index=list(range(int(self.shell.max())+1)),
                         energy=np.bincount(self.shell.ravel(),weights=modalE.ravel()).tolist(),
                         enstrophy=np.bincount(self.shell.ravel(),weights=modalZ.ravel()).tolist(),
                         component_cutoff=self.n//3,radial_support_max=float(np.sqrt(self.k2[self.keep].max())/(2*np.pi/m.L))))

def initial(m,n,shape=0):
    x,y,z,kx,ky,kz,k2,keep,dx=m.grids(n);ops=(kx,ky,kz,k2,keep,dx)
    psi=np.exp(-12*((x-.9)**2+y*y))*(x>0)
    base=np.asarray([m.deriv(psi,ky,keep),-m.deriv(psi,kx,keep),np.zeros_like(psi)])
    base=np.asarray(m.project((base+mirror(base))*.5,kx,ky,kz,k2,keep))
    params=[(.55,.4,0.,8.,1),(.7,-.3,.2,6.,2),(.4,.5,-.3,10.,1)]
    cx,cy,cz,beta,direction=params[shape]
    blob=.25*np.exp(-beta*((x-cx)**2+(y-cy)**2+(z-cz)**2))
    inp=np.zeros_like(base);inp[direction]=blob
    c=np.asarray(m.project(inp,kx,ky,kz,k2,keep));phi=(c-mirror(c))*.5;phi/=norm(phi)
    return base,phi,ops,c
