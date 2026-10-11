"""Paired 3D NS pilot; Fourier operators reused from the pinned source project.

Fluid: dealiased rotational Navier-Stokes + shared solenoidal forcing, RK4.
Markers: one-way periodic cubic-spline interpolation, joint RK4 stages.
No marker forces or molecular interpretation. Length/time are model units.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy import ndimage
from vendor.numerics import Flow, L

ROOT = Path(__file__).resolve().parent


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wrap(x):
    return (x + L/2) % L - L/2


def sample(field, x, order=3):
    """Interpolate any leading batch dimensions on periodic [-3,3)^3."""
    n = field.shape[-1]
    coords = ((x + L/2) % L * n/L).T
    flat = field.reshape((-1,n,n,n))
    vals = np.stack([ndimage.map_coordinates(a, coords, order=order,
                      mode="grid-wrap", prefilter=order>1) for a in flat], axis=-1)
    return vals.reshape((len(x),) + field.shape[:-3])


def initial(f, seed, speed=1., kind="random"):
    """Same low-mode random construction as upstream history-prediction/simulate.py.

    Beltrami alternative is an ABC curl eigenfield with random phases.
    Both constructions are smooth and exactly periodic on every grid.
    """
    rng = np.random.default_rng(seed)
    axis = np.arange(f.n)*f.dx-L/2
    xyz = np.meshgrid(axis,axis,axis,indexing="ij",sparse=True)
    raw = np.zeros((3,f.n,f.n,f.n))
    if kind == "beltrami":
        x,y,z = [2*np.pi*a/L+b for a,b in zip(xyz,rng.uniform(0,2*np.pi,3))]
        raw[0] = np.sin(z)+np.cos(y)
        raw[1] = np.sin(x)+np.cos(z)
        raw[2] = np.sin(y)+np.cos(x)
    elif kind == "random":
        modes = [(1,1,0),(1,0,2),(0,2,1),(2,-1,1),(1,2,-2),(2,2,1),
                 (1,0,0),(0,1,0),(0,0,1),(1,-1,1),(3,1,0),(0,1,3)]
        for mode in modes:
            phase = sum(2*np.pi/L*m*x for m,x in zip(mode,xyz))+rng.uniform(0,2*np.pi)
            raw += rng.normal(size=3)[:,None,None,None]*np.cos(phase)/np.linalg.norm(mode)
    else:
        raise ValueError(kind)
    h = f.project(f.hat(raw))
    h[:,0,0,0] = 0
    h *= speed/np.sqrt(f.inner(h,h)/L**3)
    return h


def impulse(f, amplitude, width=.7):
    """Curl of a smooth periodic localized vector potential, fixed modes |ki|<=5.

    Potential coefficients are computed at a fixed reference grid, so each fluid
    grid receives the same Fourier polynomial, independent of sampling/aliasing.
    Amplitude specifies global RMS velocity, not peak velocity or added energy.
    """
    nref=64
    axis=np.arange(nref)*L/nref-L/2
    x,y,z=np.meshgrid(axis,axis,axis,indexing="ij",sparse=True)
    kappa=(L/(2*np.pi*width))**2
    a=np.exp(kappa*(np.cos(2*np.pi*x/L)+np.cos(2*np.pi*y/L)+np.cos(2*np.pi*z/L)-3))
    ah=np.fft.rfftn(a)/nref**3
    ph=np.zeros((f.n,f.n,f.n//2+1),complex)
    for i in range(-5,6):
        for j in range(-5,6):
            ph[i%f.n,j%f.n,:6]=ah[i%nref,j%nref,:6]*f.n**3
    h=f.project(np.stack((1j*f.k[1]*ph,-1j*f.k[0]*ph,np.zeros_like(ph))))
    return h*(amplitude/np.sqrt(f.inner(h,h)/L**3))


class PilotFlow(Flow):
    def __init__(self,n,nu,nonlinear=True,forcing=0.):
        super().__init__(n,nu,workers=1)
        # Strict two-thirds rule, valid even on grids divisible by three.
        modes=np.fft.fftfreq(n)*n
        pos=np.fft.rfftfreq(n)*n
        self.keep=(abs(modes[:,None,None]) < n/3) & (abs(modes[None,:,None]) < n/3) & (pos[None,None,:] < n/3)
        self.high=self.high & self.keep
        self.nonlinear=nonlinear
        self.force=initial(self,99001)*forcing

    def rhs_only(self,h,u=None):
        u=self.real(h) if u is None else u
        if self.nonlinear:
            w=self.real(self.curl(h))
            cross=np.cross(u,w,axisa=0,axisb=0,axisc=0)
            adv=self.project(self.hat(cross))
        else:
            adv=0.
        return adv-self.nu*self.k2*h+self.force

    def gradient(self,h):
        return np.stack([self.real(1j*k*h) for k in self.k],axis=1)


def energy(f,h):
    return .5*f.inner(h,h)/L**3


def stage(f,h,x,order):
    u=f.real(h)
    r=f.rhs_only(h,u)
    v=sample(u,x,order)
    diss=f.nu*f.inner(f.curl(h),f.curl(h))/L**3
    work=f.inner(h,f.force)/L**3
    speed=float(np.sqrt(np.sum(u*u,axis=0)).max())
    return r,v,np.array([diss,work]),speed


def advance(f,h,x,dt,order=3):
    k1,v1,b1,s1=stage(f,h,x,order)
    hmid=h+dt*k1/2
    k2,v2,b2,s2=stage(f,hmid,x+dt*v1/2,order)
    k3,v3,b3,s3=stage(f,h+dt*k2/2,x+dt*v2/2,order)
    k4,v4,b4,s4=stage(f,h+dt*k3,x+dt*v3,order)
    hn=f.project(h+dt*(k1+2*k2+2*k3+k4)/6)
    xn=x+dt*(v1+2*v2+2*v3+v4)/6
    return hn,xn,dt*(b1+2*b2+2*b3+b4)/6,max(s1,s2,s3,s4)


def perturbation_rates(f,hb,hp):
    """Exact Galerkin difference-energy budget, evaluated without time differencing.

    P=-<delta_i delta_j partial_j ub_i>; D=nu<|curl(delta)|^2>.
    Shared forcing cancels. Computing P through the nonlinear RHS is equivalent
    on the dealiased divergence-free grid and avoids nine extra FFTs per stage.
    """
    delta=hp-hb
    d=f.nu*f.inner(f.curl(delta),f.curl(delta))/L**3
    rate=f.inner(delta,f.rhs_only(hp)-f.rhs_only(hb))/L**3
    return np.array([rate+d,d])


def advance_pair(f,hb,hp,xb,xp,dt,order=3):
    hs=[hb,hp]; xs=[xb,xp]; kh=[]; kx=[]; kb=[]; pd=[]; speeds=[]
    weights=[0,.5,.5,1.]
    for j,a in enumerate(weights):
        hh=[hs[q]+a*dt*kh[-1][q] if j else hs[q] for q in range(2)]
        xx=[xs[q]+a*dt*kx[-1][q] if j else xs[q] for q in range(2)]
        stages=[stage(f,hh[q],xx[q],order) for q in range(2)]
        kh.append([s[0] for s in stages]);kx.append([s[1] for s in stages])
        kb.append(np.array([s[2] for s in stages]));speeds += [s[3] for s in stages]
        delta=hh[1]-hh[0]
        diss=f.nu*f.inner(f.curl(delta),f.curl(delta))/L**3
        rate=f.inner(delta,kh[-1][1]-kh[-1][0])/L**3
        pd.append(np.array([rate+diss,diss]))
    ww=[1,2,2,1]
    hn=[f.project(hs[q]+dt/6*sum(ww[j]*kh[j][q] for j in range(4))) for q in range(2)]
    xn=[xs[q]+dt/6*sum(ww[j]*kx[j][q] for j in range(4)) for q in range(2)]
    return *hn,*xn,dt/6*sum(ww[j]*kb[j] for j in range(4)),dt/6*sum(ww[j]*pd[j] for j in range(4)),max(speeds)


def diagnostic(f,hb,hp,t,xb,xp,wh0):
    delta=hp-hb;du=f.real(delta);dw=f.real(f.curl(delta))
    eb=energy(f,hb);ep=energy(f,hp);ed=energy(f,delta)
    weight=.5*np.sum(du*du,axis=0)
    axis=np.arange(f.n)*f.dx-L/2
    xyz=np.meshgrid(axis,axis,axis,indexing="ij",sparse=True)
    r2=sum(a*a for a in xyz)
    radius=np.sqrt(np.sum(weight*r2)/max(np.sum(weight),1e-300))
    participation=np.sum(weight)**2/max(np.sum(weight*weight)*f.n**3,1e-300)
    wb=f.curl(hb);wp=f.curl(hp)
    ub=f.real(hb);w=f.real(wb);G=f.gradient(hb)
    stretch=float(np.mean(np.einsum('iabc,jabc,ijabc->abc',w,w,G)))
    div=float(abs(f.real(1j*sum(k*c for k,c in zip(f.k,hp)))).max())
    corr=f.inner(wb,wh0)/max(np.sqrt(f.inner(wb,wb)*f.inner(wh0,wh0)),1e-300)
    paircorr=f.inner(wb,wp)/max(np.sqrt(f.inner(wb,wb)*f.inner(wp,wp)),1e-300)
    return {"t":float(t),"baseline_energy":eb,"disturbed_energy":ep,"perturbation_energy":ed,
            "radius_from_source":float(radius),"participation_fraction":float(participation),
            "outside_radius_1_5_energy_fraction":float(np.sum(weight[r2>2.25])/max(np.sum(weight),1e-300)),
            "vorticity_difference_rms":float(np.sqrt(np.mean(np.sum(dw*dw,axis=0)))),
            "velocity_difference_rms":float(np.sqrt(2*ed)),
            "marker_difference_rms":float(np.sqrt(np.mean(np.sum(wrap(xp-xb)**2,axis=1)))),
            "vorticity_pattern_correlation_initial":float(corr),"vorticity_pair_correlation":float(paircorr),
            "enstrophy_stretching_rate":stretch,"enstrophy":.5*f.inner(wb,wb)/L**3,
            "divergence_max":div,"high_band_energy_fraction":f.inner(hp*f.high,hp*f.high)/max(2*ep*L**3,1e-300)}


def run(seed,config,destination,name=None):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    name=name or f"seed{seed}"
    npz=destination/(name+".npz");meta=destination/(name+".json")
    sources={p.name:sha(p) for p in [ROOT/"core.py",ROOT/"vendor/numerics.py",ROOT/"protocol.json"]}
    if meta.exists() and npz.exists():
        old=json.loads(meta.read_text())
        if old["sources"]!=sources or old["config"]!=config or old["sha256"]!=sha(npz):
            raise ValueError(f"Stale result {name}; choose another output directory")
        return npz
    start=time.perf_counter();c=config;f=PilotFlow(c["n"],c["nu"],c["nonlinear"],c["forcing"])
    hb=initial(f,seed,c["speed"],c["initial"])
    xb=np.random.default_rng(seed+80000).uniform(-L/2,L/2,(c["markers"],3))
    dt=c["dt"];steps=round(c["end"]/dt);kick=round(c["kick_time"]/dt);stride=round(c["save_dt"]/dt)
    if not all(np.isclose(a*dt,b) for a,b in [(steps,c["end"]),(kick,c["kick_time"]),(stride,c["save_dt"])]):
        raise ValueError("Times must align with the step")
    diffusion=c["nu"]*float(f.k2[f.keep].max())*dt
    if diffusion>.5:raise ValueError("Diffusion stability gate")
    hp=hb.copy();xp=xb.copy();wh0=f.curl(hb)
    budgets=np.zeros((2,2));pb=np.zeros(2);maxcfl=0.
    times=[];bx=[];px=[];diag=[];pre=[];anchor={}; e0=energy(f,hb)
    for step in range(steps+1):
        t=step*dt
        if step==kick:
            hp=hb+impulse(f,c["amplitude"],c["width"]);xp=xb.copy()
            e_bk=energy(f,hb);e_pk=energy(f,hp);ed0=energy(f,hp-hb)
            budgets[:]=0.;pb[:]=0.
            anchor={"velocity_base":sample(f.real(hb),xb,c["interp_order"]),
                    "velocity_delta":sample(f.real(hp-hb),xb,c["interp_order"]),
                    "gradient_base":sample(f.gradient(hb),xb,c["interp_order"]),
                    "gradient_delta":sample(f.gradient(hp-hb),xb,c["interp_order"]),
                    "initial_hb":hb.copy(),"initial_hp":hp.copy()}
        if step%stride==0:
            times.append(t);bx.append(xb.copy());px.append(xp.copy())
            if step>=kick:
                row=diagnostic(f,hb,hp,t,xb,xp,wh0)
                row.update({"baseline_energy_residual":(row["baseline_energy"]-e_bk+budgets[0,0]-budgets[0,1])/max(e_bk,1e-30),
                            "disturbed_energy_residual":(row["disturbed_energy"]-e_pk+budgets[1,0]-budgets[1,1])/max(e_pk,1e-30),
                            "perturbation_budget_residual":(row["perturbation_energy"]-ed0-pb[0]+pb[1])/max(ed0,1e-30),
                            "integrated_production":float(pb[0]),"integrated_perturbation_dissipation":float(pb[1]),
                            "baseline_dissipation":float(budgets[0,0]),"baseline_forcing_work":float(budgets[0,1])})
                diag.append(row)
            else:
                pre.append({"t":t,"energy_budget_residual":float((energy(f,hb)-e0+budgets[0,0]-budgets[0,1])/e0)})
        if step==steps:break
        if step<kick:
            hb,xb,b,s=advance(f,hb,xb,dt,c["interp_order"])
            hp=hb.copy();xp=xb.copy();budgets+=b
        else:
            hb,hp,xb,xp,b,p,s=advance_pair(f,hb,hp,xb,xp,dt,c["interp_order"])
            budgets+=b;pb+=p
        maxcfl=max(maxcfl,s*dt/f.dx)
        if maxcfl>.4 or not np.isfinite(hp).all() or not np.isfinite(xp).all():
            raise FloatingPointError("CFL or finite-state gate failed")
    np.savez_compressed(npz,time=times,positions_base=bx,positions_perturbed=px,
                        final_hb=hb,final_hp=hp,**anchor)
    save_json(meta,{"name":name,"seed":seed,"config":config,"sources":sources,"sha256":sha(npz),
                    "seconds":time.perf_counter()-start,"max_cfl":maxcfl,"diffusion_number":diffusion,
                    "pre_kick_checks":pre,"diagnostics":diag})
    print(f"{name}: completed in {time.perf_counter()-start:.1f}s; E_delta ratio={diag[-1]['perturbation_energy']/max(ed0,1e-30):.4f}",flush=True)
    return npz
