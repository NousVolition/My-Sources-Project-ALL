"""Dilute heavy-particle observations of the air runs; no feedback or phase change.

dx/dt=v, dv/dt=(u-v)/tau+g. Spherical Stokes drag, heavy-particle
approximation; no Basset memory, added mass, Brownian diffusion or collisions.
The gas kick leaves particle velocities continuous. All species start from an
identical copied particle state at the gas kick (no pre-kick particle history).
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import ndimage
from core import ROOT,PilotFlow,advance_pair,sample,wrap,save_json,L,sha


def coefficients(u):
    return np.stack([ndimage.spline_filter(a,order=3,mode='grid-wrap') for a in u])


def interp_coeff(c,x):
    coords=((x+L/2)%L*c.shape[-1]/L).T
    return np.stack([ndimage.map_coordinates(a,coords,order=3,mode='grid-wrap',prefilter=False) for a in c],axis=-1)


def particle_step(x,v,tau,g,c0,c1,dt,max_ratio=.25):
    count=max(1,int(np.ceil(dt/(max_ratio*np.min(tau)))));h=dt/count
    def rhs(xx,vv,a):
        u=interp_coeff(c0,xx)*(1-a)+interp_coeff(c1,xx)*a
        return vv,(u-vv)/tau[:,None]+g
    for j in range(count):
        a=j/count
        x1,v1=rhs(x,v,a)
        x2,v2=rhs(x+h*x1/2,v+h*v1/2,a+.5/count)
        x3,v3=rhs(x+h*x2/2,v+h*v2/2,a+.5/count)
        x4,v4=rhs(x+h*x3,v+h*v3,a+1/count)
        x=x+h/6*(x1+2*x2+2*x3+x4)
        v=v+h/6*(v1+2*v2+2*v3+v4)
    return x,v


def run_particles(path,output,suffix='',ratio=.25):
    z=np.load(path);m=json.loads(path.with_suffix('.json').read_text());c=m['config'];seed=m['seed'];M=c['markers']
    f=PilotFlow(c['n'],c['nu'],c['nonlinear'],c['forcing'])
    hb=z['initial_hb'];hp=z['initial_hp'];ix=np.argmin(abs(z['time']-c['kick_time']));x0=z['positions_base'][ix]
    vel=sample(f.real(hb),x0);xb=x0.copy();xp=x0.copy()
    species=['fog_gravity','fog_zero_g','smoke_gravity','smoke_zero_g']
    tau=np.repeat([.02,.02,.002,.002],M);rho=np.repeat([1000.,1000.,1800.,1800.],M)
    g=np.zeros((4*M,3));g[:M,2]=-.981;g[2*M:3*M,2]=-.981
    diameter=np.sqrt(18*1.8e-5*tau*.01/rho)
    xx=[np.tile(x0,(4,1)),np.tile(x0,(4,1))];vv=[np.tile(vel,(4,1)),np.tile(vel,(4,1))]
    rows=[];tracks=[];maxre=np.zeros(4);maxslip=np.zeros(4)
    dt=c['dt'];nt=round((c['end']-c['kick_time'])/dt);stride=round(c['save_dt']/dt)
    old=[coefficients(f.real(hb)),coefficients(f.real(hp))]
    for step in range(nt+1):
        if step%stride==0:
            for j,name in enumerate(species):
                sl=slice(j*M,(j+1)*M)
                row={'seed':seed,'species':name,'t':float(step*dt),'tau':float(tau[j*M]),'diameter_m':float(diameter[j*M]),
                     'paired_displacement_rms':float(np.sqrt(np.mean(np.sum(wrap(xx[1][sl]-xx[0][sl])**2,axis=1)))),
                     'distance_from_passive_rms':float(np.sqrt(np.mean(np.sum(wrap(xx[1][sl]-xp)**2,axis=1))))}
                rows.append(row)
            tracks.append(np.stack(xx))
        for q in range(2):
            slip=np.linalg.norm(interp_coeff(old[q],xx[q])-vv[q],axis=1)
            for j in range(4):
                sl=slice(j*M,(j+1)*M)
                maxslip[j]=max(maxslip[j],float(slip[sl].max()))
                maxre[j]=max(maxre[j],float(np.max(slip[sl]*.1*diameter[sl]/1.5e-5)))
        if step==nt:break
        hb,hp,xb,xp,*_=advance_pair(f,hb,hp,xb,xp,dt,c['interp_order'])
        new=[coefficients(f.real(hb)),coefficients(f.real(hp))]
        for q in range(2):xx[q],vv[q]=particle_step(xx[q],vv[q],tau,g,old[q],new[q],dt,ratio)
        old=new
    name=f'air_particles_{seed}{suffix}';output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/(name+'.npz'),positions=tracks,final_velocities=vv)
    save_json(output/(name+'.json'),{'seed':seed,'source':str(path.name),'source_sha256':sha(path),
              'source_code_sha256':sha(__file__),'config':c,'max_step_over_tau':ratio,'rows':rows,
              'max_particle_Re_by_species':maxre.tolist(),'max_slip_by_species':maxslip.tolist(),
              'gas_reproduction_max_error':float(np.max(abs(hp-z['final_hp']))),
              'assumptions':'L0=1 mm, U0=0.1 m/s, t0=0.01 s, nu_air=1.5e-5 m2/s, mu=1.8e-5 Pa s; heavy noninteracting spherical particles; isothermal droplets in saturated air; no phase change, collisions or feedback; smoke represented by idealized sphere response time, not real aggregate chemistry.'})
    print(name+' complete; max particle Re '+str(maxre.max()),flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');a=ap.parse_args()
    p=json.loads((ROOT/'protocol.json').read_text())
    for seed in p['regime_seeds']:run_particles(a.data/'regimes'/f'air_{seed}.npz',a.data/'particles')
    run_particles(a.data/'regimes'/'air_7000.npz',a.data/'particles','_half_particle_step',.125)
    run_particles(a.data/'refinement'/'air_dt_7000.npz',a.data/'particles','_half_fluid_step')
    run_particles(a.data/'refinement'/'air_grid_7000.npz',a.data/'particles','_fine_grid')


if __name__=='__main__':main()
