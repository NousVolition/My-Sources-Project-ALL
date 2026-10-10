"""Classical/PIMD structure pilot. Thermostatted centroid mobility is diagnostic,
not a physical quantum diffusion measurement. Keep number density identical.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit
from qtip4pf import make_system,initial_positions,MASSES,PARAMS,numpy_energy

def run(out,n=216,isotope='H2O',beads=1,T=300.,seed=1,dt=.00025,eq=10.,production=20.,sample=.2,platform='OpenCL',mobility=0.):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    # Same number density for both isotopes: isolates mass and quantum effects.
    side=(n*18.01528/6.02214076e23/.997*1e21)**(1/3)
    cutoff=min(.9,.49*side)
    system=make_system(n,isotope,side,cutoff)
    (out/'system.xml').write_text(mm.XmlSerializer.serialize(system))
    p=mm.Platform.getPlatformByName(platform)
    props={'Precision':'mixed'} if platform=='OpenCL' else {'Threads':'2'} if platform=='CPU' else {}
    minimize_int=mm.VerletIntegrator(dt)
    ctx=mm.Context(system,minimize_int,p,props)
    xyz=initial_positions(n,side,seed)
    ctx.setPositions(xyz);ctx.computeVirtualSites()
    mm.LocalEnergyMinimizer.minimize(ctx,10.,1000)
    minimized=ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    del ctx,minimize_int
    it=mm.RPMDIntegrator(beads,T,1.,dt)
    it.setRandomNumberSeed(seed+1700)
    ctx=mm.Context(system,it,p,props)
    rng=np.random.default_rng(seed+3000)
    mass=np.tile(MASSES[isotope],n)
    sigma=np.sqrt(np.divide(beads*.008314462618*T,mass,out=np.zeros_like(mass),where=mass>0))
    for k in range(beads):
        it.setPositions(k,minimized)
        vel=rng.normal(size=(n*4,3))*sigma[:,None]
        vel-=np.sum(mass[:,None]*vel,axis=0)/sum(mass)
        vel[mass==0]=0
        it.setVelocities(k,vel)
    config=dict(n=n,isotope=isotope,beads=beads,T_K=T,seed=seed,dt_ps=dt,equilibration_ps=eq,production_ps=production,sample_ps=sample,mobility_NVE_ps=mobility,side_nm=side,cutoff_nm=cutoff,platform=platform,parameters=PARAMS,thermostat='PILE 1/ps; disabled for mobility branch',density_control='same number density, no barostat',contractions='none',openmm_version=mm.__version__)
    (out/'config.json').write_text(json.dumps(config,indent=2))
    tick=time.perf_counter()
    for j in range(int(eq)):
        it.step(round(1/dt))
        print(f'{out.name} equilibration {j+1}/{eq} ps, {time.perf_counter()-tick:.1f}s',flush=True)
    remaining=eq-int(eq)
    if remaining>1e-9: it.step(round(remaining/dt))
    frames=[];energies=[]
    for j in range(round(production/sample)):
        it.step(round(sample/dt))
        frame=[];ep=[]
        for k in range(beads):
            state=it.getState(k,getPositions=True,getEnergy=True)
            frame.append(np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(n,4,3)[:,:3])
            ep.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)/n)
        frames.append(frame);energies.append(ep)
        if (j+1)%10==0: print(f'{out.name} production {(j+1)*sample:.1f}/{production} ps, {time.perf_counter()-tick:.1f}s',flush=True)
    arr=np.array(frames,dtype=np.float32)
    if not np.all(np.isfinite(arr)): raise RuntimeError('Nonfinite positions')
    np.savez_compressed(out/'trajectory.npz',positions_nm=arr,potential_kJmol_per_molecule=energies,times_ps=sample*np.arange(1,len(arr)+1))
    if mobility:
        it.setApplyThermostat(False)
        mobile=[];ring_energy=[]
        for j in range(round(mobility/sample)+1):
            if j: it.step(round(sample/dt))
            pos=[np.array(it.getState(k,getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(n,4,3)[:,:3] for k in range(beads)]
            mobile.append(np.mean(pos,axis=0));ring_energy.append(it.getTotalEnergy().value_in_unit(unit.kilojoule_per_mole))
        np.savez_compressed(out/'mobility.npz',centroid_nm=mobile,times_ps=sample*np.arange(len(mobile)),ring_energy_kJmol=ring_energy)
    (out/'timing.json').write_text(json.dumps({'seconds':time.perf_counter()-tick}))
    del ctx,it

def validate(out):
    xyz=initial_positions(2,1.5,7)
    xyz=xyz.reshape(2,4,3);xyz[1]+=np.array([.32,0,0])-xyz[1,0]+xyz[0,0]
    system=make_system(2,'H2O')
    it=mm.VerletIntegrator(.0001)
    ctx=mm.Context(system,it,mm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(xyz.reshape(-1,3));ctx.computeVirtualSites()
    state=ctx.getState(getEnergy=True,getForces=True)
    energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    force=np.array(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).reshape(2,4,3)
    errors=[]
    for i in range(2):
        for a in range(3):
            for d in range(3):
                plus=xyz.copy();minus=xyz.copy();plus[i,a,d]+=1e-6;minus[i,a,d]-=1e-6
                fd=-(numpy_energy(plus)-numpy_energy(minus))/2e-6
                errors.append(abs(force[i,a,d]-fd))
    result=dict(energy_OpenMM=energy,energy_numpy=numpy_energy(xyz),energy_error=abs(energy-numpy_energy(xyz)),max_force_error=max(errors))
    assert result['energy_error']<1e-6 and result['max_force_error']<.01,result
    Path(out).write_text(json.dumps(result,indent=2))
    print(result)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);ap.add_argument('--validate',action='store_true')
    ap.add_argument('--beads',type=int,default=1);ap.add_argument('--n',type=int,default=216);ap.add_argument('--isotope',default='H2O')
    ap.add_argument('--T',type=float,default=300.);ap.add_argument('--seed',type=int,default=1);ap.add_argument('--dt',type=float,default=.00025)
    ap.add_argument('--eq',type=float,default=10.);ap.add_argument('--production',type=float,default=20.);ap.add_argument('--sample',type=float,default=.2);ap.add_argument('--platform',default='OpenCL')
    ap.add_argument('--mobility',type=float,default=0.)
    a=vars(ap.parse_args());v=a.pop('validate')
    if v: validate(a['out'])
    else: run(**a)
