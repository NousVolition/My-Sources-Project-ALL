"""Identical rigid TIP3P molecules under a uniform electric-field pulse.

Run from the parent experiment directory. All inputs come from the saved water
experiment; the field couples to the existing charges, never to molecule IDs.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
import numpy as np
import openmm as mm
from openmm import unit
from scipy.stats import spearmanr

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
sys.path.insert(0,str(HERE.parent/'response-extension'))
from experiment import force_sensitivity, coordination
from impulse_response import Probe, TIMES

C=96.48533212331002  # e * volt -> kJ/mol
KB=.00831446261815324
SOURCES=np.array([0,43,86,129,172,215])


def save(path,value):
    path.write_text(json.dumps(value,indent=2,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()),encoding='utf-8')


def field_system(xml,side,interactions=True):
    system=mm.XmlSerializer.deserialize(xml)
    nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
    charges=np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                      for i in range(system.getNumParticles())]).reshape(-1,3)
    assert np.max(abs(charges-charges[0]))<1e-15 and np.max(abs(charges.sum(axis=1)))<1e-15
    if not interactions:
        for i in reversed(range(system.getNumForces())):system.removeForce(i)
    # Minimum-image O-H displacements make dipole energy invariant to periodic
    # coordinate copies. Rigid bond lengths stay far from the branch at L/2.
    force=mm.CustomCompoundBondForce(3,'-conversion*E*qH*(dz2+dz3); dz2=z2-z1-L*floor((z2-z1)/L+0.5); dz3=z3-z1-L*floor((z3-z1)/L+0.5)')
    for name,value in [('conversion',C),('E',0.),('qH',float(charges[0,1])),('L',side)]:
        force.addGlobalParameter(name,value)
    for i in range(len(charges)):force.addBond([3*i,3*i+1,3*i+2],[])
    force.setForceGroup(1);system.addForce(force)
    return system,charges


def structural(x,side):
    n=len(x);oh=x[:,1:,:]-x[:,0:1,:];oh-=side*np.rint(oh/side)
    dip=oh.sum(axis=1);direction=dip/np.linalg.norm(dip,axis=1)[:,None]
    coord,dist=coordination(x,side)
    oo=x[None,:,0,:]-x[:,None,0,:];oo-=side*np.rint(oo/side)
    bonded=np.zeros((n,n),bool)
    for k in range(2):
        cosine=np.sum(oh[:,k,None,:]*oo,axis=2)/(np.linalg.norm(oh[:,k,:],axis=1)[:,None]*np.maximum(dist,1e-10))
        bonded|=(dist<.35)&(dist>1e-8)&(cosine>np.cos(np.pi/6))
    neighbors=(dist<.35)&(dist>1e-8)
    # A descriptive connected local orientation statistic, not causality.
    mean=direction.mean(axis=0);dot=direction@direction.T
    local=float(dot[neighbors].mean()-np.dot(mean,mean)) if neighbors.any() else 0.
    return direction,coord,(bonded|bonded.T).sum(axis=1),local


def simulate(system,x,v,side,field,seed,platform_name,dt=.001,on=10.,off=10.):
    integ=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,dt*unit.picosecond)
    integ.setConstraintTolerance(1e-8);integ.setRandomNumberSeed(seed)
    platform=mm.Platform.getPlatformByName(platform_name)
    context=mm.Context(system,integ,platform,{'Precision':'double'} if platform_name=='OpenCL' else {})
    context.setPositions(x.reshape(-1,3)*unit.nanometer)
    context.setVelocities(v.reshape(-1,3)*unit.nanometer/unit.picosecond)
    context.applyConstraints(1e-8);context.applyVelocityConstraints(1e-8)
    context.setParameter('E',field)
    directions=[];coords=[];hb=[];connected=[];temperature=[];potential=[];positions=[]
    times=np.arange(round((on+off)/.1)+1)*.1
    final_on=None
    for k,t in enumerate(times):
        if k:
            # The sample at t=on includes the entire on interval. Remove the field afterward.
            if times[k-1]>=on-1e-10:context.setParameter('E',0.)
            integ.step(round(.1/dt))
        state=context.getState(getPositions=True,getEnergy=True)
        xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(-1,3,3)
        direction,coord,bond,local=structural(xyz,side)
        directions.append(direction);coords.append(coord);hb.append(bond);connected.append(local)
        kinetic=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        temperature.append(2*kinetic/(KB*6*len(x)))
        # Save only intermolecular energy, excluding the field coupling energy.
        potential.append(context.getState(getEnergy=True,groups={0}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)/len(x))
        if k in [0,round(on/.1),len(times)-1]:positions.append(xyz.copy())
        if abs(t-on)<1e-8:final_on=xyz.copy()
    result={'time_ps':times,'dipole_directions':np.array(directions,dtype=np.float32),
            'coordination_scores':np.array(coords,dtype=np.float32),'hbond_counts':np.array(hb,dtype=np.int16),
            'connected_neighbor_orientation':np.array(connected),'temperature_K':np.array(temperature),
            'water_potential_kJ_mol_per_molecule':np.array(potential),'snapshot_positions_nm':np.array(positions)}
    del context,integ
    return result,final_on


def force_control(system,x,side,charges,platform):
    integ=mm.VerletIntegrator(.001)
    context=mm.Context(system,integ,mm.Platform.getPlatformByName(platform),{'Precision':'double'} if platform=='OpenCL' else {})
    context.setPositions(x.reshape(-1,3));context.setParameter('E',.5)
    state=context.getState(getForces=True,getEnergy=True,groups={1})
    force=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).reshape(-1,3,3)
    expected=np.zeros_like(force);expected[:,:,2]=C*.5*charges
    err=float(np.max(abs(force-expected)))
    print(f'Field force check: max error={err:.12g}; first molecule={force[0].tolist()}; expected={expected[0].tolist()}',flush=True)
    assert err<1e-8  # GPU force-buffer accumulation is not exact real arithmetic.
    energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    # Translate individual coordinate copies as well as neutral whole molecules.
    shifted=x.copy();shifted[0,1,2]+=side;shifted[1,:,2]-=side
    context.setPositions(shifted.reshape(-1,3))
    e2=context.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert abs(e2-energy)<1e-9
    del context,integ
    return {'force_max_error_kJ_mol_nm':err,'periodic_image_energy_difference_kJ_mol':e2-energy,
            'net_field_force_per_molecule':float(np.max(abs(force.sum(axis=1)))),
            'all_molecules_have_same_charges':True}


def run(out,platform,replicas):
    out.mkdir(parents=True,exist_ok=True);start=time.time()
    source=HERE.parent/'response-extension'/'data'
    xml=(source/'probe_system.xml').read_text()
    side=json.loads((HERE.parent/'data'/'md_protocol.json').read_text())['box_side_nm']
    system,charges=field_system(xml,side);free,_=field_system(xml,side,False)
    first=np.load(source/'response_replica_0.npz')['positions_nm']
    checks=force_control(system,first,side,charges,platform)
    protocol={'question':'Can an external condition change collective alignment, local structure, and influence among identical model-water molecules?',
       'model':'Existing OpenMM rigid nonpolarizable TIP3P; 216 identical waters; PME; unchanged charges and pair interactions.',
       'box_side_nm':side,'density_g_cm3':.997,'temperature_K':300,'friction_ps_inverse':1.,
       'dt_ps':.001,'sample_interval_ps':.1,'field_on_ps':10.,'field_off_ps':10.,
       'field_V_nm':[0.,.1,.5,-.5],'replicas':replicas,'platform':platform,'precision':'double',
       'openmm':mm.__version__,'numpy':np.__version__,'field_energy':'-96.48533212331002*E*qH*(minimum_image_OH1_z+minimum_image_OH2_z) kJ/mol',
       'field_note':'Applied uniformly to existing partial charges. Each neutral molecule has zero net direct field force, but experiences orientation torque.',
       'statistics':'Use each replica window mean as one observation; frames and molecules are correlated. On window 6-10 ps; recovery window 16-20 ps.',
       'control':'Remove all intermolecular forces, keep each rigid molecule, the field and thermal bath; these overlapping noninteracting rotors are not liquid water.',
       'noninteracting_fields_V_nm':[0.,.5],'main_seeds':'20261020+100*replica+condition_index; noninteracting adds 50.',
       'impulse_sources':SOURCES,'impulse_times_ps':TIMES,'impulse_delta_v_nm_ps':.01,
       'impulse_note':'At each 10 ps main endpoint with E>=0, use a fresh constrained Maxwell velocity draw at 300 K (seed 20261220+replica), subtract COM velocity, then deterministic field-on RATTLE probes without thermostat. Only the six prespecified sources are sampled; no global finite-time leader claim.',
       'force_scores':'Full translational force Jacobian at E=0 and E=0.5 endpoints, epsilon=1e-4 nm; score=sum of off-diagonal 3x3 Frobenius block norms. Snapshot score, not a permanent or directional leader.',
       'source_files':{f'response_replica_{rep}.npz':hashlib.sha256((source/f'response_replica_{rep}.npz').read_bytes()).hexdigest() for rep in replicas},
       'source_system_sha256':hashlib.sha256((source/'probe_system.xml').read_bytes()).hexdigest(),
       'scope':'Short conditional model experiment, not proof of real-water switching, memory, or agency. Strong fields stress a rigid fixed-charge model; electronic polarization, bond breaking, interfaces and reactions are absent.'}
    save(out/'protocol.json',protocol);save(out/'controls.json',checks)
    field_values=[0.,.1,.5,-.5]
    for rep in replicas:
        data=np.load(source/f'response_replica_{rep}.npz');x,v=data['positions_nm'],data['velocities_nm_per_ps']
        for index,field in enumerate(field_values):
            name=f'water_r{rep}_E{field:g}';path=out/(name+'.npz')
            if path.exists():
                print('Resume: '+name,flush=True);continue
            print(f'{name}: starting 10 ps on + 10 ps off',flush=True)
            result,end=simulate(system,x,v,side,field,20261020+100*rep+index,platform)
            if field>=0:
                probe=Probe(system,platform)
                probe.context.setParameter('E',field);probe.context.setPositions(end.reshape(-1,3))
                probe.context.setVelocitiesToTemperature(300*unit.kelvin,20261220+rep)
                probe.context.applyVelocityConstraints(1e-8)
                vel=np.asarray(probe.context.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)).reshape(-1,3,3)
                vel-=np.sum(probe.mass[:,:,None]*vel,axis=(0,1))/probe.mass.sum()
                result['probe_velocities_nm_ps']=vel
                pairs=probe.sources(end,vel,SOURCES)
                result['impulse_response_norms_ps']=pairs
                if field in [0.,.5]:
                    score,jac,recip=force_sensitivity(probe.context,end,1e-4,unit)
                    result['force_scores']=score;result['force_jacobian_reciprocity_error']=np.array(recip)
                if rep==0 and field==.5:
                    half=Probe(system,platform,.0005);half.context.setParameter('E',field)
                    result['half_step_scores_ps']=half.sources(end,vel,SOURCES).sum(axis=2)
                    result['half_kick_scores_ps']=probe.sources(end,vel,SOURCES,.005).sum(axis=2)
                    perm=np.random.default_rng(20261300).permutation(len(end));inv=np.argsort(perm)
                    result['permutation']=perm
                    result['permuted_impulse_scores_ps']=probe.sources(end[perm],vel[perm],inv[SOURCES]).sum(axis=2)
                    pscore,_,_=force_sensitivity(probe.context,end[perm],1e-4,unit)
                    result['permuted_force_scores']=pscore
                    fp=Probe(free,platform);fp.context.setParameter('E',field)
                    result['free_impulse_scores_ps']=fp.sources(end,vel,SOURCES).sum(axis=2)
                    del fp.context,fp.integrator,half.context,half.integrator
                del probe.context,probe.integrator
            assert all(np.isfinite(value).all() for value in result.values())
            np.savez_compressed(path,**result)
            print(f'{name}: saved; elapsed {time.time()-start:.1f}s',flush=True)
        for index,field in enumerate([0.,.5]):
            name=f'rotors_r{rep}_E{field:g}';path=out/(name+'.npz')
            if path.exists():continue
            result,_=simulate(free,x,v,side,field,20261070+100*rep+index,platform)
            np.savez_compressed(path,**result)
            print(f'{name}: saved; elapsed {time.time()-start:.1f}s',flush=True)
    print('Driven-water measurements complete.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=HERE/'data')
    parser.add_argument('--platform',choices=['OpenCL','Reference'],default='Reference')
    parser.add_argument('--replicas',type=int,nargs='+',default=[0,1,2],choices=[0,1,2])
    args=parser.parse_args();run(args.out,args.platform,args.replicas)
