"""Prospectively specified exploratory battery for identical TIP3P water.

External controllers change the environment, never molecular identity. Reduced
observable trajectories and checkpoints are saved; this is not a full atomistic
trajectory dump. Run from any working directory. Existing files require an exact
configuration hash match before resumption.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time, traceback
import numpy as np
import openmm as mm
from openmm import unit

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'driven-water-extension'))
from drive_water import structural, C, KB, save


def cases():
    rows=[dict(name='baseline',mode='zero',duration=80.),
          dict(name='steady',mode='constant',duration=40.,amplitude=.5)]
    for a in [.1,.5]:
        for f,d in [(.1,40.),(.5,20.),(1.,12.)]:
            rows.append(dict(name=f'rotate_a{a:g}_f{f:g}',mode='rotate',amplitude=a,frequency=f,duration=d))
    rows += [
        dict(name='rocking',mode='rock',amplitude=.5,frequency=.5,duration=30.),
        dict(name='competing',mode='compete',amplitude=.25,frequency=.2,frequency2=.35,duration=40.),
        dict(name='chirp',mode='chirp',amplitude=.5,duration=40.,f0=.1,f1=1.),
        dict(name='sequence_xyz',mode='sequence',amplitude=.5,duration=60.,training=30.,order=[0,1,2]*5),
        dict(name='sequence_scrambled',mode='sequence',amplitude=.5,duration=60.,training=30.,order=[2,0,1,1,2,0,0,2,1,2,1,0,1,0,2]),
        dict(name='ramp_fast',mode='ramp',amplitude=.5,duration=20.),
        dict(name='ramp_slow',mode='ramp',amplitude=.5,duration=60.),
        dict(name='history_positive',mode='history',amplitude=.5,final_field=.1,prepare=10.,duration=60.),
        dict(name='history_negative',mode='history',amplitude=-.5,final_field=.1,prepare=10.,duration=60.),
        dict(name='pulses',mode='pulse',amplitude=.5,duration=40.),
        dict(name='modulated_rotation',mode='am_rotate',amplitude=.5,frequency=.5,duration=30.),
        dict(name='spatial_rotation',mode='rotate',amplitude=.5,frequency=.5,duration=30.,spatial=1),
        dict(name='spatial_release',mode='constant',amplitude=.5,duration=60.,stop=10.,spatial=1),
        dict(name='thermal_quench',mode='zero',duration=60.,hot_until=10.,hot_temperature=350.),
        dict(name='low_drag',mode='rotate',amplitude=.5,frequency=.5,duration=20.,friction=.2),
        dict(name='high_drag',mode='rotate',amplitude=.5,frequency=.5,duration=20.,friction=5.),
        dict(name='feedback_align',mode='feedback',amplitude=.5,gain=2.,duration=40.,feedback_rotation=False),
        dict(name='feedback_turn',mode='feedback',amplitude=.5,gain=2.,duration=40.,feedback_rotation=True),
    ]
    for i,row in enumerate(rows):
        mode=row['mode'];f=row.get('frequency',0.)
        update=(.01 if mode in ['chirp','am_rotate'] or f>=1 else .02 if mode=='compete' or f>=.5 else .05 if mode in ['rotate','rock'] else .1)
        row.update(index=i,dt=.001,update=update,sample=.1,friction=row.get('friction',1.),spatial=row.get('spatial',0))
    return rows


def schedule(case,t,feedback=None):
    mode=case['mode'];a=case.get('amplitude',0.);e=np.zeros(3);phase=0.
    if t>=case.get('stop',1e30):return e,phase
    if mode in ['rotate','am_rotate']:
        phase=2*np.pi*case['frequency']*t
        if mode=='am_rotate':a*=.55+.45*np.cos(2*phase)
        e[:2]=a*np.array([np.cos(phase),np.sin(phase)])
    elif mode=='constant':e[0]=a
    elif mode=='rock':e[0]=a*np.cos(2*np.pi*case['frequency']*t)
    elif mode=='compete':
        p=2*np.pi*case['frequency']*t;q=-2*np.pi*case['frequency2']*t
        e[:2]=a*np.array([np.cos(p)+np.cos(q),np.sin(p)+np.sin(q)])
        phase=np.arctan2(e[1],e[0])
    elif mode=='chirp':
        half=case['duration']/2;f0=case['f0'];f1=case['f1'];s=(f1-f0)/half
        cycles=f0*t+s*t*t/2 if t<=half else (f0+f1)*half/2+f1*(t-half)-s*(t-half)**2/2
        phase=2*np.pi*cycles;e[:2]=a*np.array([np.cos(phase),np.sin(phase)])
    elif mode=='sequence':
        if t<case['training']:e[case['order'][min(int(t/2),len(case['order'])-1)]]=a
    elif mode=='ramp':
        # One complete loop: 0 -> +A -> -A -> 0. A finite-rate loop is not equilibrium hysteresis.
        q=t/case['duration'];e[0]=a*(4*q if q<.25 else 2-4*q if q<.75 else 4*q-4)
    elif mode=='history':e[0]=a if t<case['prepare'] else case['final_field']
    elif mode=='pulse':e[0]=a if t%8<1 else 0.
    elif mode=='feedback':
        p=np.zeros(3) if feedback is None else np.array(feedback)
        if case['feedback_rotation']:p=np.array([-p[1],p[0],0.])
        e=case['gain']*p;e*=min(1.,a/max(np.linalg.norm(e),1e-30))
    elif mode!='zero':raise ValueError(mode)
    return e,phase


def make_system(xml,side,interacting=True):
    system=mm.XmlSerializer.deserialize(xml)
    if not interacting:
        for i in reversed(range(system.getNumForces())):system.removeForce(i)
    expr='-conversion*qH*(Ex*dx+Ey*dy+Ez*dz)*profile; profile=1-spatial+spatial*(1+cos(2*pi*x1/L))/2; '
    for axis in 'xyz':
        expr+=f'd{axis}=({axis}2-{axis}1-L*floor(({axis}2-{axis}1)/L+0.5))+({axis}3-{axis}1-L*floor(({axis}3-{axis}1)/L+0.5)); '
    expr=expr.rstrip('; ')
    force=mm.CustomCompoundBondForce(3,expr)
    for name,val in [('conversion',C),('qH',.417),('Ex',0.),('Ey',0.),('Ez',0.),('spatial',0.),('L',side),('pi',np.pi)]:force.addGlobalParameter(name,val)
    for i in range(system.getNumParticles()//3):force.addBond([3*i,3*i+1,3*i+2],[])
    force.setForceGroup(1);system.addForce(force)
    return system


def context_for(system,case,seed,platform):
    integ=mm.LangevinMiddleIntegrator(300,case['friction'],case['dt'])
    integ.setRandomNumberSeed(seed);integ.setConstraintTolerance(1e-8)
    ctx=mm.Context(system,integ,mm.Platform.getPlatformByName(platform),{'Precision':'double'} if platform=='OpenCL' else {})
    return ctx,integ


def configuration(case,rep,kind,source_hash,xml_hash):
    return dict(case=case,replica=rep,kind=kind,seed=20263000+1000*rep+case['index']+100*(kind=='rotors'),input_sha256=source_hash,system_sha256=xml_hash)


def run_one(system,x,v,side,case,seed,platform):
    context,integ=context_for(system,case,seed,platform)
    context.setPositions(x.reshape(-1,3));context.setVelocities(v.reshape(-1,3))
    context.applyConstraints(1e-8);context.applyVelocityConstraints(1e-8);context.setParameter('spatial',case['spatial'])
    sample=case['sample'];update=case['update'];dt=case['dt'];n=len(x)
    assert abs(round(sample/update)*update-sample)<1e-10 and abs(round(update/dt)*dt-update)<1e-10
    times=np.arange(round(case['duration']/sample)+1)*sample
    p=[];temperature=[];potential=[];local=[];hb=[];coord_cv=[];cells=[];counts=[];nematic=[];fields=[];phases=[]
    orientations=[];orientation_times=[];snapshots=[];snapshot_times=[];feedback=np.zeros(3)
    last_field=None;last_temperature=300.
    for k,t in enumerate(times):
        if k:
            left=times[k-1]
            for sub in range(round(sample/update)):
                mid=left+(sub+.5)*update
                e,_=schedule(case,mid,feedback)
                for i,(axis,value) in enumerate(zip('xyz',e)):
                    if last_field is None or value!=last_field[i]:context.setParameter('E'+axis,float(value))
                last_field=e.copy()
                target=case.get('hot_temperature',300.) if mid<case.get('hot_until',0.) else 300.
                if target!=last_temperature:integ.setTemperature(target);last_temperature=target
                integ.step(round(update/dt))
        e,phase=schedule(case,t,feedback);fields.append(e);phases.append(phase)
        state=context.getState(getPositions=True,getEnergy=True)
        xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(n,3,3)
        u,coord,bonds,connected=structural(xyz,side)
        mean=u.mean(axis=0);feedback=mean.copy() # zero-order hold for the next 0.1 ps, defined feedback delay
        p.append(mean);local.append(connected);hb.append(bonds.mean());coord_cv.append(coord.std()/coord.mean())
        nematic.append(np.linalg.eigvalsh((3*u.T@u/n-np.eye(3))/2)[-1])
        bins=np.floor((xyz[:,0]%side)/side*2).astype(int);ids=bins[:,0]+2*bins[:,1]+4*bins[:,2]
        count=np.bincount(ids,minlength=8);means=np.zeros((8,3))
        for i in range(8):
            if count[i]:means[i]=u[ids==i].mean(axis=0)
        cells.append(means);counts.append(count)
        temperature.append(2*state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(6*n*KB))
        potential.append(context.getState(getEnergy=True,groups={0}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)/n)
        if k%5==0:orientations.append(u.astype(np.float32));orientation_times.append(t)
        if k in {0,len(times)//4,len(times)//2,3*len(times)//4,len(times)-1}:snapshots.append(xyz.copy());snapshot_times.append(t)
    last=context.getState(getPositions=True,getVelocities=True)
    result=dict(time_ps=times,polarization=np.array(p),field_V_nm=np.array(fields),reference_phase=np.array(phases),
        temperature_K=np.array(temperature),water_potential_kJ_mol=np.array(potential),connected_neighbor_orientation=np.array(local),
        hbond_neighbors=np.array(hb),coordination_CV=np.array(coord_cv),nematic_order=np.array(nematic),cell_polarization=np.array(cells),
        cell_counts=np.array(counts),dipole_directions=np.array(orientations),orientation_times_ps=np.array(orientation_times),
        snapshot_positions_nm=np.array(snapshots),snapshot_times_ps=np.array(snapshot_times),
        final_positions_nm=np.asarray(last.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(n,3,3),
        final_velocities_nm_ps=np.asarray(last.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)).reshape(n,3,3))
    assert all(np.isfinite(a).all() for a in result.values()),'Nonfinite trajectory'
    del context,integ
    return result


def controls(xml,x,v,side,platform):
    c=cases()[0];system=make_system(xml,side);ctx,integ=context_for(system,c,87263,platform)
    ctx.setPositions(x.reshape(-1,3));e=np.array([.2,-.3,.4])
    for axis,val in zip('xyz',e):ctx.setParameter('E'+axis,float(val))
    force=np.asarray(ctx.getState(getForces=True,groups={1}).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).reshape(-1,3,3)
    expected=C*np.array([-.834,.417,.417])[None,:,None]*e
    err=float(np.max(abs(force-expected)));assert err<1e-8
    ctx.setParameter('spatial',1.)
    force=np.asarray(ctx.getState(getForces=True,groups={1}).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).reshape(-1,3,3)
    differences=[]
    for atom in [0,1,2]:
        for axis in range(3):
            xp=x.copy();xm=x.copy();xp[0,atom,axis]+=1e-5;xm[0,atom,axis]-=1e-5
            ctx.setPositions(xp.reshape(-1,3));ep=ctx.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            ctx.setPositions(xm.reshape(-1,3));em=ctx.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            differences.append(abs(-(ep-em)/2e-5-force[0,atom,axis]))
    assert max(differences)<1e-6
    ctx.setPositions(x.reshape(-1,3));energy=ctx.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    shift=x.copy();shift[0,1,0]+=side;shift[1,:,:]+=np.array([side,0,0]);ctx.setPositions(shift.reshape(-1,3))
    shifted=ctx.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert abs(shifted-energy)<1e-9
    perm=np.random.default_rng(901).permutation(len(x));ctx.setPositions(x[perm].reshape(-1,3))
    pf=np.asarray(ctx.getState(getForces=True,groups={1}).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).reshape(-1,3,3)
    permutation_error=float(np.max(abs(pf-force[perm])));assert permutation_error<1e-8
    del ctx,integ
    return dict(uniform_qE_max_error=err,spatial_force_finite_difference_max_error=max(differences),
                periodic_energy_error=abs(shifted-energy),spatial_force_permutation_error=permutation_error,passed=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--platform',default='OpenCL');p.add_argument('--only',nargs='*');p.add_argument('--replicas',nargs='+',type=int,default=[0,1,2]);p.add_argument('--kinds',nargs='+',default=['water','rotors']);p.add_argument('--refinements',action='store_true');p.add_argument('--worker',action='store_true');args=p.parse_args()
    out=HERE/'data';out.mkdir(exist_ok=True,parents=True)
    src=HERE.parent/'response-extension'/'data';xml=(src/'probe_system.xml').read_text();xmlhash=hashlib.sha256(xml.encode()).hexdigest()
    side=json.loads((HERE.parent/'data'/'md_protocol.json').read_text())['box_side_nm'];planned=cases()
    protocol=dict(cases=planned,replicas=[0,1,2],kinds=['water','rotors'],molecules=216,model='Rigid fixed-charge TIP3P, unchanged identical molecule parameters',
       design='Exploratory fixed battery; no statistical discovery claim; conditions fixed before outcomes',
       total_main_trajectories=len(planned)*6,total_main_ps=sum(c['duration'] for c in planned)*6,
       field_update='Piecewise constant midpoint field, 1 fs MD steps. Case-defined update: 0.01 ps fast/chirped/modulated, 0.02 ps intermediate/competing, 0.05 ps slow rotation, 0.1 ps piecewise-constant protocols. Compare update and MD step refinement separately for the fastest rotation.',
       spatial_profile='(1+cos(2*pi*oxygen_x/L))/2; dipole coupling includes profile-gradient translational force; identical spatial rule for all molecules',
       feedback='Engineered external controller: E=clip_norm(2*P,0.5) or clip_norm(2*[-Py,Px,0],0.5), sampled and held every 0.1 ps. No intrinsic adaptive molecular interactions.',
       phase_rule='Window starts after max(one cycle, one-quarter run). Tracking screen: |Pxy|>=0.15 on >=95% frames, phase concentration>=0.9, relative phase slope<=0.05 drive rate, valid contiguous phase range<pi. Sensitivity thresholds 0.1 and 0.2. Finite-time screen, not proof of asymptotic locking.',
       history_rule='Compare late same-field windows after positive/negative preparation; compare fast/slow loop area; finite-rate loops do not establish equilibrium bistability.',
       sequence_rule='Measure teacher-off polarization and cyclic directional transitions against scrambled training and undriven windows; forcing-period coherence is not learning.',
       recurrence_rule='State discovery on baseline replica 0 only, held-out replicas 1/2; compare coarse-state persistence and ordered transitions with surrogates, no assumed number of water saddles.',
       saved='0.1 ps collective and local observables; 0.5 ps per-molecule dipole directions; five positions snapshots and final x/v. Not every atomistic integration step.',
       uncertainty='Three correlated-start configurations with independent thermostat streams; replica-level means, SD, exploratory t intervals. Broad screening increases false-positive risk; no selection-adjusted significance claimed.',
       platform=args.platform,openmm=mm.__version__,box_side_nm=side,system_sha256=xmlhash,
       limitation='Small periodic box; finite time; fixed charges/rigid molecules; no chemical reactions, electronic polarization, electrodes or superconductivity. Ghost rotor controls are not water.')
    if not args.worker:
        save(out/'protocol.json',protocol)
        first=np.load(src/'response_replica_0.npz');check=controls(xml,first['positions_nm'],first['velocities_nm_per_ps'],side,args.platform);save(out/'numerical_controls.json',check)
    if args.refinements:
        base=next(c for c in planned if c['name']=='rotate_a0.5_f1')
        planned=[]
        for label,dt,update in [('half_step',.0005,.01),('half_update',.001,.005)]:
            planned.append(dict(base,name=base['name']+'_'+label,dt=dt,update=update))
    if args.only:planned=[c for c in planned if c['name'] in args.only]
    failures=[];start=time.time()
    for kind in args.kinds:
        system=make_system(xml,side,kind=='water')
        for rep in args.replicas:
            file=src/f'response_replica_{rep}.npz';data=np.load(file);source_hash=hashlib.sha256(file.read_bytes()).hexdigest()
            for c in planned:
                name=f'{kind}_r{rep}_{c["name"]}';path=out/(name+'.npz');config=configuration(c,rep,kind,source_hash,xmlhash)
                encoded=json.dumps(config,sort_keys=True);digest=hashlib.sha256(encoded.encode()).hexdigest()
                if path.exists():
                    with np.load(path) as old:assert str(old['configuration_sha256'])==digest,'Existing result configuration mismatch'
                    print('RESUME '+name,flush=True);continue
                print('START '+name+' '+str(c['duration'])+' ps',flush=True)
                try:
                    result=run_one(system,data['positions_nm'],data['velocities_nm_per_ps'],side,c,config['seed'],args.platform)
                    result['configuration_json']=np.array(encoded);result['configuration_sha256']=np.array(digest)
                    temporary=path.with_suffix('.partial.npz');np.savez_compressed(temporary,**result);temporary.replace(path)
                    print(f'DONE {name} elapsed={time.time()-start:.1f}s',flush=True)
                except Exception as exc:
                    failures.append(dict(name=name,error=repr(exc),traceback=traceback.format_exc()))
                    save(out/f'failed_runs_{kind}_r{rep}.json',failures);print('FAILED '+name+' '+repr(exc),flush=True)
    tag='_'.join(args.kinds)+'_r'+''.join(map(str,args.replicas))+('_refinement' if args.refinements else '')
    save(out/f'completion_{tag}.json',dict(elapsed_seconds=time.time()-start,failed_runs=failures))
    if failures:raise RuntimeError(f'{len(failures)} failed runs retained in failed_runs.json')
    print('BATTERY COMPLETE',flush=True)


if __name__=='__main__':main()
