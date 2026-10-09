"""Reproducible positional-influence experiment; see report.html for interpretation.

python -m pip install -r requirements.txt
python experiment.py --mode all --platform CPU --out data
python experiment.py --mode analyze --out data  # reanalyze saved trajectories

The recorded run uses OpenCL mixed precision; trajectories need not be bitwise
identical across platforms. All molecule labels are zero-based.
"""
from pathlib import Path
import argparse, csv, json, time, sys, platform as pyplatform
import numpy as np

SEED = 20261008


def save_json(path, obj):
    def native(v):
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, np.generic): return v.item()
        raise TypeError(type(v).__name__)
    path.write_text(json.dumps(obj, indent=2, default=native), encoding='utf-8')


def csv_write(path, columns, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f); writer.writerow(columns); writer.writerows(rows)


def weights(points, ell=.18, periodic=True):
    delta = points[:, None, :] - points[None, :, :]
    if periodic: delta -= np.rint(delta)
    w = np.exp(-np.sum(delta**2, axis=2)/(2*ell**2))
    np.fill_diagonal(w, 0)
    return w/w.sum(axis=1).mean()


def response(w, t=.5):
    lap = np.diag(w.sum(axis=1))-w
    lam, vec = np.linalg.eigh(lap)
    heat = (vec * np.exp(-t*np.maximum(lam, 0))) @ vec.T
    # Off-source share of an equal unit impulse. Total sum remains one.
    return heat.sum(axis=0)-np.diag(heat), heat


def stats(x):
    return dict(mean=np.mean(x), sd=np.std(x), cv=np.std(x)/np.mean(x)
                if np.mean(x) else 0., minimum=np.min(x), maximum=np.max(x),
                max_over_min=np.max(x)/np.min(x) if np.min(x)>0 else None,
                leader=int(np.argmax(x)))


def run_toy(out):
    print('Toy experiment', flush=True)
    rng = np.random.default_rng(SEED)
    pos = rng.random((64, 3))
    w = weights(pos); score, heat = response(w)
    leader, low = int(np.argmax(score)), int(np.argmin(score))
    np.savez_compressed(out/'toy_arrays.npz', positions=pos, weights=w,
                        response=heat, scores=score)
    csv_write(out/'toy_molecules.csv', ['id','x','y','z','weighted_degree','pulse_share'],
              [[i,*pos[i],w[i].sum(),score[i]] for i in range(64)])
    errors, changed = [], 0
    for _ in range(1000):
        perm = rng.permutation(64) # new label i occupies old site perm[i]
        new, _ = response(weights(pos[perm]))
        errors.append(np.max(np.abs(new-score[perm])))
        changed += int(np.argmax(new) != leader)
    swap = np.arange(64); swap[leader],swap[low]=low,leader
    swapped, _ = response(weights(pos[swap]))
    lattice = np.array(np.meshgrid(*([np.arange(4)/4]*3),indexing='ij')).reshape(3,-1).T
    uniform = (np.ones((64,64))-np.eye(64))/63
    controls = {}
    for name, matrix in [('periodic_cubic_lattice',weights(lattice)),
                         ('equal_all_to_all',uniform), ('no_interactions',np.zeros((64,64)))]:
        s, _ = response(matrix)
        controls[name] = dict(**stats(s), spread=np.ptp(s),
                              tied_leaders=int(np.isclose(s,s.max(),rtol=0,atol=1e-12).sum()))
    # Fresh independently sampled configurations, all retained without selection.
    ensemble = []
    for repetition in range(200):
        points = rng.random((64,3))
        for periodic in [True,False]:
            for ell in [.12,.18,.30]:
                s, _ = response(weights(points,ell,periodic))
                ensemble.append([repetition,periodic,ell,np.std(s)/np.mean(s),s.max()/s.min(),np.argmax(s)])
    csv_write(out/'toy_ensemble.csv',['repetition','periodic','ell','cv','max_over_min','leader'],ensemble)
    ensembles = {}
    for periodic in [True,False]:
        for ell in [.12,.18,.30]:
            vals=np.array([row[3:5] for row in ensemble if row[1]==periodic and row[2]==ell],float)
            ensembles[f'{"periodic" if periodic else "open"}_ell_{ell}'] = {
                'n':len(vals),'cv_mean':vals[:,0].mean(),'cv_sd':vals[:,0].std(ddof=1),
                'cv_mean_se':vals[:,0].std(ddof=1)/np.sqrt(len(vals)),
                'cv_2.5_50_97.5_percentiles':np.percentile(vals[:,0],[2.5,50,97.5]),
                'ratio_2.5_50_97.5_percentiles':np.percentile(vals[:,1],[2.5,50,97.5])}
    time_scores={str(t):stats(response(w,t)[0]) for t in [.05,.5,5.,50.]}
    checks={
        'symmetric_response_max_error':np.max(np.abs(heat-heat.T)),
        'conserved_impulse_max_error':np.max(np.abs(heat.sum(axis=0)-1)),
        'permutation_max_absolute_error':max(errors),
        'permutation_count':1000,'changed_leader_count':changed,
        'expected_change_fraction':63/64,
        'swap_old_leader':leader,'swap_old_lowest':low,'swap_new_leader':int(np.argmax(swapped)),
        'uniform_long_time_share':63/64,
        'long_time_max_error':np.max(np.abs(response(w,50.)[0]-63/64))}
    assert checks['permutation_max_absolute_error']<1e-12
    assert checks['conserved_impulse_max_error']<1e-12
    assert int(np.argmax(swapped)) == low
    assert all(v['spread']<1e-12 for v in controls.values())
    result={'seed':SEED,'N':64,'ell':.18,'t':.5,'primary':stats(score),
            'checks':checks,'controls':controls,'time_sensitivity':time_scores,'ensembles':ensembles}
    save_json(out/'toy_results.json',result)
    print(json.dumps(result['primary'],default=lambda x:float(x)),flush=True)


def force_sensitivity(context, xyz, eps, unit):
    """Full intermolecular translational force Jacobian; no integration here.

    J[j,a,i,b] = d(net force on molecule j, component a)/d(R_i,b).
    Move all three atoms of i together, preserving its orientation and geometry.
    """
    n = len(xyz)
    jac = np.empty((n,3,n,3))
    for i in range(n):
        for b in range(3):
            fs=[]
            for sign in [1.,-1.]:
                moved=xyz.copy(); moved[i,:,b] += sign*eps
                context.setPositions(moved.reshape(-1,3)*unit.nanometer)
                f=context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
                fs.append(np.asarray(f).reshape(n,3,3).sum(axis=1))
            jac[:,:,i,b]=(fs[0]-fs[1])/(2*eps)
    context.setPositions(xyz.reshape(-1,3)*unit.nanometer)
    pair=np.sqrt(np.sum(jac**2,axis=(1,3)))
    np.fill_diagonal(pair,0)
    score=pair.sum(axis=0)
    reciprocity=np.linalg.norm(jac-jac.transpose(2,3,0,1))/np.linalg.norm(jac)
    return score,jac,reciprocity


def run_md(out, platform_name):
    import openmm as mm
    from openmm import app,unit
    ff=app.ForceField('tip3p.xml')
    model=app.Modeller(app.Topology(),[])
    model.addSolvent(ff,numAdded=216)
    xyz=np.array(model.positions.value_in_unit(unit.nanometer)).reshape(216,3,3)
    old_box=np.array(model.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
    # 18.015324 g/mol is the molecular mass in the installed TIP3P XML.
    molar_mass=18.015324
    density=.997
    side=(216*molar_mass/6.02214076e23/density*1e21)**(1/3)
    oxygen=xyz[:,0,:].copy()
    xyz=xyz-oxygen[:,None,:]+oxygen[:,None,:]*(side/old_box[0,0])
    model.positions=xyz.reshape(-1,3)*unit.nanometer
    model.topology.setPeriodicBoxVectors(tuple(mm.Vec3(*(np.eye(3)[i]*side)) for i in range(3))*unit.nanometer)
    system=ff.createSystem(model.topology,nonbondedMethod=app.PME,
        nonbondedCutoff=.8*unit.nanometer,constraints=app.HBonds,rigidWater=True,
        ewaldErrorTolerance=1e-5)
    actual_mass=sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(648))
    actual_density=actual_mass/6.02214076e23/(side**3*1e-21)
    (out/'water_system.xml').write_text(mm.XmlSerializer.serialize(system),encoding='utf-8')
    with (out/'initial_water.pdb').open('w') as f: app.PDBFile.writeFile(model.topology,model.positions,f)
    platform=mm.Platform.getPlatformByName(platform_name)
    properties=({'Threads':'2','DeterministicForces':'true'} if platform_name=='CPU'
                else {'Precision':'mixed'} if platform_name=='OpenCL' else {})
    # Use a separate force-evaluation context so perturbations cannot alter trajectories.
    force_integrator=mm.VerletIntegrator(.001*unit.picoseconds)
    force_platform=platform if platform_name=='OpenCL' else mm.Platform.getPlatformByName('Reference')
    force_properties={'Precision':'double'} if platform_name=='OpenCL' else {}
    force_context=mm.Context(system,force_integrator,force_platform,force_properties)
    provenance={'seed_base':SEED,'replicas':3,'molecules':216,'box_side_nm':side,
        'density_g_cm3':actual_density,'model':'OpenMM tip3p.xml, rigid, nonpolarizable',
        'ensemble':'NVT','temperature_K':300,'time_step_fs':2,'friction_per_ps':1,
        'equilibration_ps':20,'production_ps':50,'sample_interval_ps':.1,
        'force_score_interval_ps':5,'force_displacement_nm':1e-4,
        'cutoff_nm':.8,'electrostatics':'PME','ewald_error_tolerance':1e-5,
        'platform':platform_name,'properties':properties,
        'force_platform':force_platform.getName(),'force_properties':force_properties,'openmm':mm.__version__,
        'numpy':np.__version__,'python':sys.version,'os':pyplatform.platform(),
        'replica_independence':'Same initial oxygen sites; independently randomized orientations, thermostat and velocity seeds.'}
    save_json(out/'md_protocol.json',provenance)
    rows=[]
    for rep in range(3):
        start=time.time(); seed=SEED+rep
        integrator=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,.002*unit.picoseconds)
        integrator.setRandomNumberSeed(seed)
        integrator.setConstraintTolerance(1e-6)
        sim=app.Simulation(model.topology,system,integrator,platform,properties)
        # Independently orient each intact water while preserving its oxygen site.
        from scipy.spatial.transform import Rotation
        rotations=Rotation.random(216,random_state=np.random.default_rng(seed)).as_matrix()
        local=xyz-xyz[:,0:1,:]
        initial=np.einsum('nij,naj->nai',rotations,local)+xyz[:,0:1,:]
        sim.context.setPositions(initial.reshape(-1,3)*unit.nanometer)
        sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=2000)
        sim.context.setVelocitiesToTemperature(300*unit.kelvin,seed+100)
        print(f'MD replica {rep}: equilibrating 20 ps',flush=True)
        sim.step(10000)
        (out/f'equilibrated_state_{rep}.xml').write_text(mm.XmlSerializer.serialize(sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True)),encoding='utf-8')
        frames=[]; energies=[]; temperatures=[]; fscores=[]; recips=[]
        for frame in range(500):
            sim.step(50)
            state=sim.context.getState(getPositions=True,getEnergy=True)
            frame_xyz=np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(216,3,3)
            frames.append(frame_xyz)
            energies.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            kinetic=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            temperatures.append(2*kinetic/(.00831446261815324*(6*216-3)))
            if (frame+1)%50==0:
                s,jac,recip=force_sensitivity(force_context,frame_xyz,1e-4,unit)
                fscores.append(s);recips.append(recip)
                rows.append([rep,(frame+1)/10,int(np.argmax(s)),s.min(),s.max(),s.mean(),s.std()/s.mean(),recip])
                print(f'MD replica {rep}: production {(frame+1)/10:g} ps; force leader {np.argmax(s)}; elapsed {time.time()-start:.1f}s',flush=True)
        np.savez_compressed(out/f'md_replica_{rep}.npz',positions_nm=frames,
             potential_kJ_mol=energies,temperature_K=temperatures,force_scores=fscores,
             jacobian_reciprocity_errors=recips,box_side_nm=side,seed=seed)
        if rep==0:
            # Independent recomputation at half the perturbation size.
            refined,_,_=force_sensitivity(force_context,frames[-1],5e-5,unit)
            perm=np.random.default_rng(SEED+1000).permutation(216)
            pscore,pjac,_=force_sensitivity(force_context,frames[-1][perm],1e-4,unit)
            base=np.array(fscores[-1])
            force_context.setPositions(frames[-1].reshape(-1,3)*unit.nanometer)
            energy_before=force_context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            force_context.setPositions(frames[-1][perm].reshape(-1,3)*unit.nanometer)
            energy_after=force_context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            from scipy.stats import spearmanr
            control={'base_leader':int(np.argmax(base)),
                'half_step_leader':int(np.argmax(refined)),
                'half_step_relative_L2_difference':np.linalg.norm(refined-base)/np.linalg.norm(base),
                'half_step_max_relative_difference':np.max(np.abs(refined-base)/base),
                'half_step_spearman':spearmanr(refined,base).statistic,
                'permutation_new_leader':int(np.argmax(pscore)),
                'permutation_expected_new_leader':int(np.flatnonzero(perm==np.argmax(base))[0]),
                'permutation_relative_L2_difference':np.linalg.norm(pscore-base[perm])/np.linalg.norm(base),
                'permutation_spearman':spearmanr(pscore,base[perm]).statistic,
                'permutation_energy_difference_kJ_mol':energy_after-energy_before}
            save_json(out/'md_permutation_control.json',control)
            np.savez_compressed(out/'md_control_arrays.npz',base_scores=base,
                refined_scores=refined,permuted_scores=pscore,permutation=perm,
                base_jacobian=jac,permuted_jacobian=pjac)
            assert control['base_leader']==control['half_step_leader']
            assert control['permutation_new_leader']==control['permutation_expected_new_leader']
            assert control['half_step_relative_L2_difference']<.01
            assert control['permutation_relative_L2_difference']<.01
        del sim,integrator
    csv_write(out/'md_force_snapshots.csv',['replica','production_ps','leader','min','max','mean','cv','reciprocity_error'],rows)
    del force_context,force_integrator


def coordination(xyz,side):
    delta=xyz[:,None,0,:]-xyz[None,:,0,:]
    delta-=side*np.rint(delta/side)
    dist=np.linalg.norm(delta,axis=2)
    # Continuous oxygen-neighborhood proxy, not a force or hydrogen-bond score.
    w=np.exp(-.5*((dist-.28)/.05)**2)
    np.fill_diagonal(w,0)
    return w.sum(axis=1),dist


def analyze_md(out):
    from scipy.stats import spearmanr,rankdata
    results=[]; leader_rows=[]; rdf_counts=np.zeros(90); edges=np.linspace(0,.9,91); total_frames=0
    for rep in range(3):
        data=np.load(out/f'md_replica_{rep}.npz')
        xyz=data['positions_nm']; side=float(data['box_side_nm']); n=xyz.shape[1]
        coord=[]; hb=[]
        for k,frame in enumerate(xyz):
            score,dist=coordination(frame,side); coord.append(score)
            oo=frame[None,:,0,:]-frame[:,None,0,:]
            oo-=side*np.rint(oo/side)
            bonded=np.zeros((n,n),bool)
            # donor O->H within 30 degrees of donor O->acceptor O; OO < .35 nm
            for h in [1,2]:
                oh=frame[:,h,:]-frame[:,0,:]
                oh-=side*np.rint(oh/side)
                cosine=np.sum(oh[:,None,:]*oo,axis=2)/(np.linalg.norm(oh,axis=1)[:,None]*np.maximum(dist,1e-10))
                bonded |= (dist<.35)&(dist>1e-8)&(cosine>np.cos(np.pi/6))
            undirected=bonded|bonded.T
            hb.append(undirected.sum(axis=1))
            rdf_counts+=np.histogram(dist[np.triu_indices(n,1)],bins=edges)[0]
            total_frames+=1
        coord=np.array(coord); hb=np.array(hb); leaders=coord.argmax(axis=1)
        force=data['force_scores']; fleaders=force.argmax(axis=1)
        ranks=rankdata(coord,axis=1)
        # Rank correlation follows the same molecular labels over time.
        lags={}
        for lag in [1,10,50,100]:
            a=ranks[:-lag]-ranks[:-lag].mean(axis=1,keepdims=True)
            b=ranks[lag:]-ranks[lag:].mean(axis=1,keepdims=True)
            lags[str(lag*.1)]=np.mean(np.sum(a*b,axis=1)/np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1)))
        split=spearmanr(coord[:250].mean(axis=0),coord[250:].mean(axis=0)).statistic
        # Account for conservation of sample count; do not treat frames as independent.
        result={'replica':rep,'temperature_mean_K':data['temperature_K'].mean(),
            'temperature_sd_K':data['temperature_K'].std(),
            'potential_per_molecule_mean_kJ_mol':data['potential_kJ_mol'].mean()/n,
            'potential_first_half_per_molecule':data['potential_kJ_mol'][:250].mean()/n,
            'potential_second_half_per_molecule':data['potential_kJ_mol'][250:].mean()/n,
            'coordination_cv_frame_mean':np.mean(coord.std(axis=1)/coord.mean(axis=1)),
            'coordination_unique_leaders':len(np.unique(leaders)),
            'coordination_leader_change_fraction_0.1ps':np.mean(leaders[1:]!=leaders[:-1]),
            'coordination_max_leader_frame_share':np.bincount(leaders,minlength=n).max()/len(leaders),
            'coordination_rank_correlation_by_lag_ps':lags,
            'coordination_first_vs_second_half_mean_spearman':split,
            'hydrogen_bond_unique_neighbor_mean':hb.mean(),
            'force_cv_snapshot_mean':np.mean(force.std(axis=1)/force.mean(axis=1)),
            'force_max_min_ratio_snapshot_median':np.median(force.max(axis=1)/force.min(axis=1)),
            'force_unique_leaders':len(np.unique(fleaders)),
            'force_leader_change_fraction_5ps':np.mean(fleaders[1:]!=fleaders[:-1]),
            'force_leaders':fleaders,
            'force_neighbor_score_spearman_mean':np.mean([spearmanr(force[i],coord[(i+1)*50-1]).statistic for i in range(10)]),
            'force_jacobian_reciprocity_error_max':data['jacobian_reciprocity_errors'].max()}
        results.append(result)
        csv_write(out/f'md_molecule_time_means_{rep}.csv',['id','mean_coordination','first_half','second_half','leader_frame_count'],
            [[i,coord[:,i].mean(),coord[:250,i].mean(),coord[250:,i].mean(),np.sum(leaders==i)] for i in range(n)])
        for frame in range(500): leader_rows.append([rep,(frame+1)*.1,leaders[frame],coord[frame].min(),coord[frame].max(),hb[frame].mean()])
        np.savez_compressed(out/f'md_analysis_{rep}.npz',coordination_scores=coord,
                             hydrogen_bond_neighbor_counts=hb,force_scores=force)
    shell=4*np.pi/3*(edges[1:]**3-edges[:-1]**3)
    rdf=rdf_counts/(total_frames*n*(n-1)/2*shell/side**3)
    csv_write(out/'oxygen_rdf.csv',['r_nm','g_OO'],zip((edges[1:]+edges[:-1])/2,rdf))
    csv_write(out/'md_structural_frames.csv',['replica','production_ps','leader','min_score','max_score','mean_hbond_neighbors'],leader_rows)
    save_json(out/'md_results.json',{'replicas':results,
        'rdf_first_peak_nm':float(((edges[1:]+edges[:-1])/2)[np.argmax(rdf)]),
        'rdf_first_peak_height':rdf.max(),
        'uncertainty_note':'Frames and molecules are correlated; no iid-frame p-values or confidence intervals are claimed.'})
    print('MD analysis complete',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=['all','toy','md','analyze'],default='all')
    p.add_argument('--out',type=Path,default=Path('data'))
    p.add_argument('--platform',default='CPU',choices=['CPU','OpenCL','Reference'])
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.mode in ['all','toy']:run_toy(args.out)
    if args.mode in ['all','md']:run_md(args.out,args.platform)
    if args.mode in ['all','md','analyze']:analyze_md(args.out)
