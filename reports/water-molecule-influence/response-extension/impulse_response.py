"""Finite-time follow-up to the positional water-influence experiment.

From the parent experiment directory:
  python response-extension/impulse_response.py --platform OpenCL
Use --platform Reference for a portable, slower double-precision calculation.
The saved inputs come from the original three TIP3P trajectories; no new water
model or independent long equilibration is claimed.
"""
from pathlib import Path
import argparse, csv, hashlib, json, time
import numpy as np
from scipy.stats import spearmanr
import openmm as mm
from openmm import unit

HERE=Path(__file__).resolve().parent
TIMES=np.array([.005,.02,.1,.2])
SEED=20261011


def native(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)


def save(path,value):path.write_text(json.dumps(value,indent=2,default=native),encoding='utf-8')


class Probe:
    def __init__(self,system,platform_name,dt=.001):
        self.dt=dt
        self.steps=np.rint(TIMES/dt).astype(int)
        assert np.allclose(self.steps*dt,TIMES)
        # RATTLE velocity Verlet stores velocities at the same time as positions.
        # This makes equal-state comparisons across timesteps well defined.
        self.integrator=mm.CustomIntegrator(dt*unit.picoseconds)
        self.integrator.addPerDofVariable('unconstrained_position',0)
        self.integrator.addUpdateContextState()
        self.integrator.addComputePerDof('v','v+0.5*dt*f/m')
        self.integrator.addComputePerDof('x','x+dt*v')
        self.integrator.addComputePerDof('unconstrained_position','x')
        self.integrator.addConstrainPositions()
        self.integrator.addComputePerDof('v','v+0.5*dt*f/m+(x-unconstrained_position)/dt')
        self.integrator.addConstrainVelocities()
        self.integrator.setConstraintTolerance(1e-8)
        platform=mm.Platform.getPlatformByName(platform_name)
        props={'Precision':'double'} if platform_name=='OpenCL' else {}
        self.context=mm.Context(system,self.integrator,platform,props)
        self.mass=np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                           for i in range(system.getNumParticles())]).reshape(-1,3)
        self.n=len(self.mass)
        self.weights=self.mass/self.mass.sum(axis=1,keepdims=True)

    def reset(self,x,v):
        self.context.setPositions(x.reshape(-1,3)*unit.nanometer)
        self.context.setVelocities(v.reshape(-1,3)*unit.nanometer/unit.picosecond)
        self.context.setTime(0)
        self.context.setStepCount(0)

    def energy(self):
        s=self.context.getState(getEnergy=True)
        return (s.getKineticEnergy()+s.getPotentialEnergy()).value_in_unit(unit.kilojoule_per_mole)

    def trajectory(self,x,v,energy=False):
        self.reset(x,v)
        positions=[];energies=[self.energy()] if energy else []
        last=0
        for step in self.steps:
            self.integrator.step(int(step-last));last=step
            state=self.context.getState(getPositions=True,getEnergy=energy)
            xyz=np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(self.n,3,3)
            positions.append(np.einsum('na,nac->nc',self.weights,xyz))
            if energy:energies.append((state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoule_per_mole))
        return np.array(positions),np.array(energies)

    def sources(self,x,v,source_ids,delta_v=.01,progress=False):
        """One score per source and horizon, with direction-neutral 3x3 blocks.

        Add +dv to source and -dv/(N-1) to every other molecule; reverse all
        kicks for the minus trajectory. B is the central positional response
        after subtracting the direct compensating ballistic displacement, then
        rescaling by (N-1)/N to match a single-source impulse derivative.
        """
        pairs=[];start=time.time()
        for count,i in enumerate(source_ids):
            response=np.empty((len(TIMES),self.n,3,3))
            for axis in range(3):
                kick=np.zeros_like(v)
                kick[:,:,axis]=-delta_v/(self.n-1)
                kick[i,:,axis]=delta_v
                plus,_=self.trajectory(x,v+kick)
                minus,_=self.trajectory(x,v-kick)
                derivative=(plus-minus)/(2*delta_v)
                derivative[:,:,axis]+=TIMES[:,None]/(self.n-1)
                response[:,:,:,axis]=derivative*(self.n-1)/self.n
            norm=np.linalg.norm(response,axis=(2,3))
            norm[:,i]=0
            pairs.append(norm)
            if progress and ((count+1)%36==0 or count+1==len(source_ids)):
                print(f'  {count+1}/{len(source_ids)} sources; {time.time()-start:.1f}s',flush=True)
        # dimensions: source, time, target
        return np.array(pairs)


def rank_controls(base,other):
    return {'relative_L2_difference':np.linalg.norm(other-base)/np.linalg.norm(base),
            'max_relative_difference':np.max(np.abs(other-base)/np.maximum(base,1e-30)),
            'spearman_by_time':[spearmanr(base[:,k],other[:,k]).statistic for k in range(base.shape[1])]}


def run(data,out,platform_name):
    out.mkdir(parents=True,exist_ok=True)
    original_xml=(data/'water_system.xml').read_text()
    system=mm.XmlSerializer.deserialize(original_xml)
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i),mm.CMMotionRemover):system.removeForce(i)
    (out/'probe_system.xml').write_text(mm.XmlSerializer.serialize(system),encoding='utf-8')
    probe=Probe(system,platform_name)
    smaller_step=Probe(system,platform_name,.0005)
    free_system=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
    for i in reversed(range(free_system.getNumForces())):free_system.removeForce(i)
    free_probe=Probe(free_system,platform_name)
    n=probe.n;mass=probe.mass.sum(axis=1)
    assert np.allclose(mass,mass[0])
    protocol={'question':'Does snapshot force sensitivity predict finite-time influence after an equal velocity impulse?',
        'molecules':n,'snapshots':'50 ps production snapshot from each of the original three runs',
        'velocity_initialization':'New constrained Maxwell velocities at 300 K; seeds 20261011, 20261012, 20261013; remove overall center-of-mass velocity.',
        'positions':'Original saved positions projected to 1e-8 constraint tolerance; correction recorded. Molecular centers of mass are measured without wrapping.',
        'time_step_ps':.001,'constraint_tolerance':1e-8,'times_ps':TIMES,'delta_v_nm_per_ps':.01,
        'dynamics':'Short deterministic NVE probes with RATTLE velocity Verlet and on-step velocities; thermostat absent and CMMotionRemover removed.',
        'perturbation':'All atoms of source receive the same kick; all other molecules receive the opposite kick divided by N-1. Paired positive/negative kicks along each Cartesian axis.',
        'response':'Central difference of target COM positions; subtract known compensating ballistic displacement; multiply by (N-1)/N; sum Frobenius norms over targets excluding source.',
        'response_units':'ps (position / velocity); not an energy-transfer fraction',
        'controls':'Half kick on 10 force-rank-spaced sources per snapshot; half timestep on weakest, median and strongest force-score sources; zero interactions on same 3; same-trajectory replay; full-state label permutation on same 3 in snapshot 0.',
        'platform':platform_name,'precision':'double','openmm':mm.__version__,'numpy':np.__version__,
        'source_system_sha256':hashlib.sha256((data/'water_system.xml').read_bytes()).hexdigest(),
        'seeds':[SEED+i for i in range(3)],
        'caution':'Three conditional microstates, not a longer trajectory study, new force-field test, or test of oxygen positions independently of orientations.'}
    save(out/'protocol.json',protocol)
    summary=[];controls=[];rows=[]
    for rep in range(3):
        source=np.load(data/f'md_replica_{rep}.npz')
        x=source['positions_nm'][-1].copy()
        force_score=source['force_scores'][-1].copy()
        probe.context.setPositions(x.reshape(-1,3)*unit.nanometer)
        probe.context.applyConstraints(1e-8)
        corrected=np.array(probe.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(n,3,3)
        constraint_correction=np.max(np.abs(corrected-x))
        x=corrected
        probe.context.setVelocitiesToTemperature(300*unit.kelvin,SEED+rep)
        probe.context.applyVelocityConstraints(1e-8)
        v=np.array(probe.context.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)).reshape(n,3,3)
        v-=np.sum(probe.mass[:,:,None]*v,axis=(0,1))/probe.mass.sum()
        baseline,energies=probe.trajectory(x,v,energy=True)
        replay,_=probe.trajectory(x,v)
        print(f'Run {rep+1}: all {n} source molecules, paired impulses',flush=True)
        pairs=probe.sources(x,v,np.arange(n),progress=True)
        scores=pairs.sum(axis=2)
        assert np.isfinite(scores).all() and np.all(scores>0)
        order=np.argsort(force_score)
        selected=order[np.rint(np.linspace(0,n-1,10)).astype(int)]
        minimal=order[[0,n//2,n-1]]
        print(f'Run {rep+1}: convergence and noninteracting controls',flush=True)
        half_kick=probe.sources(x,v,selected,.005).sum(axis=2)
        half_time=smaller_step.sources(x,v,minimal).sum(axis=2)
        free=free_probe.sources(x,v,minimal).sum(axis=2)
        control={'replica':rep,'initial_constraint_correction_max_nm':constraint_correction,
            'same_initial_state_replay_max_COM_difference_nm':np.max(abs(replay-baseline)),
            'energy_change_kJ_mol':energies-energies[0],
            'max_energy_drift_kJ_mol_per_water':np.max(abs(energies-energies[0]))/n,
            'half_kick_sources':selected,'half_kick':rank_controls(scores[selected],half_kick),
            'half_timestep_sources':minimal,'half_timestep':rank_controls(scores[minimal],half_time),
            'no_interactions_max_score_ps':np.max(free),
            'no_interactions_max_relative_to_smallest_physical_score':np.max(free)/np.min(scores),
            'initial_total_momentum_norm':np.linalg.norm(np.sum(probe.mass[:,:,None]*v,axis=(0,1)))}
        if rep==0:
            permutation=np.random.default_rng(SEED+100).permutation(n)
            inverse=np.argsort(permutation)
            permuted=probe.sources(x[permutation],v[permutation],inverse[minimal]).sum(axis=2)
            control['permuted_sources']=inverse[minimal]
            control['permutation']=rank_controls(scores[minimal],permuted)
            np.savez_compressed(out/'permutation_control.npz',permutation=permutation,
                                sources=minimal,permuted_sources=inverse[minimal],permuted_scores=permuted)
        for k,t in enumerate(TIMES):
            s=scores[:,k];leader=int(s.argmax());force_leader=int(force_score.argmax())
            rank=int(np.flatnonzero(np.argsort(-s)==force_leader)[0]+1)
            predicted=force_score*t**3/(6*mass[0])
            row={'replica':rep,'time_ps':t,'response_mean_ps':s.mean(),
                'response_cv':s.std()/s.mean(),'response_max_min_ratio':s.max()/s.min(),
                'response_leader':leader,'initial_force_leader':force_leader,
                'initial_force_leader_response_rank':rank,
                'spearman_force_vs_response':spearmanr(force_score,s).statistic,
                'median_response_over_short_time_prediction':np.median(s/predicted)}
            rows.append(row)
        corrs=[[spearmanr(scores[:,i],scores[:,j]).statistic for j in range(len(TIMES))] for i in range(len(TIMES))]
        summary.append({'replica':rep,'response_rank_correlation_across_times':corrs})
        controls.append(control)
        np.savez_compressed(out/f'response_replica_{rep}.npz',positions_nm=x,
             velocities_nm_per_ps=v,baseline_COM_nm=baseline,baseline_energies_kJ_mol=energies,
             pair_response_norms_ps=pairs,scores_ps=scores,force_scores=force_score,
             half_kick_sources=selected,half_kick_scores_ps=half_kick,
             half_timestep_sources=minimal,half_timestep_scores_ps=half_time,free_scores_ps=free)
        save(out/'controls.json',controls)
        save(out/'results.json',{'rows':rows,'replicas':summary})
        print(f'Run {rep+1} saved; force/response correlations: '+', '.join(f'{r["spearman_force_vs_response"]:.3f}' for r in rows[-4:]),flush=True)
    with (out/'results.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('Finite-time response experiment complete.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-data',type=Path,default=HERE.parent/'data')
    parser.add_argument('--out',type=Path,default=HERE/'data')
    parser.add_argument('--platform',choices=['Reference','OpenCL'],default='Reference')
    args=parser.parse_args();run(args.source_data,args.out,args.platform)
