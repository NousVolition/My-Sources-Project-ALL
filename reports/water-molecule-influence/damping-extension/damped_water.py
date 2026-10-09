"""Uniform drag around constrained velocity Verlet, with matched initial states.

This water implementation stores on-step velocities (RATTLE), not staggered
velocities. damped_leapfrog.py separately implements literal half-step leapfrog.
No thermal noise: gamma > 0 intentionally removes energy from the water.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
import numpy as np
import openmm as mm
from openmm import unit
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'response-extension'))
from impulse_response import Probe, TIMES, save


def damped_integrator(dt, gamma):
    if gamma < 0 or dt <= 0:
        raise ValueError('gamma >= 0 and dt > 0 required')
    integ = mm.CustomIntegrator(dt*unit.picosecond)
    integ.addGlobalVariable('gamma', gamma)
    integ.addPerDofVariable('unconstrained_position', 0)
    integ.addUpdateContextState()
    # Exact half drag, one constrained velocity-Verlet step, exact half drag.
    integ.addComputePerDof('v', 'exp(-0.5*gamma*dt)*v')
    integ.addComputePerDof('v', 'v+0.5*dt*f/m')
    integ.addComputePerDof('x', 'x+dt*v')
    integ.addComputePerDof('unconstrained_position', 'x')
    integ.addConstrainPositions()
    integ.addComputePerDof('v', 'v+0.5*dt*f/m+(x-unconstrained_position)/dt')
    integ.addConstrainVelocities()
    integ.addComputePerDof('v', 'exp(-0.5*gamma*dt)*v')
    integ.setConstraintTolerance(1e-8)
    return integ


class DampedProbe(Probe):
    def __init__(self, system, platform_name, gamma, dt=.001):
        self.dt, self.gamma = dt, gamma
        self.steps = np.rint(TIMES/dt).astype(int)
        assert np.allclose(self.steps*dt, TIMES)
        self.integrator = damped_integrator(dt, gamma)
        platform = mm.Platform.getPlatformByName(platform_name)
        props = {'Precision': 'double'} if platform_name == 'OpenCL' else {}
        self.context = mm.Context(system, self.integrator, platform, props)
        self.mass = np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                              for i in range(system.getNumParticles())]).reshape(-1, 3)
        self.n = len(self.mass)
        self.weights = self.mass/self.mass.sum(axis=1, keepdims=True)

    def free_response(self):
        if self.gamma == 0:
            return TIMES.copy()
        # Exact discrete COM response of this split algorithm. The continuum
        # result (1-exp(-gamma*t))/gamma differs by O(dt^2).
        return self.dt*np.exp(-self.gamma*self.dt/2)*np.expm1(-self.gamma*TIMES)/np.expm1(-self.gamma*self.dt)

    def sources(self, x, v, source_ids, delta_v=.01, progress=False):
        pairs = []
        start = time.time()
        for count, i in enumerate(source_ids):
            response = np.empty((len(TIMES), self.n, 3, 3))
            for axis in range(3):
                kick = np.zeros_like(v)
                kick[:, :, axis] = -delta_v/(self.n-1)
                kick[i, :, axis] = delta_v
                plus, _ = self.trajectory(x, v+kick)
                minus, _ = self.trajectory(x, v-kick)
                derivative = (plus-minus)/(2*delta_v)
                derivative[:, :, axis] += self.free_response()[:, None]/(self.n-1)
                response[:, :, :, axis] = derivative*(self.n-1)/self.n
            norm = np.linalg.norm(response, axis=(2, 3))
            norm[:, i] = 0
            pairs.append(norm)
            if progress and ((count+1) % 36 == 0 or count+1 == len(source_ids)):
                print(f'gamma={self.gamma:g}: {count+1}/{len(source_ids)} sources, {time.time()-start:.1f}s', flush=True)
        return np.array(pairs)

    def close(self):
        del self.context
        del self.integrator


def relative(base, other):
    return {'L2_by_time': np.linalg.norm(other-base, axis=0)/np.linalg.norm(base, axis=0),
            'max_individual_relative_by_time': np.max(abs(other-base)/np.maximum(abs(base), 1e-30), axis=0)}


def run(out, platform, gammas, replica):
    out.mkdir(parents=True, exist_ok=True)
    source_dir = HERE.parent/'response-extension'/'data'
    input_path = source_dir/f'response_replica_{replica}.npz'
    source = np.load(input_path)
    x, v = source['positions_nm'], source['velocities_nm_per_ps']
    force, undamped = source['force_scores'], source['scores_ps']
    system = mm.XmlSerializer.deserialize((source_dir/'probe_system.xml').read_text())
    n = len(x)
    selected = np.unique(np.r_[np.argsort(force)[[0, n//2, n-1]], undamped[:, -1].argmax()])
    protocol = {'gamma_ps_inverse': gammas, 'replica': replica, 'molecules': n, 'times_ps': TIMES,
                'dt_ps': .001, 'delta_v_nm_ps': .01, 'constraint_tolerance': 1e-8,
                'integrator': 'Symmetric exact half drag / RATTLE velocity Verlet / exact half drag; on-step velocities.',
                'equation': 'm dv/dt = F - m gamma v, with rigid constraints; no random force',
                'platform': platform, 'precision': 'double', 'openmm': mm.__version__,
                'source_sha256': hashlib.sha256(input_path.read_bytes()).hexdigest(),
                'input': f'../response-extension/data/response_replica_{replica}.npz',
                'initial_states': 'Exact same saved positions and on-step velocities for every gamma and dt.',
                'baseline': 'Previously saved gamma=0 all-source response, checked again on selected sources.',
                'controls_source_ids': selected,
                'scope': 'One conditional TIP3P state, external damping without thermal noise; not equilibrium water.',
                'response': 'Same paired COM derivative as previous study, replacing ballistic t with discrete free-drag displacement G_gamma,h(t).'}
    save(out/'protocol.json', protocol)
    baseline_probe = DampedProbe(system, platform, 0)
    replay = baseline_probe.sources(x, v, selected).sum(axis=2)
    baseline_control = relative(undamped[selected], replay)
    assert np.max(baseline_control['L2_by_time']) < 1e-4
    baseline_probe.close()
    results = []
    for gamma in gammas:
        probe = DampedProbe(system, platform, gamma)
        baseline, energies = probe.trajectory(x, v, energy=True)
        pairs = probe.sources(x, v, np.arange(n), progress=True)
        scores = pairs.sum(axis=2)
        assert np.isfinite(pairs).all() and (pairs >= 0).all()
        half_kick = probe.sources(x, v, selected, .005).sum(axis=2)
        permutation = np.random.default_rng(20261015).permutation(n)
        permuted = probe.sources(x[permutation], v[permutation], np.argsort(permutation)[selected]).sum(axis=2)
        probe.close()
        small = DampedProbe(system, platform, gamma, .0005)
        half_step = small.sources(x, v, selected).sum(axis=2)
        small.close()
        free_system = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
        for i in reversed(range(free_system.getNumForces())):
            free_system.removeForce(i)
        free = DampedProbe(free_system, platform, gamma)
        free_scores = free.sources(x, v, selected).sum(axis=2)
        free.close()
        np.savez_compressed(out/f'water_gamma_{gamma:g}.npz', positions_nm=x, velocities_nm_per_ps=v,
                            force_scores=force, undamped_scores_ps=undamped,
                            pair_response_norms_ps=pairs, scores_ps=scores, baseline_COM_nm=baseline,
                            baseline_energies_kJ_mol=energies, control_sources=selected,
                            half_kick_scores_ps=half_kick, half_step_scores_ps=half_step,
                            free_scores_ps=free_scores, permutation=permutation, permuted_scores_ps=permuted)
        rows = []
        for k, t in enumerate(TIMES):
            rows.append({'time_ps': t, 'mean_score_ps': scores[:, k].mean(),
                         'mean_over_undamped': scores[:, k].mean()/undamped[:, k].mean(),
                         'spearman_vs_undamped': spearmanr(scores[:, k], undamped[:, k]).statistic,
                         'spearman_vs_force': spearmanr(scores[:, k], force).statistic,
                         'leader': int(scores[:, k].argmax()),
                         'force_leader_rank': int(np.flatnonzero(np.argsort(-scores[:, k]) == force.argmax())[0]+1)})
        record = {'gamma_ps_inverse': gamma, 'rows': rows,
                  'energy_change_kJ_mol_per_water': (energies-energies[0])/n,
                  'half_kick': relative(scores[selected], half_kick),
                  'half_timestep': relative(scores[selected], half_step),
                  'permutation': relative(scores[selected], permuted),
                  'no_interactions_max_score_ps': free_scores.max()}
        results.append(record)
        save(out/'results.json', {'zero_damping_replay': baseline_control, 'runs': results})
        print(json.dumps(record, default=lambda a: a.tolist(), indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gamma', type=float, nargs='+', default=[5, 50], help='Drag rates in inverse ps (nonnegative).')
    parser.add_argument('--replica', type=int, choices=[0, 1, 2], default=0)
    parser.add_argument('--platform', choices=['Reference', 'OpenCL'], default='Reference')
    parser.add_argument('--out', type=Path, default=HERE/'data')
    args = parser.parse_args()
    if any(not np.isfinite(g) or g < 0 for g in args.gamma) or len(set(args.gamma)) != len(args.gamma):
        parser.error('gamma values must be finite, nonnegative and distinct')
    run(args.out, args.platform, args.gamma, args.replica)
