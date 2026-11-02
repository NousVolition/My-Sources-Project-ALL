"""Resumable matched-stretch convergence, perturbation and parameter suite."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import os
import time
import traceback
import numpy as np
from numerics import Flow

HERE = Path(__file__).resolve().parent
END = .40
OUTPUT_DT = .01


def save(path, value):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    # Windows readers/indexers may briefly deny replacement of an open JSON file.
    # Retry the metadata write; no integration state or numerical setting changes.
    for attempt in range(100):
        try:
            tmp.replace(path)
            return
        except PermissionError:
            if attempt == 99:
                raise
            time.sleep(.1)


def jobs():
    all_jobs = []
    for n in (64, 128, 256):
        for half in (False, True):
            all_jobs.append(dict(id=f'baseline-n{n}-'+('half' if half else 'base'),
                                 n=n, family='convergence', nu=.001, half=half))
    for n in (64, 128):
        for amplitude in (.001, .005, .01):
            for seed in (101, 202, 303):
                all_jobs.append(dict(id=f'perturb-n{n}-a{amplitude}-s{seed}', n=n,
                    family='perturbation', nu=.001, amplitude=amplitude, seed=seed,
                    dependency=f'baseline-n{n}-base'))
        for nu in (.01, .002, .0005, 0.):
            all_jobs.append(dict(id=f'viscosity-n{n}-nu{nu}', n=n, family='viscosity', nu=nu))
        for name, values in (('radius', (.18, .22)), ('spin_factor', (.9, 1.1)),
                             ('strain_factor', (.9, 1.1)), ('angle_degrees', (-5., 5.))):
            for value in values:
                all_jobs.append(dict(id=f'{name}-n{n}-v{value}', n=n, family=name,
                                     nu=.001, **{name: value}))
    return all_jobs


def hashes():
    return {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('numerics.py', 'run_suite.py')}


def run_job(job):
    folder = HERE/'runs'/job['id']; folder.mkdir(parents=True, exist_ok=True)
    result_path = folder/'result.json'
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if result.get('status') == 'complete': return result
    f = Flow(job['n'], job['nu'], workers=4 if job['n'] == 256 else 2)
    params = {name: job[name] for name in ('radius', 'spin_factor', 'strain_factor', 'angle_degrees') if name in job}
    saved_initial = HERE/f'initial-n{job["n"]}.npy'
    h = np.load(saved_initial) if not params else f.initial(**params)
    if job['family'] == 'perturbation': h = f.perturb(h, job['amplitude'], job['seed'])
    cached = f.rhs(h)
    raw_dt = min(.04*f.dx/max(cached[2], 1e-12),
                 .1*f.dx*f.dx/(6*f.nu) if f.nu > 0 else np.inf)
    # All trajectories at a grid use exactly the same output times. dt/2 doubles
    # the number of steps between outputs, not the physical simulation duration.
    substeps = int(np.ceil(OUTPUT_DT/raw_dt))*(2 if job.get('half') else 1)
    dt = OUTPUT_DT/substeps
    initial = f.observe(h, 0)
    normalizer = np.sqrt(2*initial['energy'])
    start_output = 0; accum = np.zeros(5); series = []; max_cfl = 0.
    if result_path.exists():
        previous = json.loads(result_path.read_text())
        if previous.get('source_hashes') != hashes(): raise RuntimeError('Saved solver hashes changed; refusing silent resume')
        if previous['dt'] != dt: raise RuntimeError('Saved timestep changed')
        series = previous['series']; accum = np.array(previous['accum'])
        start_output = previous['last_output']+1
        h = np.load(folder/previous['field'])
        cached = f.rhs(h); max_cfl = previous['max_speed_cfl']
    start_wall = time.perf_counter()
    for output_index in range(start_output, 41):
        t = output_index*OUTPUT_DT
        if output_index:
            for sub in range(substeps):
                before_W = cached[1][0]
                updated, stage_integral, stage_speed = f.step(h, dt, cached)
                if not np.isfinite(updated).all():
                    raise FloatingPointError(f'Nonfinite field at t={(output_index-1)*OUTPUT_DT+(sub+1)*dt}')
                h = updated
                cached = f.rhs(h)
                stage_integral[0] = .5*dt*(before_W+cached[1][0])
                accum += stage_integral
                max_cfl = max(max_cfl, stage_speed*dt/f.dx)
                if max_cfl > .5: raise FloatingPointError('Speed CFL exceeded predeclared 0.5 stop threshold')
                if sub % 50 == 0:
                    save(folder/'progress.json', {'pid': os.getpid(), 'status': 'running', 'job': job['id'],
                         't': (output_index-1)*OUTPUT_DT+(sub+1)*dt, 'end': END,
                         'step': (output_index-1)*substeps+sub+1, 'total_steps': 40*substeps,
                         'updated_epoch': time.time(), 'Wmax': float(cached[1][0]),
                         'elapsed_this_process_seconds': time.perf_counter()-start_wall})
        row = f.observe(h, t, accum, initial)
        if job['family'] == 'perturbation':
            baseline = np.load(HERE/'runs'/job['dependency']/f'field-{output_index:03d}.npy')
            delta = h-baseline
            baseline_initial = json.loads((HERE/'runs'/job['dependency']/'result.json').read_text())['series'][0]
            baseline_norm = np.sqrt(2*baseline_initial['energy'])
            row['D_relative_initial_baseline_l2'] = np.sqrt(f.inner(delta, delta))/baseline_norm
            del baseline, delta
        row['resolved_screen'] = bool(row['high_band_energy_fraction'] < .001 and
            row['high_band_enstrophy_fraction'] < .01 and
            (row['width_at_global_peak']['minimum_chord_cells'] or 0) >= 6)
        scaled_div = row['divergence_max']*f.dx/max(row['peak_speed'], 1e-12)
        if scaled_div > 1e-10: raise FloatingPointError(f'Incompressibility check failed: {scaled_div}')
        if abs(row['energy_budget_relative_residual']) > .005:
            # Preserve the failed row and field before ending this member; other jobs continue.
            row['hard_failure'] = 'Energy balance residual exceeds 0.5% of initial energy'
        field_name = f'field-{output_index:03d}.npy'
        temp_field = folder/(field_name+'.tmp')
        with temp_field.open('wb') as handle: np.save(handle, h)
        temp_field.replace(folder/field_name)
        series.append(row)
        status = 'failed_check' if 'hard_failure' in row else ('complete' if output_index == 40 else 'running')
        result = {'status': status, 'job': job, 'dt': dt, 'steps_per_output': substeps,
                  'end_time': END, 'last_output': output_index, 'field': field_name,
                  'accum': accum.tolist(), 'source_hashes': hashes(), 'series': series,
                  'max_speed_cfl': max_cfl, 'updated_epoch': time.time(),
                  'elapsed_this_process_seconds': time.perf_counter()-start_wall}
        save(result_path, result)
        print(f'{job["id"]}: t={t:.2f} W={row["Wmax"]:.6g} I={row["I"]:.6g} '
              f'tailZ={row["high_band_enstrophy_fraction"]:.3%} status={status}', flush=True)
        if status == 'failed_check': return result
    save(folder/'progress.json', {'pid': os.getpid(), 'status': 'complete', 'job': job['id'], 't': END, 'updated_epoch': time.time()})
    return result


def protocol():
    return {'created_epoch': time.time(), 'end_time': END, 'output_dt': OUTPUT_DT,
            'domain_side': 6, 'baseline_viscosity': .001, 'external_force': 0,
            'initial_source': '../fixed-cutoff-run/source/grids_compared.py:build',
            'initial_audit': 'initial-audit.json', 'implementation_validation': 'implementation-check.json',
            'time_method': 'Heun, same as supplied; output intervals are divided into equal steps with starting speed CFL <= 0.04',
            'integrals': 'Wmax is trapezoid-integrated at every accepted step; energy/enstrophy rates use both Heun stages',
            'perturbation_norm': 'L2(delta u0)/L2(u0); zero-mean divergence-free modes, seeds 101,202,303 at n64 and n128',
            'parameter_design': 'One factor at a time at n64 and n128; radius +/-10% at fixed analytic tube spin, spin +/-10%, strain +/-10%, strain-axis tilt +/-5 degrees',
            'resolution_screen': 'Flags, not proof: outer 20% of retained component modes carry <0.1% energy and <1% enstrophy; minimum transverse half-peak chord >=6 cells',
            'convergence_screen': 'Report actual relative errors; exploratory screens: <5% finest grid curve difference and <1% dt-halving difference; require errors decrease with refinement',
            'hard_stops_per_job': 'Nonfinite fields, speed CFL>0.5, scaled divergence>1e-10, energy budget residual>0.5%; preserve evidence and continue unrelated jobs',
            'publication': 'Local calculations and factual reports. No conjectures or speculative project connections are authorized for GitHub.',
            'jobs': jobs(), 'source_hashes': hashes()}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if json.loads((HERE/'implementation-check.json').read_text())['status'] != 'passed':
        raise RuntimeError('Implementation validation must pass before evolution')
    state_path = HERE/'task-state.json'
    lock = HERE/'controller.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError('Controller lock exists. Inspect its process before resuming; never run a duplicate.')
    os.write(descriptor, str(os.getpid()).encode()); os.close(descriptor)
    try:
        if not (HERE/'protocol.json').exists(): save(HERE/'protocol.json', protocol())
        else:
            if json.loads((HERE/'protocol.json').read_text())['source_hashes'] != hashes():
                raise RuntimeError('Protocol source hashes changed')
        pending = jobs(); completed = {}; active = {}; failures = {}
        for job in list(pending):
            path = HERE/'runs'/job['id']/'result.json'
            if path.exists():
                result = json.loads(path.read_text())
                if result['status'] == 'complete': completed[job['id']] = result['series'][-1]['t']; pending.remove(job)
                elif result['status'] == 'failed_check': failures[job['id']] = 'Saved material check failure'; pending.remove(job)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            while pending or active:
                while len(active) < args.workers:
                    choices = [j for j in pending if not j.get('dependency') or j['dependency'] in completed]
                    # One convergence member and one smaller auxiliary member when available.
                    # This lets the perturbation/parameter work advance during expensive refinement.
                    if any(j['family'] == 'convergence' for j in active.values()):
                        auxiliary = [j for j in choices if j['family'] != 'convergence']
                        if auxiliary: choices = auxiliary
                    if not choices: break
                    job = choices[0]; pending.remove(job)
                    active[pool.submit(run_job, job)] = job
                    print('START '+job['id'], flush=True)
                for future in list(active):
                    if not future.done(): continue
                    job = active.pop(future)
                    try:
                        result = future.result()
                        if result['status'] == 'complete': completed[job['id']] = result['series'][-1]['t']
                        else: failures[job['id']] = result['status']
                    except Exception:
                        failures[job['id']] = traceback.format_exc()
                        print(f'FAILED {job["id"]}: {failures[job["id"]]}', flush=True)
                        save(HERE/'runs'/job['id']/'failure.json', {'error': failures[job['id']], 'updated_epoch': time.time()})
                blocked = [j for j in pending if j.get('dependency') in failures]
                for j in blocked:
                    failures[j['id']] = 'Baseline failed; comparison cannot run'; pending.remove(j)
                state = {'status': 'running' if pending or active else ('complete_with_failures' if failures else 'complete'),
                         'pid': os.getpid(), 'updated_epoch': time.time(), 'planned': len(jobs()),
                         'completed': completed, 'active': [j['id'] for j in active.values()],
                         'pending': [j['id'] for j in pending], 'failures': failures}
                save(state_path, state)
                if pending or active: time.sleep(2)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__': main()
