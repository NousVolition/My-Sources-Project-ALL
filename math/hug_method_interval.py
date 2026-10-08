"""Longer Euler/Heun comparison from one identical saved 384^3 state.

Uses the verified time integrators without altering the original solver.
Each branch records energy at common times and saves reusable checkpoints.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from box_experiment import diagnostics
from continue_hug_refinement import file_hash, numerical_source_hash
from hug_refinement import compare_fields
from hug_shape import energy_shape
from hug_time_method import difference, energy, evolve
from navier import PurePythonNavierStokes3D


def observe(u, box, nu, dt):
    """Retain the established velocity, gradient and shape definitions."""
    n=u.shape[1]
    sim=PurePythonNavierStokes3D(N=4)
    sim.N,sim.dx=n,box/n
    sim.u,sim.v,sim.w=u
    measured=diagnostics(sim)
    measured.update(energy_shape(u,box))
    speed2=np.sum(u*u,axis=0)
    measured.update(max_speed=float(np.sqrt(speed2.max())),
        advective_cfl=float(dt/sim.dx*np.max(np.sum(np.abs(u),axis=0))),
        diffusion_number=6*nu*dt/sim.dx**2)
    if not all(np.isfinite(v) for v in measured.values()) or measured['max_div']>1e-10:
        raise RuntimeError('Nonfinite field or incompatible divergence at checkpoint')
    return measured


def write_json(path, value):
    temporary=path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def read_point(path, provenance):
    saved=json.loads(path.read_text())
    if saved['provenance']!=provenance:
        raise ValueError('Existing checkpoint has different inputs or source')
    arrays=path.with_suffix('.npz')
    if file_hash(arrays)!=saved['arrays_sha256']:
        raise ValueError('Checkpoint array hash changed')
    return saved,arrays


def run_branch(initial,start,stops,dt,method,cache,source,box=6.,nu=.01,sample_interval=.001):
    """Resume matching completed checkpoints; do not re-project on restart."""
    stops=tuple(stops)
    if (method not in ('euler','heun') or not np.isfinite(dt) or dt<=0
            or not np.isfinite(sample_interval) or sample_interval<=0
            or not np.isfinite(start) or not stops
            or not np.isfinite(box) or box<=0 or not np.isfinite(nu) or nu<0
            or abs(round(sample_interval/dt)*dt-sample_interval)>1e-12
            or sample_interval<dt
            or any(not np.isfinite(t) or t<=a or abs(round((t-start)/sample_interval)*sample_interval-(t-start))>1e-12
                   for a,t in zip((start,)+stops[:-1],stops))):
        raise ValueError('Require increasing aligned stops and a positive aligned time step')
    initial=np.asarray(initial)
    n=initial.shape[1]
    provenance=dict(source=source,N=n,start_time=start,dt=dt,method=method,box=box,nu=nu,
        external_force=0.,sample_interval=sample_interval,
        runner_sha256=file_hash(__file__),
        integrator_sha256=file_hash(Path(__file__).with_name('hug_time_method.py')))
    cache=Path(cache);cache.mkdir(parents=True,exist_ok=True)
    def path_for(t):return cache/f'{method}-dt{dt:g}-t{t:g}.json'
    # Verify every retained point, even when the final one already exists.
    existing={}
    for stop in stops:
        path=path_for(stop)
        if path.exists():existing[stop]=read_point(path,provenance)
    if stops[-1] in existing:
        print(f'Reuse {method} dt={dt:g} through {stops[-1]:g}',flush=True)
        return [dict(file=path_for(t).name,**existing[t][0]) for t in stops]
    if existing:
        current_time=max(existing)
        saved,arrays_path=existing[current_time]
        with np.load(arrays_path,allow_pickle=False) as arrays:
            current=arrays['final']
        rows=list(saved['rows'])
        if not np.isclose(energy(current,box),rows[-1]['energy'],rtol=2e-12,atol=1e-13):
            raise ValueError('Restored energy differs from checkpoint')
        cumulative={k:saved[k] for k in ('work_quadrature','temporal_energy_defect','stable_energy_change',
                                        'max_identity_residual','max_stage_divergence')}
        print(f'Resume {method} dt={dt:g} at {current_time:g}',flush=True)
    else:
        current_time=start;current=initial
        measured=observe(current,box,nu,dt)
        cumulative=dict(work_quadrature=0.,temporal_energy_defect=0.,stable_energy_change=0.,
                        max_identity_residual=0.,max_stage_divergence=measured['max_div'])
        rows=[dict(time=start,energy=measured['energy'],**cumulative)]
    baseline=rows[0]['energy']
    began=time.perf_counter()
    for stop in stops:
        if stop<=current_time:continue
        chunks=round((stop-current_time)/sample_interval)
        for _ in range(chunks):
            current,ledger=evolve(current,sample_interval,round(sample_interval/dt),method,box,nu)
            current_time=round(current_time+sample_interval,12)
            for k in ('work_quadrature','temporal_energy_defect','stable_energy_change'):
                cumulative[k]+=ledger[k]
            for k in ('max_identity_residual','max_stage_divergence'):
                cumulative[k]=max(cumulative[k],ledger[k])
            row=dict(time=current_time,energy=ledger['energy_after'],**cumulative)
            if abs(row['energy']-baseline-row['work_quadrature']-row['temporal_energy_defect'])>5e-12*max(1.,baseline):
                raise RuntimeError('Accumulated energy identity failed')
            rows.append(row)
            print(f'{method} dt={dt:g} t={current_time:.3f}: E={row["energy"]:.9f} '
                  f'({time.perf_counter()-began:.1f}s)',flush=True)
        measured=observe(current,box,nu,dt)
        path=path_for(stop);arrays_path=path.with_suffix('.npz')
        temporary=arrays_path.with_suffix('.npz.tmp')
        with temporary.open('wb') as handle:
            np.savez_compressed(handle,final=current)
        temporary.replace(arrays_path)
        saved=dict(provenance=provenance,end_time=stop,steps=round((stop-start)/dt),
            rows=list(rows),**cumulative,observations=measured,
            energy_change=measured['energy']-baseline,
            total_identity_residual=measured['energy']-baseline-cumulative['work_quadrature']-cumulative['temporal_energy_defect'],
            arrays_sha256=file_hash(arrays_path))
        write_json(path,saved)
        existing[stop]=(saved,arrays_path)
        print(f'Saved {path.name}; finite fields, divergence={measured["max_div"]:.3e}',flush=True)
    return [dict(file=path_for(t).name,**existing[t][0]) for t in stops]


def study(input_path,checkpoint_cache,cache,out):
    input_path,checkpoint_cache,cache,out=map(Path,(input_path,checkpoint_cache,cache,out))
    source=json.loads(input_path.read_text())
    fingerprint=numerical_source_hash()
    if source['external_force']!=0 or source['end_time']!=.4 or source['numerical_source_sha256']!=fingerprint:
        raise ValueError('Require the verified saved-flow record at 0.40')
    origin=next(r for r in source['row_sources'] if (r['N'],r['dt'],r['time'])==(384,.0005,.4))
    path=checkpoint_cache/origin['file']
    if file_hash(path)!=origin['summary_sha256'] or file_hash(path.with_suffix('.npz'))!=origin['arrays_sha256']:
        raise ValueError('Original checkpoint hash changed')
    saved=json.loads(path.read_text())
    if (saved['N'],saved['box'],saved['nu'],saved['sigma'],saved['P_U'],saved['source_sha256'])!=(384,6.,.01,0.,0.,fingerprint):
        raise ValueError('Original equation settings changed')
    with np.load(path.with_suffix('.npz'),allow_pickle=False) as arrays:
        initial=arrays['final']
    initial_observations=observe(initial,6.,.01,.001)
    reference=next(r for r in source['rows'] if (r['N'],r['dt'],r['time'])==(384,.0005,.4))
    if not np.isclose(initial_observations['energy'],reference['energy'],rtol=2e-12,atol=1e-13):
        raise ValueError('Starting energy differs from prior record')
    common=dict(checkpoint=origin,input_file=input_path.name,input_sha256=file_hash(input_path),
                numerical_source_sha256=fingerprint)
    configurations=(('euler',.001),('euler',.0005),('heun',.001),('heun',.0005))
    runs=[]
    for method,dt in configurations:
        points=run_branch(initial,.4,(.405,.410),dt,method,cache,common)
        runs.append(dict(method=method,dt=dt,checkpoints=points,**{k:v for k,v in points[-1].items()
            if k not in ('provenance','file')}))
    del initial
    comparisons=[]
    for method in ('euler','heun'):
        full,half=[r for r in runs if r['method']==method]
        for j,stop in enumerate((.405,.410)):
            with np.load((cache/full['checkpoints'][j]['file']).with_suffix('.npz'),allow_pickle=False) as a:
                coarse=a['final']
            with np.load((cache/half['checkpoints'][j]['file']).with_suffix('.npz'),allow_pickle=False) as b:
                fine=b['final']
            comparison=dict(method=method,time=stop,**difference(coarse,fine,6.),**compare_fields(coarse,fine,6.))
            comparisons.append(comparison)
            del coarse,fine
            print(f'Compared {method} at {stop:g}: energy difference={comparison["energy_difference"]:.6e}, '
                  f'velocity difference={comparison["velocity_relative_l2"]:.6e}',flush=True)
    result=dict(start_time=.4,end_time=.41,grid=384,input_history_dt=.0005,external_force=0.,box=6.,nu=.01,
        original_trajectory_preserved_through=.4,common_source=common,
        numerical_source_sha256=fingerprint,analysis_source_sha256=file_hash(__file__),
        integrator_sha256=file_hash(Path(__file__).with_name('hug_time_method.py')),
        initial_observations=initial_observations,runs=runs,comparisons=comparisons,
        limits=['All four branches start from the identical saved 384^3 half-step Euler state at time 0.40; prior histories were not rerun.',
            'These separate saved method controls cover 0.40 to 0.41. The earlier three-run grid-comparison series is preserved through 0.40.',
            'Work quadrature includes projected spatial transport and viscosity. Temporal defect is not the earlier viscosity-only interval gap.',
            'Full-step versus half-step differences measure sensitivity, not an exact error or rigorous bound.',
            'One grid and one additional finite interval do not establish long-time stability, continuum convergence or a Clay proof.'])
    out.parent.mkdir(parents=True,exist_ok=True);write_json(out,result)
    rows=[dict(method=r['method'],dt=r['dt'],**row) for r in runs for row in r['rows']]
    with out.with_suffix('.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('Longer method comparison complete.',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('math/results/hug-resolution.json'))
    parser.add_argument('--checkpoint-cache',type=Path,default=Path('scratch/hug-refinement'))
    parser.add_argument('--cache',type=Path,default=Path('scratch/hug-method-interval'))
    parser.add_argument('--out',type=Path,default=Path('math/results/hug-method-interval.json'))
    a=parser.parse_args();study(a.input,a.checkpoint_cache,a.cache,a.out)
