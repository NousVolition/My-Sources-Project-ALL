"""Trace one unchanged Euler/projection step on disposable checkpoint copies.

The actual solver supplies the predictor and projected velocities. Independent
slab sums account for transport work, viscous work, the combined quadratic
Euler contribution, and the measured projection change. Probe outputs never
replace saved checkpoints or extend the recorded trajectory.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import simpson, trapezoid

from continue_hug_refinement import file_hash, numerical_source_hash
from navier import PurePythonNavierStokes3D


def shifted(a, first, last, axis, offset):
    if axis == 0:
        return a[(np.arange(first,last)+offset) % a.shape[0]]
    return np.roll(a[first:last], -offset, axis=axis)


def step_energy(velocity, dt, box=6., nu=.01, return_velocity=False):
    """Actual one-step probe; E_after-E_before = transport + viscous + Q + P.

    Q=dt^2 ||F||_h^2/2 includes cross terms between transport and viscosity.
    P is the measured energy change of the existing orthogonal projection.
    Positive contributions increase energy, negative contributions decrease it.
    No sign-based rejection is made for transport work.
    """
    u=np.asarray(velocity)
    if (u.ndim != 4 or u.shape[0] != 3 or u.shape[1] < 4
            or u.shape[1:] != (u.shape[1],)*3 or not np.isrealobj(u)
            or not np.isfinite(dt) or dt <= 0 or not np.isfinite(box) or box <= 0
            or not np.isfinite(nu) or nu < 0):
        raise ValueError('Require a real cubic velocity, positive dt/box and nonnegative viscosity')
    n=u.shape[1]
    for c in u:
        for first in range(0,n,8):
            if not np.isfinite(c[first:first+8]).all():
                raise ValueError('Velocity must be finite')
    dx=box/n
    # Initialize defaults on a small grid; install the exact saved arrays as
    # inputs. The NumPy solver creates new outputs without changing these.
    sim=PurePythonNavierStokes3D(N=4)
    sim.N,sim.dx,sim.nu,sim.sigma=n,dx,nu,0.
    sim.u,sim.v,sim.w=u
    sim.S=np.zeros((n,n,n),dtype=float)
    sim.rho_epsilon=0.  # P_U=0, so this inactive scalar source is exactly zero.
    sim.step(dt,0.,backend='numpy')
    predictor=(sim.u,sim.v,sim.w)
    sums=dict(energy_before=0.,energy_predictor=0.,transport_work_rate=0.,
              viscous_work_rate=0.,rhs_norm_squared=0.)
    predictor_error,initial_div=0.,0.
    for first in range(0,n,8):
        last=min(first+8,n)
        divergence=np.zeros_like(u[0,first:last],dtype=float)
        for c in range(3):
            center=u[c,first:last]
            adv=np.zeros_like(center,dtype=float)
            lap=np.zeros_like(center,dtype=float)
            for d in range(3):
                plus=shifted(u[c],first,last,d,1)
                minus=shifted(u[c],first,last,d,-1)
                derivative=(plus-minus)/(2*dx)
                adv += u[d,first:last]*derivative
                lap += (plus+minus-2*center)/dx**2
                if c == d:
                    divergence += derivative
            rhs=-adv+nu*lap
            w=predictor[c][first:last]
            sums['energy_before'] += .5*dx**3*float(np.sum(center**2))
            sums['energy_predictor'] += .5*dx**3*float(np.sum(w**2))
            sums['transport_work_rate'] += dx**3*float(np.sum(-center*adv))
            sums['viscous_work_rate'] += nu*dx**3*float(np.sum(center*lap))
            sums['rhs_norm_squared'] += dx**3*float(np.sum(rhs**2))
            predictor_error=max(predictor_error,float(np.max(np.abs(w-(center+dt*rhs)))))
        initial_div=max(initial_div,float(np.max(np.abs(divergence))))
    projection=sim.project(backend='fft')
    after=(sim.u,sim.v,sim.w)
    energy_after,projection_change,correction_norm2,orthogonal_cross=0.,0.,0.,0.
    for c in range(3):
        for first in range(0,n,8):
            v=after[c][first:first+8]
            w=predictor[c][first:first+8]
            removed=w-v
            energy_after += .5*dx**3*float(np.sum(v**2))
            projection_change += .5*dx**3*float(np.sum((v-w)*(v+w)))
            correction_norm2 += dx**3*float(np.sum(removed**2))
            orthogonal_cross += dx**3*float(np.sum(v*removed))
    transport=dt*sums['transport_work_rate']
    viscous=dt*sums['viscous_work_rate']
    quadratic=.5*dt**2*sums['rhs_norm_squared']
    change=energy_after-sums['energy_before']
    accounted=transport+viscous+quadratic+projection_change
    metrics=dict(**sums,energy_after=energy_after,energy_change=change,
        transport_change=transport,viscous_change=viscous,
        quadratic_euler_change=quadratic,projection_change=projection_change,
        combined_quadratic_projection_change=quadratic+projection_change,
        excess_over_start_viscous_loss=change-viscous,
        accounted_energy_change=accounted,ledger_residual=change-accounted,
        predictor_ledger_residual=sums['energy_predictor']-sums['energy_before']-transport-viscous-quadratic,
        projection_ledger_residual=energy_after-sums['energy_predictor']-projection_change,
        projection_orthogonality_residual=orthogonal_cross,
        projection_norm_identity_residual=projection_change+.5*correction_norm2,
        removed_velocity_norm_squared=correction_norm2,
        predictor_max_difference=predictor_error,initial_max_div=initial_div,
        predictor_max_div=projection['before'],after_max_div=projection['after'])
    if not all(np.isfinite(v) for v in metrics.values()):
        raise RuntimeError('Nonfinite step diagnostic')
    tolerance=5e-12*max(1.,sums['energy_before'])
    for key in ('ledger_residual','predictor_ledger_residual','projection_ledger_residual',
                'projection_orthogonality_residual','projection_norm_identity_residual'):
        if abs(metrics[key])>tolerance:
            raise RuntimeError(f'Energy identity failed: {key}={metrics[key]}')
    if predictor_error>2e-12*max(1.,np.sqrt(2*sums['energy_before']/box**3)):
        raise RuntimeError('Independent stencil differs from actual predictor')
    if projection_change>tolerance:
        raise RuntimeError('Orthogonal projection increased energy')
    if return_velocity:
        return metrics,np.asarray(after)
    return metrics


def sampled_accounting(rows, prior_budget):
    """Estimate contributions over time; do not label sparse probes a full ledger.

    Summed one-step work uses left-endpoint time samples at EVERY solver step.
    The saved checkpoints are sparse. Integrating rates below is an estimate.
    The leading left-sum correction of a smooth viscous rate is included and
    labelled as such; it is not an exact quadrature-error correction.
    """
    results=[]
    for n,dt in ((256,.001),(384,.001),(384,.0005)):
        group=[r for r in rows if (r['N'],r['dt'])==(n,dt)]
        t=np.array([r['time'] for r in group])
        transport=np.array([r['transport_work_rate'] for r in group])
        quadratic=np.array([r['quadratic_euler_change']/dt for r in group])
        pressure=np.array([r['projection_change']/dt for r in group])
        dissipation=-np.array([r['viscous_work_rate'] for r in group])
        leading=.5*dt*(dissipation[-1]-dissipation[0])
        old=next(r for r in prior_budget['runs'] if (r['N'],r['dt'])==(n,dt))['summaries']['all_saved']
        for method in ('simpson','trapezoid'):
            integrate=(lambda y:float(simpson(y,x=t))) if method=='simpson' else (lambda y:float(trapezoid(y,t)))
            a,q,p=map(integrate,(transport,quadratic,pressure))
            predicted=a+q+p+leading
            prior=old['residual_'+method]
            results.append(dict(N=n,dt=dt,method=method,transport=a,quadratic_euler=q,
                projection=p,quadratic_plus_projection=q+p,
                leading_viscous_time_sampling_correction=float(leading),
                estimated_gap=predicted,previous_measured_gap=prior,
                remaining_difference=prior-predicted,
                remaining_percent_observed_loss=100*(prior-predicted)/old['observed_energy_loss']))
    return results


def study(input_path,budget_path,cache,probe_cache,out):
    input_path,budget_path,cache,probe_cache,out=map(Path,(input_path,budget_path,cache,probe_cache,out))
    source=json.loads(input_path.read_text())
    budget=json.loads(budget_path.read_text())
    fingerprint,analysis_hash=numerical_source_hash(),file_hash(__file__)
    if (source['numerical_source_sha256']!=fingerprint or budget['numerical_source_sha256']!=fingerprint
            or budget['source_sha256']!=file_hash(input_path) or source['external_force']!=0
            or budget['external_force']!=0 or source['end_time']!=.4):
        raise ValueError('Require the existing unforced saved-flow studies through 0.40')
    flow_path=input_path.parent/source['flow_file']
    if file_hash(flow_path)!=source['flow_sha256']:
        raise ValueError('Flow record changed')
    probe_cache.mkdir(parents=True,exist_ok=True)
    verified,rows=[],[]
    for reference,origin in zip(source['rows'],source['row_sources']):
        n,dt,t=reference['N'],reference['dt'],reference['time']
        path=cache/origin['file']
        if path.name not in verified:
            if file_hash(path)!=origin['summary_sha256'] or file_hash(path.with_suffix('.npz'))!=origin['arrays_sha256']:
                raise ValueError('Checkpoint changed: '+path.name)
            saved=json.loads(path.read_text())
            if (saved['N'],saved['box'],saved['nu'],saved['sigma'],saved['P_U'],saved['source_sha256'])!=(n,6.,.01,0.,0.,fingerprint):
                raise ValueError('Checkpoint settings changed')
            verified.append(path.name)
        provenance=dict(checkpoint=origin,probe_dt=dt,box=6.,nu=.01)
        result_path=probe_cache/f'n{n}-dt{dt:g}-t{t:g}.json'
        old=json.loads(result_path.read_text()) if result_path.exists() else {}
        began=time.perf_counter()
        if old.get('analysis_source_sha256')==analysis_hash and old.get('source')==provenance:
            metrics=old['metrics']
        else:
            with np.load(path.with_suffix('.npz'),allow_pickle=False) as arrays:
                velocity=arrays[origin['array_key']]
            metrics=step_energy(velocity,dt)
            del velocity
            result_path.write_text(json.dumps(dict(source=provenance,analysis_source_sha256=analysis_hash,
                metrics=metrics),allow_nan=False)+'\n')
        if not np.isclose(metrics['energy_before'],reference['energy'],rtol=2e-12,atol=1e-13):
            raise ValueError('Probe starting energy differs from saved measurement')
        expected=.01*6**3*reference['neighbor_gradient_rms']**2
        if not np.isclose(-metrics['viscous_work_rate'],expected,rtol=2e-12,atol=1e-13):
            raise ValueError('Viscous work disagrees with prior gradient measurement')
        if metrics['initial_max_div']>1e-10:
            raise ValueError('Saved field has lost its projected divergence condition')
        row=dict(N=n,dt=dt,time=t,probe_end_time=t+dt,**metrics)
        rows.append(row)
        print(f'N={n} dt={dt:g} t={t:g}: transport={metrics["transport_change"]:.3e} '
              f'viscous={metrics["viscous_change"]:.3e} quadratic={metrics["quadratic_euler_change"]:.3e} '
              f'projection={metrics["projection_change"]:.3e} closure={metrics["ledger_residual"]:.2e} '
              f'({time.perf_counter()-began:.1f}s)',flush=True)
    result=dict(start_time=0.,end_time=.4,external_force=0.,saved_trajectory_extension_steps=0,
        independent_one_step_probes=len(rows),numerical_source_sha256=fingerprint,analysis_source_sha256=analysis_hash,
        input_file=input_path.name,input_sha256=file_hash(input_path),budget_file=budget_path.name,budget_sha256=file_hash(budget_path),
        source_records=source['sources'],row_sources=source['row_sources'],rows=rows,
        sampled_accounting=sampled_accounting(rows,budget),
        limits=['Each output is one step from a saved state on disposable arrays; it is not the next saved checkpoint or a new continued trajectory.',
            'The four contributions form an exact algebraic identity for each computed step, up to floating-point roundoff. Transport and viscosity are simultaneous parts of one predictor; the quadratic term includes their cross terms.',
            'Sparse probes do not reconstruct every intervening step. Integrated contributions and the leading viscous left-sampling correction assume time smoothness and remain estimates.',
            'The positive quadratic Euler term is an accounting term of the discrete update, not an external force or a proof of a continuous energy source.',
            'No numerical-scheme, equation, starting-field, viscosity, grid, or force changes were made. This is not a Clay proof.'])
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with out.with_suffix('.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=Path('math/results/hug-resolution.json'))
    p.add_argument('--budget',type=Path,default=Path('math/results/hug-energy-balance.json'))
    p.add_argument('--cache',type=Path,default=Path('scratch/hug-refinement'))
    p.add_argument('--probe-cache',type=Path,default=Path('scratch/hug-step-energy'))
    p.add_argument('--out',type=Path,default=Path('math/results/hug-step-energy.json'))
    a=p.parse_args();study(a.input,a.budget,a.cache,a.probe_cache,a.out)
