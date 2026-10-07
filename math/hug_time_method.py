"""Local Euler/Heun controls from saved hug states; no trajectory replacement.

Heun is composed from two existing projected Euler maps. Spatial operators,
viscosity, pressure projection, initial arrays and zero force are preserved.
Each control covers the SAME short interval using 1, 2 or 4 substeps.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from continue_hug_refinement import file_hash, numerical_source_hash
from navier import PurePythonNavierStokes3D


def energy(u, box):
    return .5*(box/u.shape[1])**3*sum(float(np.sum(c*c)) for c in u)


def difference(a, b, box):
    """Complete-field L2 and maximum difference, without grid interpolation."""
    norm2 = maximum = 0.
    for c in range(3):
        for first in range(0, a.shape[1], 8):
            d = a[c, first:first+8] - b[c, first:first+8]
            norm2 += float(np.sum(d*d))
            maximum = max(maximum, float(np.max(np.abs(d))))
    l2 = float(np.sqrt((box/a.shape[1])**3*norm2))
    return dict(velocity_l2_difference=l2, velocity_max_difference=maximum,
                relative_velocity_l2_difference=l2/np.sqrt(2*energy(b, box)),
                energy_difference=energy(a, box)-energy(b, box))


def euler_map(u, dt, box, nu):
    """Run the actual original predictor and pressure correction."""
    n = u.shape[1]
    sim = PurePythonNavierStokes3D(N=4)
    sim.N, sim.dx, sim.nu, sim.sigma = n, box/n, nu, 0.
    sim.u, sim.v, sim.w = u
    sim.S = np.zeros((n, n, n), dtype=float)
    sim.rho_epsilon = 0.
    sim.step(dt, 0., backend='numpy')
    diagnostic = sim.project(backend='fft')
    if np.any(sim.S != 0):
        raise RuntimeError('Zero-force scalar control changed')
    return np.asarray((sim.u, sim.v, sim.w)), diagnostic['after']


def method_step(u, dt, method, box=6., nu=.01):
    """One step and its energy-work identity, for a projected input field.

Euler: v=u+h*f. Heun: y=u+h*f, v=u+h*(f+g)/2, g=G(y).
G is the unchanged projected spatial right-hand side. With Q the method's
stage quadrature of <u,G(u)>_h, E(v)-E(u)-Q is h^2||f||^2/2 for Euler
and h^2||g-f||^2/8 for Heun. These are discrete identities, not exact
continuum energy errors or substitutes for the refinement comparison.
"""
    if method not in ('euler', 'heun'):
        raise ValueError('method must be euler or heun')
    y, max_div = euler_map(u, dt, box, nu)
    if method == 'heun':
        z, div2 = euler_map(y, dt, box, nu)
        max_div = max(max_div, div2)
        v = .5*u + .5*z
    else:
        z, v = None, y
    weight = (box/u.shape[1])**3
    work = defect = stable_change = 0.
    for c in range(3):
        for first in range(0, u.shape[1], 8):
            s = slice(first, first+8)
            before, predictor, after = u[c,s], y[c,s], v[c,s]
            d0 = predictor-before
            if method == 'euler':
                work += weight*float(np.sum(before*d0))
                defect += .5*weight*float(np.sum(d0*d0))
            else:
                d1 = z[c,s]-predictor
                work += .5*weight*float(np.sum(before*d0+predictor*d1))
                defect += .125*weight*float(np.sum((d1-d0)**2))
            stable_change += .5*weight*float(np.sum((after-before)*(after+before)))
    residual = stable_change-work-defect
    if not np.isfinite(stable_change) or abs(residual) > 5e-12*max(1.,energy(u,box)):
        raise RuntimeError('Nonfinite field or failed step energy identity')
    return v, dict(work_quadrature=work, temporal_energy_defect=defect,
                   stable_energy_change=stable_change, identity_residual=residual,
                   max_stage_divergence=max_div)


def evolve(u, interval, substeps, method, box=6., nu=.01, progress=None):
    u = np.asarray(u)
    if (u.ndim != 4 or u.shape[0] != 3 or u.shape[1] < 4
            or u.shape[1:] != (u.shape[1],)*3 or not np.isrealobj(u)
            or not np.isfinite(interval) or interval <= 0
            or not isinstance(substeps, int) or substeps < 1
            or not np.isfinite(box) or box <= 0 or not np.isfinite(nu) or nu < 0
            or method not in ('euler','heun')):
        raise ValueError('Invalid field, method or integration settings')
    for c in u:
        for first in range(0,u.shape[1],8):
            if not np.isfinite(c[first:first+8]).all():
                raise ValueError('Velocity must be finite')
    dx = box/u.shape[1]
    div = sum((np.roll(u[d],-1,d)-np.roll(u[d],1,d))/(2*dx) for d in range(3))
    if np.max(np.abs(div)) > 1e-10:
        raise ValueError('Control requires an already projected input')
    del div
    before = energy(u,box)
    current = u
    metrics = dict(work_quadrature=0.,temporal_energy_defect=0.,stable_energy_change=0.,
                   max_identity_residual=0.,max_stage_divergence=0.)
    for step in range(substeps):
        current, ledger = method_step(current,interval/substeps,method,box,nu)
        for key in ('work_quadrature','temporal_energy_defect','stable_energy_change'):
            metrics[key] += ledger[key]
        metrics['max_identity_residual'] = max(metrics['max_identity_residual'],abs(ledger['identity_residual']))
        metrics['max_stage_divergence'] = max(metrics['max_stage_divergence'],ledger['max_stage_divergence'])
        if progress:
            progress(step+1,substeps)
    metrics.update(energy_before=before,energy_after=energy(current,box))
    metrics['energy_change'] = metrics['energy_after']-before
    metrics['total_identity_residual'] = metrics['energy_change']-metrics['work_quadrature']-metrics['temporal_energy_defect']
    if (not all(np.isfinite(x) for x in metrics.values())
            or abs(metrics['total_identity_residual']) > 5e-12*max(1.,before)):
        raise RuntimeError('Nonfinite diagnostic or accumulated identity failure')
    return current,metrics


def compare_methods(u, interval=.001, box=6., nu=.01, progress=None):
    rows, comparisons = [], []
    for method in ('euler','heun'):
        previous = None
        for count in (1,2,4):
            notify = (lambda i,n:progress(method,count,i,n)) if progress else None
            endpoint,metrics = evolve(u,interval,count,method,box,nu,notify)
            rows.append(dict(method=method,substeps=count,dt=interval/count,**metrics))
            if previous is not None:
                comparisons.append(dict(method=method,coarse_substeps=count//2,
                    fine_substeps=count,**difference(previous,endpoint,box)))
            previous = endpoint
        del previous,endpoint
    return dict(rows=rows,comparisons=comparisons)


def study(input_path,cache,out,interval=.001):
    input_path,cache,out = map(Path,(input_path,cache,out))
    source = json.loads(input_path.read_text())
    fingerprint,analysis = numerical_source_hash(),file_hash(__file__)
    if (source['numerical_source_sha256'] != fingerprint or source['external_force'] != 0
            or source['end_time'] != .4):
        raise ValueError('Require the verified unforced saved-flow study through 0.40')
    # Same-history full-step inputs on both grids and the finest half-step input.
    # All method/substep controls WITHIN a run start from identical arrays.
    configs = ((256,.001),(384,.001),(384,.0005))
    results = []
    probe_cache = cache.parent/'hug-time-method'
    probe_cache.mkdir(exist_ok=True,parents=True)
    for n,history_dt in configs:
        origin = next(r for r in source['row_sources'] if (r['N'],r['dt'],r['time'])==(n,history_dt,.4))
        reference = next(r for r in source['rows'] if (r['N'],r['dt'],r['time'])==(n,history_dt,.4))
        path = cache/origin['file']
        if file_hash(path)!=origin['summary_sha256'] or file_hash(path.with_suffix('.npz'))!=origin['arrays_sha256']:
            raise ValueError('Checkpoint hash changed')
        saved = json.loads(path.read_text())
        if (saved['N'],saved['box'],saved['nu'],saved['sigma'],saved['P_U'],saved['source_sha256'])!=(n,6.,.01,0.,0.,fingerprint):
            raise ValueError('Checkpoint equation settings changed')
        provenance = dict(checkpoint=origin,interval=interval,box=6.,nu=.01,
                          input_sha256=file_hash(input_path),numerical_source_sha256=fingerprint,
                          analysis_source_sha256=analysis)
        result_path = probe_cache/f'n{n}-historydt{history_dt:g}-interval{interval:g}.json'
        previous = json.loads(result_path.read_text()) if result_path.exists() else {}
        if previous.get('provenance') == provenance:
            result = previous['result']
            print(f'Reuse completed N={n} history dt={history_dt:g}',flush=True)
        else:
            with np.load(path.with_suffix('.npz'),allow_pickle=False) as arrays:
                u = arrays[origin['array_key']]
            if not np.isclose(energy(u,6.),reference['energy'],rtol=2e-12,atol=1e-13):
                raise ValueError('Starting energy disagrees with published measurement')
            began = time.perf_counter()
            def progress(method,count,i,total):
                print(f'N={n} historydt={history_dt:g} {method} {count} substeps: '
                      f'{i}/{total} ({time.perf_counter()-began:.1f}s)',flush=True)
            result = compare_methods(u,interval,progress=progress)
            del u
            result_path.write_text(json.dumps(dict(provenance=provenance,result=result),allow_nan=False)+'\n')
        results.append(dict(N=n,history_dt=history_dt,start_time=.4,end_time=.4+interval,
                            provenance=provenance,**result))
    output = dict(question='Does projected Heun reduce local time-step sensitivity and temporal energy defect?',
        saved_trajectory_end=.4,saved_trajectory_extension_steps=0,external_force=0.,box=6.,nu=.01,
        interval=interval,numerical_source_sha256=fingerprint,analysis_source_sha256=analysis,
        input_file=input_path.name,input_sha256=file_hash(input_path),runs=results,
        limits=['These are disposable short-interval controls at the final saved state, not a replacement of the Euler trajectory through 0.40.',
            'Within each control the physical equation and spatial operators are identical; only the time integrator and substep size vary.',
            'Temporal energy defect is relative to each method\'s stage quadrature of the projected spatial energy rate, including both viscosity and transport. It is not the earlier viscosity-only interval gap.',
            'Endpoint refinement differences measure time-step sensitivity, not errors against an exact solution or rigorous bounds.',
            'Spatial discretization, initial-field sampling and earlier Euler history errors remain. This supplies no Clay proof.'])
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    flat = [dict(N=r['N'],history_dt=r['history_dt'],start_time=r['start_time'],end_time=r['end_time'],**row)
            for r in results for row in r['rows']]
    with out.with_suffix('.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    print('Completed all method controls; original trajectory unchanged.',flush=True)
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('math/results/hug-resolution.json'))
    parser.add_argument('--cache',type=Path,default=Path('scratch/hug-refinement'))
    parser.add_argument('--out',type=Path,default=Path('math/results/hug-time-method.json'))
    parser.add_argument('--interval',type=float,default=.001)
    args=parser.parse_args();study(args.input,args.cache,args.out,args.interval)
