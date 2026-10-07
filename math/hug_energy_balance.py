"""Estimate the hug flow's viscous energy loss from published measurements.

No checkpoint arrays are evolved or changed. The adjacent-cell gradient norm
matches the existing nearest-neighbor Laplacian's discrete energy identity.
Time integrals use saved observations, so their differences are sensitivity
measurements, not rigorous error bounds or an exact per-step energy ledger.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson, trapezoid

from continue_hug_refinement import file_hash, numerical_source_hash


def dissipation_rate(gradient_rms, box, nu):
    """nu * integral |D+ u|^2 = nu * box^3 * RMS(D+ u)^2."""
    g = np.asarray(gradient_rms, dtype=float)
    if (not np.isfinite(box) or box <= 0 or not np.isfinite(nu) or nu < 0
            or not np.isfinite(g).all() or np.any(g < 0)):
        raise ValueError('Require finite nonnegative gradients/viscosity and a positive box')
    return nu * box**3 * g**2


def energy_budget(times, energies, gradient_rms, box=6., nu=.01):
    """Compare observed loss with two integrals over the same saved times.

    Composite Simpson uses an odd number of observations. Also report a
    trapezoidal integral after retaining every second observation, including
    both endpoints. A positive residual means LESS energy was lost than the
    viscous estimate predicts. No numerical-error threshold is imposed.
    """
    t, e, g = [np.asarray(a, dtype=float) for a in (times, energies, gradient_rms)]
    if (t.ndim != 1 or e.shape != t.shape or g.shape != t.shape
            or len(t) < 3 or len(t) % 2 != 1
            or not np.isfinite(t).all() or not np.isfinite(e).all()
            or np.any(np.diff(t) <= 0) or np.any(e < 0)):
        raise ValueError('Require odd, matching finite series with increasing times and nonnegative energy')
    q = dissipation_rate(g, box, nu)
    accumulated = cumulative_trapezoid(q, t, initial=0.)
    observed = e[0] - e
    trap, quadratic = float(accumulated[-1]), float(simpson(q, x=t))
    coarse = float(trapezoid(q[::2], t[::2]))
    loss = float(observed[-1])
    percentage = lambda x: float(100*x/abs(loss)) if loss != 0 else None
    rows = [dict(time=float(tt), energy=float(ee), neighbor_gradient_rms=float(gg),
                 viscous_loss_rate=float(qq), observed_energy_loss=float(ll),
                 cumulative_viscous_trapezoid=float(ii),
                 residual_trapezoid=float(ii-ll))
            for tt, ee, gg, qq, ll, ii in zip(t, e, g, q, observed, accumulated)]
    summary = dict(start_time=float(t[0]), end_time=float(t[-1]),
        observations=len(t), max_observation_gap=float(np.diff(t).max()),
        initial_energy=float(e[0]), final_energy=float(e[-1]), observed_energy_loss=loss,
        viscous_trapezoid=trap, viscous_simpson=quadratic, viscous_coarsened_trapezoid=coarse,
        residual_trapezoid=trap-loss, residual_simpson=quadratic-loss,
        trapezoid_gap_percent_loss=percentage(trap-loss), simpson_gap_percent_loss=percentage(quadratic-loss),
        quadrature_method_difference=abs(trap-quadratic),
        quadrature_method_difference_percent_loss=percentage(abs(trap-quadratic)),
        coarsening_difference=abs(coarse-trap),
        coarsening_difference_percent_loss=percentage(abs(coarse-trap)),
        energy_decreases_at_saved_times=bool(np.all(np.diff(e) <= 0)),
        loss_rate_decreases_at_saved_times=bool(np.all(np.diff(q) <= 0)))
    return dict(rows=rows, summary=summary)


def study(input_path, out):
    input_path, out = Path(input_path), Path(out)
    data = json.loads(input_path.read_text())
    folder = Path(__file__).parent
    history_path, flow_path = input_path.parent/data['history_file'], input_path.parent/data['flow_file']
    flow = json.loads(flow_path.read_text())
    history = json.loads(history_path.read_text())
    source_hash = numerical_source_hash()
    if (data['numerical_source_sha256'] != source_hash
            or data['analysis_source_sha256'] != file_hash(folder/'hug_resolution.py')
            or data['flow_sha256'] != file_hash(flow_path)
            or data['history_sha256'] != file_hash(history_path)
            or data['external_force'] != 0 or data['new_evolution_steps'] != 0
            or data['start_time'] != 0 or data['end_time'] != .4
            or flow['external_force'] != 0 or flow['end_time'] != .4):
        raise ValueError('Source provenance or unforced interval changed')
    identities = [(r['N'], r['dt'], r['time']) for r in data['rows']]
    if len(set(identities)) != len(identities) or len(identities) != len(history['rows']):
        raise ValueError('Duplicate or missing observations')
    for row, prior in zip(data['rows'], history['rows']):
        if (tuple(row[k] for k in ('N','dt','time')) != tuple(prior[k] for k in ('N','dt','time'))
                or not np.isclose(row['energy'], prior['energy'], rtol=2e-12, atol=1e-13)):
            raise ValueError('Energy or observation identity differs from prior measurement')
    controls = [(256,.001),(384,.001),(384,.0005)]
    groups = [[r for r in data['rows'] if (r['N'],r['dt']) == c] for c in controls]
    common_times = sorted(set.intersection(*[set(r['time'] for r in group) for group in groups]))
    runs, rows = [], []
    for (n, dt), group in zip(controls, groups):
        config = next(r for r in flow['runs'] if (r['N'],r['dt']) == (n,dt))
        if (config['box'],config['nu'],config['sigma'],config['P_U'],config['source_sha256']) != (6.,.01,0.,0.,source_hash):
            raise ValueError('Numerical settings differ from recorded study')
        summaries = {}
        schedules = dict(all_saved=group,
            common_saved=[r for r in group if r['time'] in common_times],
            uniform_later=[r for r in group if r['time'] >= .16])
        for schedule, samples in schedules.items():
            budget = energy_budget([r['time'] for r in samples], [r['energy'] for r in samples],
                [r['neighbor_gradient_rms'] for r in samples], config['box'], config['nu'])
            summaries[schedule] = budget['summary']
            rows.extend(dict(N=n,dt=dt,schedule=schedule,**r) for r in budget['rows'])
        runs.append(dict(N=n,dt=dt,box=config['box'],nu=config['nu'],summaries=summaries))
    result = dict(start_time=0.,end_time=.4,external_force=0.,new_evolution_steps=0,
        observations_reused=len(data['rows']), numerical_source_sha256=source_hash,
        analysis_source_sha256=file_hash(__file__),
        source_file=input_path.name,source_sha256=file_hash(input_path),
        source_analysis_sha256=data['analysis_source_sha256'],
        history_file=history_path.name,history_sha256=file_hash(history_path),
        flow_file=flow_path.name,flow_sha256=file_hash(flow_path),
        checkpoint_provenance_reused=data['sources'],
        common_saved_times=common_times,runs=runs,rows=rows,
        definitions=dict(
            continuum='For smooth unforced periodic flow: E(t)-E(0) = -nu * integral_0^t integral |grad u|^2 dx ds.',
            discrete='For the current nearest-neighbor Laplacian: -nu <u,L_h u>_h = nu * box^3 * RMS(D+ u)^2. This is its viscous contribution, not an identity for the entire Euler/projection update.',
            residual='Estimated viscous loss minus observed loss; positive means the saved flow retained more energy than viscosity alone predicts.',
            percent='100 * residual / abs(observed energy loss); undefined for zero observed loss.',
            trapezoid='Piecewise linear time integration of the saved dissipation rates.',
            simpson='Composite integration of quadratic interpolants through successive triples of saved rates.',
            coarsened='Trapezoidal integration after dropping every other observation, retaining both endpoints.',
            common='All three controls use the same saved observation times for the grid comparison.',
            later='The interval 0.16 to 0.40 has uniform observation spacing 0.04 for every control.'),
        limits=['No new flow simulation or remeasurement of checkpoint arrays; previously verified summaries and hashes are reused.',
            'Differences between quadrature estimates and coarsened samples are sensitivity indicators, not rigorous error bars; unsampled time behavior is unknown.',
            'The residual combines discretization and time-integration error; it does not isolate advection, pressure-projection, or forward-Euler contributions.',
            'Adjacent-cell gradients match the discrete viscous operator but are not exact continuum gradients.',
            'A smaller time-step difference is evidence of sensitivity, not a proof of convergence order or global smoothness.'])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with out.with_suffix('.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('math/results/hug-resolution.json'))
    parser.add_argument('--out',type=Path,default=Path('math/results/hug-energy-balance.json'))
    args=parser.parse_args()
    result=study(args.input,args.out)
    for run in result['runs']:
        print(json.dumps(dict(N=run['N'],dt=run['dt'],**run['summaries']['all_saved'])))
