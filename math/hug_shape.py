"""Measure the spatial spread of kinetic energy in saved hug velocities.

Read-only postprocessing: no solver steps and no imposed breathing schedule.
Widths are central second moments in the fixed periodic-box coordinate chart;
they are not material boundaries, particle displacements, or pressure measures.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from continue_hug_refinement import file_hash, numerical_source_hash


def energy_shape(velocity, box=6.0):
    """Return energy-weighted RMS widths and energy near the box faces.

    Sum in slabs to avoid additional full-grid temporary arrays. Uniform cell
    volume and the factor 1/2 cancel in normalized spatial moments.
    """
    velocity = np.asarray(velocity)
    if velocity.ndim != 4 or velocity.shape[0] != 3:
        raise ValueError('Expected three velocity components on a cubic grid')
    n = velocity.shape[1]
    if n < 2 or velocity.shape[1:] != (n, n, n):
        raise ValueError('Expected a cubic grid with at least two cells per axis')
    if not np.isfinite(box) or box <= 0:
        raise ValueError('Box length must be finite and positive')
    coords = (np.arange(n) - (n-1)/2)*box/n
    marginals = np.zeros((3, n), dtype=float)
    # Union of the outer 10% of each half-width: |x|, |y| or |z| >= .45 L.
    edge = np.abs(coords) >= .45*box
    yz_edge = edge[:, None] | edge[None, :]
    edge_weight = 0.
    for first in range(0, n, 8):
        slab = velocity[:, first:first+8].astype(float, copy=False)
        if not np.isfinite(slab).all():
            raise ValueError('Velocity must be finite')
        weight = np.sum(slab*slab, axis=0)
        marginals[0, first:first+8] = weight.sum(axis=(1, 2))
        marginals[1] += weight.sum(axis=(0, 2))
        marginals[2] += weight.sum(axis=(0, 1))
        for offset, plane in enumerate(weight):
            edge_weight += plane.sum() if edge[first+offset] else plane[yz_edge].sum()
    total = marginals[0].sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError('A finite, nonzero kinetic energy is required')
    centers = marginals @ coords / total
    variances = np.array([np.dot(m, (coords-c)**2)/total
                          for m, c in zip(marginals, centers)])
    return dict(energy=float(.5*(box/n)**3*total),
                center_x=float(centers[0]), center_y=float(centers[1]),
                center_z=float(centers[2]),
                radial_width=float(np.sqrt(variances[0]+variances[1])),
                axial_width=float(np.sqrt(variances[2])),
                rms_radius=float(np.sqrt(variances.sum())),
                edge_energy_fraction=float(edge_weight/total))


def summarize_shapes(rows, configurations, end):
    summaries, latest = [], []
    for n, dt in configurations:
        current = sorted((r for r in rows if (r['N'],r['dt'])==(n,dt)),
                         key=lambda r:r['time'])
        start, final = current[0], current[-1]
        if start['time'] != 0 or final['time'] != end:
            raise ValueError('Shape history must span zero to the requested endpoint')
        summary = dict(N=n, dt=dt, sample_times=[r['time'] for r in current])
        for name in ('radial_width', 'axial_width', 'rms_radius'):
            differences = np.diff([r[name] for r in current])
            summary[name+'_change_percent'] = 100*(final[name]/start[name]-1)
            summary[name+'_sample_trend'] = ('increasing' if np.all(differences > 0)
                else 'decreasing' if np.all(differences < 0) else 'mixed or constant')
        summaries.append(summary)
        latest.append(final)
    comparisons = {}
    for label, a, b in [('grid', latest[0], latest[1]), ('time_step', latest[1], latest[2])]:
        comparisons[label] = {name+'_relative_difference_percent': 100*abs(a[name]-b[name])/b[name]
                              for name in ('radial_width', 'axial_width', 'rms_radius')}
    return summaries, comparisons


def save_study(result, out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with out.with_suffix('.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result['rows'][0]))
        writer.writeheader()
        writer.writerows(result['rows'])


def extend_study(prior_path, cache, stops, out):
    """Reuse recorded shape observations and read only new checkpoint arrays."""
    prior_path, cache, out = Path(prior_path), Path(cache), Path(out)
    previous = json.loads(prior_path.read_text())
    if prior_path.resolve() == out.resolve():
        raise ValueError('Preserve the prior shape record; use a new output path')
    source = numerical_source_hash()
    configurations = ((256,.001),(384,.001),(384,.0005))
    start = previous['end_time']
    if previous['numerical_source_sha256'] != source or previous['external_force'] != 0:
        raise ValueError('Require the unchanged unforced numerical source')
    if (not stops or any(not math.isfinite(t) for t in stops)
            or any(b <= a for a,b in zip((start,)+tuple(stops[:-1]), stops))
            or any(abs(round(t/dt)*dt-t) > 1e-12 for _,dt in configurations for t in stops)):
        raise ValueError('Shape checkpoint times must increase after the saved endpoint and align with dt')
    # Validate the three retained histories before reading any new array.
    summarize_shapes(previous['rows'], configurations, start)
    rows, sources = list(previous['rows']), list(previous['sources'])
    for n, dt in configurations:
        for stop in stops:
            path = cache/f'n{n}-dt{dt:g}-t{stop:g}.json'
            saved = json.loads(path.read_text())
            if (saved['N'], saved['dt'], saved['end_time'], saved['box'], saved['nu'],
                saved['sigma'], saved['P_U'], saved['source_sha256']) != (
                    n, dt, stop, 6., .01, 0., 0., source):
                raise ValueError(f'Unexpected checkpoint settings: {path}')
            array_path = path.with_suffix('.npz')
            with np.load(array_path) as arrays:
                velocity = arrays['final']
                if velocity.shape != (3,n,n,n):
                    raise ValueError('Checkpoint grid differs from its metadata')
                shape = energy_shape(velocity, saved['box'])
                del velocity
            if not np.isclose(shape['energy'], saved['rows'][-1]['energy'], rtol=2e-12, atol=1e-13):
                raise ValueError('Checkpoint energy differs from saved diagnostic')
            rows.append(dict(N=n,dt=dt,time=stop,**shape))
            sources.append(dict(file=path.name,summary_sha256=file_hash(path),
                                arrays_sha256=file_hash(array_path)))
            print(f'Shape N={n} dt={dt:g} t={stop:g}: radial={shape["radial_width"]:.8f} '
                  f'axial={shape["axial_width"]:.8f}',flush=True)
    rows.sort(key=lambda r:(r['N'],-r['dt'],r['time']))
    summaries, comparisons = summarize_shapes(rows,configurations,stops[-1])
    result = dict(previous, end_time=stops[-1], rows=rows, sources=sources,
        summary=summaries,final_comparisons=comparisons,analysis_source_sha256=file_hash(__file__),
        continuation=dict(previous_file=prior_path.name,previous_sha256=file_hash(prior_path),
            previous_analysis_source_sha256=previous['analysis_source_sha256'],
            reused_rows=len(previous['rows']),new_times=list(stops),new_measurements=3*len(stops)))
    save_study(result,out)
    print(json.dumps(dict(summary=summaries,final_comparisons=comparisons),indent=2),flush=True)
    return result


def study(cache, out):
    cache, out = Path(cache), Path(out)
    fingerprint = numerical_source_hash()
    rows, sources = [], []
    configurations = ((256, .001), (384, .001), (384, .0005))
    for n, dt in configurations:
        stops = (.08, .16, .20, .24) if n == 256 else (.04, .08, .12, .16, .20, .24)
        current = []
        for index, stop in enumerate(stops):
            path = cache/f'n{n}-dt{dt:g}-t{stop:g}.json'
            saved = json.loads(path.read_text())
            if (saved['N'], saved['dt'], saved['end_time'], saved['box'],
                saved['nu'], saved['sigma'], saved['P_U'], saved['source_sha256']) != (
                    n, dt, stop, 6., .01, 0., 0., fingerprint):
                raise ValueError(f'Unexpected checkpoint settings: {path}')
            arrays_path = path.with_suffix('.npz')
            sources.append(dict(file=path.name, summary_sha256=file_hash(path),
                                arrays_sha256=file_hash(arrays_path)))
            with np.load(arrays_path) as arrays:
                for key, t, reference in ([('initial', 0., saved['rows'][0])]
                        if index == 0 else []) + [('final', stop, saved['rows'][-1])]:
                    velocity = arrays[key]
                    if velocity.shape != (3, n, n, n):
                        raise ValueError('Checkpoint grid differs from its metadata')
                    shape = energy_shape(velocity, saved['box'])
                    del velocity
                    if not np.isclose(shape['energy'], reference['energy'], rtol=2e-12, atol=1e-13):
                        raise ValueError('Checkpoint energy differs from saved diagnostic')
                    row = dict(N=n, dt=dt, time=t, **shape)
                    current.append(row)
                    print(f'N={n} dt={dt:g} t={t:g}: radial={shape["radial_width"]:.8f} '
                          f'axial={shape["axial_width"]:.8f}', flush=True)
        rows.extend(current)
    summaries, comparisons = summarize_shapes(rows,configurations,.24)
    result = dict(start_time=0., end_time=.24, external_force=0., new_evolution_steps=0,
        numerical_source_sha256=fingerprint, analysis_source_sha256=file_hash(__file__),
        method='Kinetic-energy-weighted central second moments of saved cell-centered velocity',
        definitions=dict(weight='u_x^2 + u_y^2 + u_z^2',
            center='sum(weight * coordinate) / sum(weight)',
            radial_width='sqrt(Var_E(x) + Var_E(y))', axial_width='sqrt(Var_E(z))',
            rms_radius='sqrt(Var_E(x) + Var_E(y) + Var_E(z))',
            edge_region='max(abs(x), abs(y), abs(z)) >= 0.45 * box_length',
            coordinates='cell centers in [-3, 3); fixed coordinate chart on the periodic box'),
        limits=['Energy widths describe the distribution of motion, not material boundaries or particle paths.',
            'Local dissipation, transport and redistribution can change these widths; the diagnostic does not isolate a cause.',
            'A prescribed hug breathing cycle is not applied during fluid evolution.',
            'Saved times cannot resolve changes between snapshots or establish a repeating cycle.',
            'Euclidean moments depend on the chosen periodic coordinate cut; edge energy is reported as a localization check, not an error bound.',
            'Finite-time numerical evidence has no rigorous exact-solution error bound.'],
        sources=sources, rows=rows, summary=summaries, final_comparisons=comparisons)
    save_study(result,out)
    print(json.dumps(dict(summary=summaries, final_comparisons=comparisons), indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('scratch/hug-refinement'))
    parser.add_argument('--out', type=Path, default=Path('math/results/hug-shape.json'))
    parser.add_argument('--from-study', type=Path, help='Reuse earlier shape observations')
    parser.add_argument('--stops', type=float, nargs='+', help='Only the new checkpoint times to measure')
    args = parser.parse_args()
    if (args.from_study is None) != (args.stops is None):
        parser.error('--from-study and --stops must be supplied together')
    if args.from_study is not None:
        extend_study(args.from_study,args.cache,args.stops,args.out)
    else:
        study(args.cache, args.out)
