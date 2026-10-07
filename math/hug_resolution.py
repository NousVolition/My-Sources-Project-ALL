"""Inspect saved hug velocities for centered-grid blind modes and short waves.

Read-only analysis. No filtering, pressure correction, or solver steps are applied.
The FFT is used to measure discrete energy and two discrete gradient norms.
Neither gradient is presented as the exact continuum gradient.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.fft import rfftn

from continue_hug_refinement import file_hash, numerical_source_hash


def grid_spectrum(velocity, box=6.0):
    """Measure a real, even, cubic grid, using the full energy of a real FFT.

    The reduced last axis is counted twice except at DC and Nyquist. Symbols
    sin(kh)/h and 2*sin(kh/2)/h give centered and adjacent-cell gradient norms.
    They describe the SAME saved arrays; adjacent-cell differences detect every
    nonconstant grid Fourier mode, including the seven joint blind modes.
    """
    velocity = np.asarray(velocity)
    if velocity.ndim != 4 or velocity.shape[0] != 3 or not np.isrealobj(velocity):
        raise ValueError('Expected three real velocity components on a cubic grid')
    n = velocity.shape[1]
    if n < 4 or n % 2 or velocity.shape[1:] != (n, n, n):
        raise ValueError('Expected an even cubic grid with at least four cells per axis')
    if not np.isfinite(box) or box <= 0:
        raise ValueError('Box length must be finite and positive')
    power = np.zeros((n, n, n//2+1), dtype=float)
    direct_sum = 0.
    for component in velocity:
        for first in range(0, n, 8):
            slab = component[first:first+8].astype(float, copy=False)
            if not np.isfinite(slab).all():
                raise ValueError('Velocity must be finite')
            direct_sum += float(np.sum(slab*slab))
        transformed = rfftn(component, norm='backward', workers=1)
        for first in range(0, n, 8):
            slab = transformed[first:first+8]
            power[first:first+8] += slab.real**2 + slab.imag**2
        del transformed
    if not np.isfinite(direct_sum) or direct_sum <= 0:
        raise ValueError('A finite, nonzero kinetic energy is required')
    modes = np.abs(np.rint(np.fft.fftfreq(n)*n)).astype(int)
    z_modes = np.arange(n//2+1)
    multiplicity = np.full(n//2+1, 2.)
    multiplicity[[0, -1]] = 1.
    dx = box/n
    center = (np.sin(2*np.pi*modes/n)/dx)**2
    center[modes == n//2] = 0.
    neighbor = (2*np.sin(np.pi*modes/n)/dx)**2
    center_z = (np.sin(2*np.pi*z_modes/n)/dx)**2
    center_z[-1] = 0.
    neighbor_z = (2*np.sin(np.pi*z_modes/n)/dx)**2
    totals = dict(mass=0., checkerboard=0., nyquist_planes=0., high_band=0.,
                  center=0., neighbor=0., high_band_neighbor=0.)
    shells = np.zeros(n//2+1)
    gradient_shells = np.zeros_like(shells)
    yz_max = np.maximum(modes[:, None], z_modes[None, :])
    yz_blind = ((modes[:, None] == 0) | (modes[:, None] == n//2)) & (
        (z_modes[None, :] == 0) | (z_modes[None, :] == n//2))
    for first in range(0, n, 8):
        last = min(first+8, n)
        x_modes = modes[first:last, None, None]
        shell = np.maximum(x_modes, yz_max[None, :, :])
        weighted = power[first:last]*multiplicity[None, None, :]
        joint_blind = ((x_modes == 0) | (x_modes == n//2)) & yz_blind & (shell > 0)
        high = shell >= n/4
        c_weight = center[first:last, None, None] + center[None, :, None] + center_z
        e_weight = neighbor[first:last, None, None] + neighbor[None, :, None] + neighbor_z
        edge_power = weighted*e_weight
        totals['mass'] += float(weighted.sum())
        totals['checkerboard'] += float(weighted[joint_blind].sum())
        totals['nyquist_planes'] += float(weighted[shell == n//2].sum())
        totals['high_band'] += float(weighted[high].sum())
        totals['center'] += float(np.sum(weighted*c_weight))
        totals['neighbor'] += float(edge_power.sum())
        totals['high_band_neighbor'] += float(edge_power[high].sum())
        shells += np.bincount(shell.ravel(), weights=weighted.ravel(), minlength=n//2+1)
        gradient_shells += np.bincount(shell.ravel(), weights=edge_power.ravel(), minlength=n//2+1)
    mass, edge = totals['mass'], totals['neighbor']
    fft_sum = mass/n**3
    if not np.isfinite(fft_sum) or not np.isclose(fft_sum, direct_sum, rtol=2e-12, atol=1e-13):
        raise ValueError('Real-space and FFT energies disagree')
    metrics = dict(energy=.5*dx**3*direct_sum,
        parseval_relative_difference=abs(fft_sum-direct_sum)/direct_sum,
        checkerboard_energy_fraction=totals['checkerboard']/mass,
        nyquist_plane_energy_fraction=totals['nyquist_planes']/mass,
        short_wave_energy_fraction=totals['high_band']/mass,
        short_wave_neighbor_gradient_fraction=totals['high_band_neighbor']/edge if edge else 0.,
        centered_gradient_rms=float(np.sqrt(totals['center']/n**6)),
        neighbor_gradient_rms=float(np.sqrt(edge/n**6)),
        gradient_rms_gap_percent=100*(1-np.sqrt(totals['center']/edge)) if edge else 0.)
    spectrum = dict(max_axis_mode=list(range(n//2+1)),
                    energy_fraction=(shells/mass).tolist(),
                    neighbor_gradient_fraction=(gradient_shells/edge).tolist() if edge else shells.tolist())
    if not edge:
        spectrum['neighbor_gradient_fraction'] = [0.]*len(shells)
    return metrics, spectrum


def study(history_path, flow_path, cache, analysis_cache, out):
    """Analyze each saved field once, caching results and checking provenance."""
    history_path, flow_path, cache = map(Path, (history_path, flow_path, cache))
    analysis_cache, out = Path(analysis_cache), Path(out)
    history, flow = json.loads(history_path.read_text()), json.loads(flow_path.read_text())
    source, analysis_source = numerical_source_hash(), file_hash(__file__)
    if (history['numerical_source_sha256'] != source or history['end_time'] != .4
            or history['external_force'] != 0 or flow['end_time'] != .4
            or flow['external_force'] != 0 or not flow['initial_arrays_identical']):
        raise ValueError('Require the verified unforced study through 0.40')
    analysis_cache.mkdir(parents=True, exist_ok=True)
    source_records = {s['file']:s for s in history['sources']}
    verified, rows, spectra, row_sources, initial = {}, [], [], [], {}
    for reference in history['rows']:
        n, dt, t = reference['N'], reference['dt'], reference['time']
        if t == 0 and n in initial:
            # The recorded full-field comparison establishes exact equality of
            # the 384-grid starting arrays. Reuse this one diagnostic too.
            metrics, spectrum, provenance = initial[n]
            provenance = dict(provenance, reused_initial_from_dt=.001)
        else:
            first = min(r['time'] for r in history['rows']
                        if (r['N'],r['dt']) == (n,dt) and r['time'] > 0)
            stop, key = (first, 'initial') if t == 0 else (t, 'final')
            path = cache/f'n{n}-dt{dt:g}-t{stop:g}.json'
            arrays_path = path.with_suffix('.npz')
            expected = source_records[path.name]
            if path.name not in verified:
                saved = json.loads(path.read_text())
                if (saved['N'],saved['dt'],saved['end_time'],saved['box'],saved['nu'],
                    saved['sigma'],saved['P_U'],saved['source_sha256']) != (n,dt,stop,6.,.01,0.,0.,source):
                    raise ValueError(f'Unexpected checkpoint metadata: {path}')
                if file_hash(path) != expected['summary_sha256'] or file_hash(arrays_path) != expected['arrays_sha256']:
                    raise ValueError(f'Checkpoint changed since the recorded shape analysis: {path}')
                verified[path.name] = expected
            provenance = dict(file=path.name, array_key=key, **{
                k:expected[k] for k in ('summary_sha256','arrays_sha256')})
            cached_path = analysis_cache/f'{path.stem}-{key}.json'
            cached = json.loads(cached_path.read_text()) if cached_path.exists() else {}
            if cached.get('analysis_source_sha256') == analysis_source and cached.get('source') == provenance:
                metrics, spectrum = cached['metrics'], cached['spectrum']
            else:
                with np.load(arrays_path, allow_pickle=False) as arrays:
                    velocity = arrays[key]
                    if velocity.shape != (3,n,n,n):
                        raise ValueError('Checkpoint grid differs from metadata')
                    metrics, spectrum = grid_spectrum(velocity, box=6.)
                    del velocity
                cached_path.write_text(json.dumps(dict(analysis_source_sha256=analysis_source,
                    source=provenance, metrics=metrics, spectrum=spectrum), allow_nan=False)+'\n')
            if t == 0:
                initial[n] = (metrics,spectrum,provenance)
        if not np.isclose(metrics['energy'], reference['energy'], rtol=2e-12, atol=1e-13):
            raise ValueError('Energy differs from the published saved-field measurement')
        rows.append(dict(N=n,dt=dt,time=t,**metrics))
        spectra.append(dict(N=n,dt=dt,time=t,**spectrum))
        row_sources.append(dict(N=n,dt=dt,time=t,**provenance))
        print(f'N={n} dt={dt:g} t={t:g}: blind={metrics["checkerboard_energy_fraction"]:.3e}, '
              f'short={metrics["short_wave_energy_fraction"]:.3e}, '
              f'gradient RMS gap={metrics["gradient_rms_gap_percent"]:.6f}%', flush=True)
    result = dict(start_time=0.,end_time=.4,external_force=0.,new_evolution_steps=0,
        numerical_source_sha256=source,analysis_source_sha256=analysis_source,
        history_file=history_path.name,history_sha256=file_hash(history_path),
        flow_file=flow_path.name,flow_sha256=file_hash(flow_path),
        rows=rows,spectra=spectra,sources=list(verified.values()),row_sources=row_sources,
        observations=len(rows),distinct_fields_measured=len(rows)-1,reused_identical_initial_observations=1,
        definitions=dict(checkerboard='Nonzero modes with every coordinate index in {0,N/2}: seven joint centered-derivative blind modes.',
            nyquist_planes='Union of Fourier planes with any absolute coordinate mode equal to N/2; each coefficient counted once.',
            short_wave='Any absolute coordinate mode >= N/4, i.e. wavelength <= four cells in at least one direction.',
            short_wave_gradient='Fraction of the sum of squared adjacent-cell differences in the short-wave band.',
            gradients='RMS over grid cells of all nine velocity derivatives; centered (u[i+1]-u[i-1])/(2h) and adjacent-cell (u[i+1]-u[i])/h.',
            gradient_rms_gap_percent='100*(1-centered_gradient_rms/neighbor_gradient_rms); zero for a constant field.',
            spectra='Bins by max(|m_x|,|m_y|,|m_z|), not spherical shells. Fractions are over all three velocity components.',
            normalization='Unnormalized forward real FFT; double reduced-axis interior powers, once at DC and Nyquist; Parseval checked against direct grid energy.'),
        limits=['Measures sampled arrays only; cannot detect content that aliased before sampling or lies above grid resolution.',
            'The short-wave cutoff is an inspection band, not a universal accuracy threshold; its physical wavelength differs by grid.',
            'Neither discrete gradient is an exact continuum derivative. The gap includes ordinary finite-difference attenuation of resolved modes.',
            'Joint blind modes, Nyquist planes and the broader short-wave band are nested and must not be added together.',
            'Saved times do not resolve intervening events; floating-point tails near roundoff should not be interpreted as physical signals.',
            'No new fluid evolution, filtering, or external force was applied. This supplies no rigorous error bound or Clay proof.'])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with out.with_suffix('.csv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history',type=Path,default=Path('math/results/hug-shape-0.4.json'))
    parser.add_argument('--flow',type=Path,default=Path('math/results/hug-continued-0.4.json'))
    parser.add_argument('--cache',type=Path,default=Path('scratch/hug-refinement'))
    parser.add_argument('--analysis-cache',type=Path,default=Path('scratch/hug-resolution'))
    parser.add_argument('--out',type=Path,default=Path('math/results/hug-resolution.json'))
    args = parser.parse_args()
    study(args.history,args.flow,args.cache,args.analysis_cache,args.out)
