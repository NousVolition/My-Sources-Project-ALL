"""Fixed-state spatial refinement; parent NS, Heun and interpolant unchanged."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import sys
import json
import argparse
import time
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'mirror-direction'))
from run_matrix import R, prepare, cloud_target, source_matches
from run_pairs import Flow, advance, identity, sha, save_json, wrap


def transfer(h, n):
    """Copy Fourier modes with FFT normalization; supported spectra avoid Nyquist."""
    s = h.shape[1]
    if n == s:
        return h.copy()
    modes = np.rint(np.fft.fftfreq(s)*s).astype(int)
    keep = np.flatnonzero(abs(modes) < min(n, s)/2)
    z = np.arange(min(n, s)//2)
    out = np.zeros((3, n, n, n//2+1), dtype=h.dtype)
    out[:, modes[keep][:, None, None] % n,
        modes[keep][None, :, None] % n, z[None, None, :]] = (
        h[:, keep[:, None, None], keep[None, :, None], z[None, None, :]]*(n/s)**3)
    return out


def sources():
    paths = [HERE/'run_refinement.py', HERE/'protocol.json',
             HERE.parent/'mirror-direction'/'run_matrix.py',
             HERE.parent/'velocity-reversal'/'run_pairs.py',
             HERE.parent/'simulate.py', HERE.parent/'features.py',
             HERE.parent.parent/'numerics.py']
    return {str(x.relative_to(HERE.parent.parent)).replace('\\', '/'): sha(x) for x in paths}


def load_reference(seed):
    prior = HERE.parent/'velocity-reversal'/'recorded-data'
    name = f's{seed}_n56_dt0.005'
    checkpoint = prior/(name+'_branch-state.npz')
    meta = json.loads((prior/(name+'.json')).read_text())
    if sha(checkpoint) != meta['checkpoint_sha256']:
        raise ValueError('Changed starting checkpoint')
    for path, digest in meta['sources'].items():
        if not source_matches(HERE.parent.parent/path, digest):
            raise ValueError('Changed reference source '+path)
    with np.load(checkpoint, allow_pickle=False) as z:
        return z['field'].copy(), z['positions'].copy(), checkpoint


def run(seed, n, dt, p):
    name = f's{seed}_n{n}_dt{dt}'
    out = HERE/'recorded-data'; out.mkdir(exist_ok=True)
    dest = out/(name+'.npz'); meta_path = dest.with_suffix('.json')
    field_path = out/(name+'_endpoint-spectrum.npz')
    src = sources()
    if dest.exists() and meta_path.exists():
        old = json.loads(meta_path.read_text())
        if old['sources'] != src or sha(dest) != old['sha256'] or sha(field_path) != old['spectrum_sha256']:
            raise ValueError('Changed completed cache '+name)
        return old
    if dest.exists() or field_path.exists():
        raise ValueError('Unverified partial output preserved '+name)
    h56, x0, checkpoint = load_reference(seed)
    f = Flow(n, p['nu'], workers=1); f56 = Flow(56, p['nu'], workers=1)
    h0 = transfer(h56, n)
    energy_error = abs(f.inner(h0, h0)/f56.inner(h56, h56)-1)
    back_error = float(abs(transfer(h0, 56)-h56).max())
    if max(energy_error, back_error) > 1e-9:
        raise ValueError('Fourier state transfer failed')
    idx, r0, local = prepare(x0, p['neighbors']); m = len(x0)
    arrays = {'neighbor_indices': idx, 'initial_offsets': r0,
              'initial_local_mirror_offsets': r0@R}
    diagnostics = {}; spectra = {}; start = time.perf_counter()
    scratch = HERE/'recovery'; scratch.mkdir(exist_ok=True)
    recovery = scratch/(name+'_latest.npz')
    if recovery.exists():
        raise ValueError('Partial checkpoint exists; inspect before restart '+name)
    for sign, label in [(1, 'normal'), (-1, 'reversed')]:
        h = sign*h0.copy(); x = np.concatenate([x0, x0@R, local.reshape(-1, 3)])
        M = identity(len(x)); energy0 = f.inner(h, h)/2
        values = {key: [] for key in ['time', 'positions', 'tangent', 'global_positions',
                                      'global_tangent', 'local_neighbor_positions']}
        rows = []; diss = 0.; max_cfl = 0.
        try:
            for step in range(round(p['duration']/dt)+1):
                if step % round(p['save_dt']/dt) == 0:
                    vals = [step*dt, x[:m].copy(), M[:m].copy(), x[m:2*m].copy(),
                            M[m:2*m].copy(), x[2*m:].reshape(m, p['neighbors'], 3).copy()]
                    for key, value in zip(values, vals): values[key].append(value)
                    div = f.real(1j*sum(k*c for k, c in zip(f.k, h)))
                    rows.append({'elapsed': step*dt, 'energy': f.inner(h, h)/2,
                                 'energy_budget_error': (f.inner(h, h)/2-energy0+diss)/energy0,
                                 'spectral_divergence_max': float(abs(div).max()),
                                 'probe_tangent_volume_error_max': float(abs(np.linalg.det(M)-1).max()),
                                 'high_band_energy_fraction': f.inner(h*f.high, h*f.high)/f.inner(h, h)})
                    temp = recovery.with_suffix('.tmp.npz')
                    np.savez_compressed(temp, field=h, positions=x, tangent=M, step=step,
                                        dt=dt, label=label, dissipation=diss)
                    temp.replace(recovery)
                    save_json(recovery.with_suffix('.json'), {'sources': src, 'seed': seed, 'n': n,
                              'label': label, 'elapsed': step*dt, 'checkpoint_sha256': sha(recovery)})
                    print(f'{name} {label} saved elapsed {step*dt:.3f}', flush=True)
                if step == round(p['duration']/dt): break
                h, x, M, loss, cfl = advance(f, h, x, M, dt)
                diss += loss; max_cfl = max(max_cfl, cfl)
                if cfl >= .5 or not np.isfinite(h).all() or not np.isfinite(M).all():
                    raise FloatingPointError('Nonfinite state/CFL failure')
            if max(abs(r['energy_budget_error']) for r in rows) >= .005 or max(r['spectral_divergence_max'] for r in rows) >= 1e-10:
                raise FloatingPointError('Energy/divergence failure')
        except Exception as exc:
            save_json(scratch/(name+'_FAILURE.json'), {'error': repr(exc), 'sources': src, 'diagnostics': rows})
            raise
        arrays.update({label+'_'+key: np.asarray(value) for key, value in values.items()})
        spectra[label] = transfer(h, 56)/56**3
        diagnostics[label] = {'max_cfl': max_cfl, 'rows': rows}
    np.savez_compressed(dest, **arrays)
    np.savez_compressed(field_path, **spectra)
    meta = {'name': name, 'seed': seed, 'n': n, 'dt': dt, 'sources': src,
            'sha256': sha(dest), 'spectrum_sha256': sha(field_path),
            'checkpoint_path': str(checkpoint.relative_to(HERE.parent)).replace('\\', '/'),
            'checkpoint_sha256': sha(checkpoint), 'transfer_energy_relative_error': energy_error,
            'transfer_roundtrip_coefficient_error': back_error, 'diagnostics': diagnostics,
            'seconds': time.perf_counter()-start, 'complete': True}
    save_json(meta_path, meta)
    recovery.unlink(); recovery.with_suffix('.json').unlink()
    print(name+' pair complete', flush=True)
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pilot', action='store_true'); ap.add_argument('--extension', action='store_true')
    args = ap.parse_args(); p = json.loads((HERE/'protocol.json').read_text())
    jobs = [(s, n, p['dt']) for n in p['new_grids'] for s in p['seeds']]
    jobs += [(p['time_control_seed'], n, p['half_step_dt']) for n in p['time_control_grids']]
    if args.extension:
        decision = json.loads((HERE/'extension-decision.json').read_text())
        if not decision['extend']: raise ValueError('Extension not indicated by frozen screens')
        jobs = [(s, p['extension_grid'], p['dt']) for s in p['seeds']]
        jobs += [(p['time_control_seed'], p['extension_grid'], p['half_step_dt'])]
    if args.pilot: jobs = [(p['seeds'][0], p['pilot_grid'], p['dt'])]
    lock = HERE/'running.json'
    with lock.open('x') as handle: json.dump({'pid': os.getpid(), 'jobs': jobs, 'sources': sources()}, handle)
    records = []
    try:
        for s, n, dt in jobs:
            m = run(s, n, dt, p); records.append({k: m[k] for k in ['name', 'sha256', 'seconds', 'complete']})
            record_name = 'pilot.json' if args.pilot else ('extension-execution.json' if args.extension else 'execution.json')
            save_json(HERE/record_name, {'jobs': records, 'complete': len(records) == len(jobs), 'sources': sources()})
    finally:
        lock.unlink()


if __name__ == '__main__': main()
