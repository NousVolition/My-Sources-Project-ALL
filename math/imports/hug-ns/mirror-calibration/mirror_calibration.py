"""Signed-coordinate calibration of the supplied mirror pair and Gaussian bias.

This is a new, explicitly defined calibration test, not a reconstruction of
missing uploaded experiment drivers. The existing matched-stretch suite is untouched.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import time
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE/'mirror-calibration'
OUT.mkdir(exist_ok=True)
SOURCE = HERE/'clay_hug.py'
spec = importlib.util.spec_from_file_location('mirror_flow_source', SOURCE)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
m.NU = .01


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for attempt in range(50):
        try:
            temporary.replace(path); return
        except PermissionError:
            if attempt == 49: raise
            time.sleep(.1)


def mirror(u):
    ref = (-np.arange(u.shape[1])) % u.shape[1]
    out = u[:, ref].copy(); out[0] *= -1
    return out


def inner(a, b, dx):
    return float(np.sum(a*b)*dx**3)


def build(n):
    x, y, z, kx, ky, kz, k2, keep, dx = m.grids(n)
    psi = np.exp(-12*((x-.9)**2+y*y))*(x > 0)
    raw = np.array([m.deriv(psi, ky, keep), -m.deriv(psi, kx, keep), np.zeros_like(psi)])
    symmetric = .5*(raw+mirror(raw))
    symmetric = np.array(m.project(symmetric, kx, ky, kz, k2, keep))
    symmetric = .5*(symmetric+mirror(symmetric))
    gaussian = np.zeros_like(symmetric)
    gaussian[1] = .4*np.exp(-10*((x-.55)**2+(y-.45)**2+z*z))
    projected = np.array(m.project(gaussian, kx, ky, kz, k2, keep))
    phi = .5*(projected-mirror(projected))
    phi /= np.sqrt(inner(phi, phi, dx))
    base_norm = np.sqrt(inner(symmetric, symmetric, dx))
    assert np.max(abs(symmetric-mirror(symmetric))) < 1e-12
    assert np.max(abs(phi+mirror(phi))) < 1e-12
    assert m.div_max(phi, kx, ky, kz, keep) < 1e-10
    assert abs(inner(symmetric, phi, dx)) < 1e-10
    return symmetric, phi, base_norm, (kx, ky, kz, k2, keep, dx)


def measure(u, phi, base_norm, ops, t):
    dx = ops[-1]
    asym = .5*(u-mirror(u))
    norm_u = np.sqrt(inner(u, u, dx))
    a = inner(asym, phi, dx)/inner(phi, phi, dx)
    residual = asym-a*phi
    asym_sq = inner(asym, asym, dx); residual_sq = inner(residual, residual, dx)
    E = 2*np.sqrt(asym_sq)/norm_u
    reconstruction = 2*np.sqrt(a*a*inner(phi, phi, dx)+residual_sq)/norm_u
    return {'t': t, 'a_signed': a, 's_fraction_of_initial_norm': a/base_norm,
            'E_mirror_mismatch': E, 'E_single_mode_prediction': 2*abs(a)/norm_u,
            'residual_asymmetry_fraction': float(np.sqrt(residual_sq/asym_sq)) if asym_sq > 1e-24*base_norm**2 else None,
            'exact_decomposition_error': abs(E-reconstruction),
            'flow_l2_norm': norm_u, 'energy': .5*norm_u**2,
            'divergence_max': m.div_max(u, ops[0], ops[1], ops[2], ops[4])}


def main():
    planned = [dict(n=n, half=half, amplitude=amplitude) for n in (33, 49)
               for half in (False, True) for amplitude in (0., .001, -.001, .005, -.005)]
    protocol = {'purpose': 'Define and check a signed fluid mirror coordinate', 'nu': .01, 'L': 6.,
        'end_time': .2, 'output_dt': .02, 'jobs': planned,
        'baseline': 'Exact supplied four_node.py mirror-pair start, repeated at33 and49',
        'template': 'Project supplied removal_test.py Gaussian velocity bias, take (g-Mg)/2, normalize its spatial L2 norm to1',
        'initial_condition': 'u0 = u_sym + amplitude * norm(u_sym) * phi',
        'signed_coordinate': 'a = <(u-Mu)/2,phi>/<phi,phi>; s_calibrated = a/norm(u_sym at0)',
        'unsigned_coordinate': 'E = norm(u-Mu)/norm(u)',
        'exact_identity': 'E^2=4*(a^2*norm(phi)^2+norm(residual)^2)/norm(u)^2, residual=(u-Mu)/2-a*phi',
        'initial_exact_values': 's_calibrated(0)=amplitude; E(0)=2*abs(amplitude)/sqrt(1+amplitude^2)',
        'scope': 'Calibration and reflection checks; does not fit or validate the cubic reduced evolution law',
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    existing = OUT/'protocol.json'
    if existing.exists():
        assert json.loads(existing.read_text()) == protocol, 'Calibration protocol changed'
    else: save(existing, protocol)
    completed = []
    for job in planned:
        name = f'n{job["n"]}-'+('half' if job['half'] else 'base')+f'-a{job["amplitude"]}'
        folder = OUT/name; folder.mkdir(exist_ok=True)
        path = folder/'result.json'
        if path.exists() and json.loads(path.read_text())['status'] == 'complete':
            completed.append(name); continue
        save(OUT/'task-state.json', {'pid': os.getpid(), 'status': 'running', 'completed': completed,
                                   'active': name, 'planned': len(planned), 'updated_epoch': time.time()})
        base, phi, base_norm, ops = build(job['n'])
        u = base+job['amplitude']*base_norm*phi
        maxspeed = float(np.sqrt(np.sum(base*base, axis=0)).max())
        substeps = int(np.ceil(.02/(.04*ops[-1]/maxspeed)))*(2 if job['half'] else 1)
        dt = .02/substeps
        rows = []; start_index = 0
        if path.exists():
            previous = json.loads(path.read_text())
            assert previous['dt'] == dt
            rows = previous['series']; start_index = len(rows)
            u = np.load(folder/previous['field'])
        for i in range(start_index, 11):
            if i:
                for _ in range(substeps):
                    u = np.array(m.advance(u, dt, ops))
                    assert np.isfinite(u).all()
            row = measure(u, phi, base_norm, ops, i*.02)
            assert row['divergence_max'] < 1e-9
            assert row['exact_decomposition_error'] < 1e-12
            if i == 0:
                assert abs(row['s_fraction_of_initial_norm']-job['amplitude']) < 1e-12
                assert abs(row['E_mirror_mismatch']-2*abs(job['amplitude'])/np.sqrt(1+job['amplitude']**2)) < 1e-12
            rows.append(row)
            field = f'field-{i:02d}.npy'; np.save(folder/field, u)
            save(path, {'status': 'complete' if i == 10 else 'running', 'job': job, 'dt': dt,
                        'series': rows, 'field': field, 'updated_epoch': time.time()})
            print(f'{name} t={i*.02:.2f} signed={row["s_fraction_of_initial_norm"]:.7g} E={row["E_mirror_mismatch"]:.7g}', flush=True)
        completed.append(name)
    save(OUT/'task-state.json', {'pid': os.getpid(), 'status': 'complete', 'completed': completed,
                               'active': None, 'planned': len(planned), 'updated_epoch': time.time()})


if __name__ == '__main__': main()
