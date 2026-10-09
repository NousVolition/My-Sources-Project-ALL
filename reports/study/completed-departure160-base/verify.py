"""Audit completed 160-grid departure runs from saved fields; perform no evolution."""
from pathlib import Path
from datetime import datetime, timezone
import sys, json, hashlib, gc
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'adversarial-vortex-study'
sys.path.insert(0, str(STUDY))
from solver import Solver
from protocol import source_hashes

IDS = ('exodus-fourier-n160-base',)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    checks = []
    hashes = source_hashes()
    for rid in IDS:
        folder = STUDY/'runs'/rid
        record = json.loads((folder/'result.json').read_text(encoding='utf-8'))
        assert record['status'] == 'complete' and len(record['rows']) == 21
        assert record['source_hashes'] == hashes
        assert np.allclose([r['t'] for r in record['rows']], np.arange(21)*.02, rtol=0, atol=1e-12)
        assert record['maximum_stage_cfl'] < .75
        assert max(abs(r['energy_balance_relative']) for r in record['rows']) < .005
        assert np.all(np.diff([r['I'] for r in record['rows']]) >= 0)
        n = record['job']['n']
        solver = Solver(n, record['job']['method'], record['job']['nu'], workers=1)
        rows = {round(r['t'], 8): r for r in record['rows']}
        fields = []
        for t in (0., .1, .2, .3, .4):
            field = folder/f'field-t{t:.2f}.npz'
            with np.load(field, allow_pickle=False) as a:
                assert float(a['t']) == t
                assert json.loads(str(a['job_json'])) == record['job']
                assert json.loads(str(a['source_hashes'])) == hashes
                h, markers, accum, initial = (a[k] for k in ('h', 'markers', 'accum', 'initial'))
            assert np.isfinite(h).all() and np.isfinite(markers).all()
            # Physical-space summation checks the spectral normalization used by the runner.
            u = np.fft.irfftn(h, s=(n,)*3, axes=(-3,-2,-1))
            w = np.fft.irfftn(solver.curl(h, physical=True), s=(n,)*3, axes=(-3,-2,-1))
            native = np.fft.irfftn(solver.curl(h), s=(n,)*3, axes=(-3,-2,-1))
            measured = {'Wmax': float(np.linalg.norm(w, axis=0).max()),
                        'energy': float(.5*np.sum(u*u)*solver.dx**3),
                        'spectral_enstrophy': float(.5*np.sum(w*w)*solver.dx**3),
                        'native_enstrophy': float(.5*np.sum(native*native)*solver.dx**3)}
            row = rows[t]
            errors = {k: abs(v-row[k])/max(abs(row[k]), 1) for k,v in measured.items()}
            assert max(errors.values()) < 1e-11
            assert np.max(abs(u.mean(axis=(1,2,3)))) < 1e-10
            assert accum[0] == row['I']
            assert np.isclose((measured['energy']-initial[0]+accum[2])/initial[0], row['energy_balance_relative'], atol=1e-12, rtol=0)
            assert np.isclose((measured['native_enstrophy']-initial[1]-accum[5]+accum[4])/initial[1], row['enstrophy_balance_relative'], atol=1e-12, rtol=0)
            fields.append({'file': field.name, 'sha256': sha(field), 'relative_errors': errors})
            del u,w,native
            if t == .4:
                final = solver.observe(h,t,markers,accum,initial)
                keys = ('high_band_enstrophy_fraction','direct_spectral_enstrophy_production','W_central_roi',
                        'origin_axial_strain','central_marker_parallel_stretch_spin_weighted')
                endpoint_errors = {k: abs(final[k]-row[k])/max(abs(row[k]),1) for k in keys}
                assert max(endpoint_errors.values()) < 1e-10
                for k in ('core_width','global_width'):
                    assert np.isclose(final[k]['minimum_chord_cells'],row[k]['minimum_chord_cells'],atol=1e-10,rtol=0)
            del h,markers,accum,initial
            gc.collect()
        checks.append({'id':rid,'result_sha256':sha(folder/'result.json'),
                       'fields':fields,'endpoint_errors':endpoint_errors,
                       'maximum_stage_cfl':record['maximum_stage_cfl'],
                       'maximum_absolute_energy_budget_relative':max(abs(r['energy_balance_relative']) for r in record['rows']),
                       'first_spectral_warning':next((r['t'] for r in record['rows'] if r['high_band_enstrophy_fraction']>.01),None),
                       'first_width_warning':next((r['t'] for r in record['rows'] if r['global_width']['minimum_chord_cells']<6),None)})
        print('Verified '+rid+' at all five retained field times.', flush=True)
        del solver
        gc.collect()
    assert hashes == source_hashes()
    result = {'status':'passed','checked_utc':datetime.now(timezone.utc).isoformat(),
              'new_simulations':0,'completed_runs':list(IDS),'saved_fields_checked':5*len(IDS),
              'source_hashes':hashes,'checks':checks,
              'scope':'Physical W, energy and both enstrophies independently summed from saved arrays at 0, .1, .2, .3, .4. Endpoint widths, strain and spectra remeasured with frozen diagnostic code. Integrated I and budget rates retain recorded RK-stage accumulations; trajectories were not rerun. Numerical completion does not imply spatial convergence.'}
    (HERE/'vortex-verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

if __name__ == '__main__':
    main()

