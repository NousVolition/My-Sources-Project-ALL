"""Audit completed matrix provenance and matched initial observations, without integration."""
from pathlib import Path
import json
import numpy as np
from run_refinement import sha, source_matches

HERE = Path(__file__).resolve().parent


def main():
    p = json.loads((HERE/'protocol.json').read_text())
    done = json.loads((HERE/'completion.json').read_text())
    r = json.loads((HERE/'results.json').read_text())
    assert done['complete'] and done['results_sha256'] == sha(HERE/'results.json')
    assert done['report_sha256'] == sha(HERE/'README.md')
    assert done['analysis_sha256'] == sha(HERE/'analyze.py')
    assert p == r['protocol']
    assert len(r['runs']) == (23 if done['extended'] else 18)
    expected = {(seed, n, p['dt']) for seed in p['seeds'] for n in r['grids']}
    expected |= {(p['time_control_seed'], n, p['half_step_dt']) for n in p['time_control_grids']}
    if done['extended']:
        expected.add((p['time_control_seed'], p['extension_grid'], p['half_step_dt']))
    actual = {(v['seed'], v['n'], v['dt']) for v in r['runs'].values()}
    assert actual == expected
    reference = {}
    for seed in p['seeds']:
        path = HERE.parent/'velocity-reversal'/'recorded-data'/f's{seed}_n56_dt0.005_branch-state.npz'
        with np.load(path, allow_pickle=False) as z:
            reference[seed] = z['positions'].copy()
    centers = {}; geometry = {}; records = []
    R = np.diag([-1., 1., 1.])
    times = np.arange(round(p['duration']/p['save_dt'])+1)*p['save_dt']
    for name, record in r['runs'].items():
        path = HERE/'recorded-data'/(name+'.npz')
        m = json.loads(path.with_suffix('.json').read_text())
        assert m['complete'] and m['sha256'] == record['data_sha256'] == sha(path)
        assert m['spectrum_sha256'] == sha(path.with_name(name+'_endpoint-spectrum.npz'))
        assert m['checkpoint_sha256'] == sha(HERE.parent/m['checkpoint_path'])
        assert m['sources'] == done['sources']
        for rel, digest in m['sources'].items():
            assert source_matches(HERE.parent.parent/rel, digest), rel
        with np.load(path, allow_pickle=False) as z:
            assert all(np.isfinite(z[k]).all() for k in z.files)
            x = reference[m['seed']]
            assert x.shape == (p['markers'], 3)
            assert np.array_equal(z['normal_positions'][0], x)
            assert z['neighbor_indices'].shape == (p['markers'], p['neighbors'])
            assert np.array_equal(z['initial_local_mirror_offsets'], z['initial_offsets']@R)
            for label in ['normal', 'reversed']:
                np.testing.assert_allclose(z[label+'_time'], times, atol=2e-16, rtol=0)
                assert np.array_equal(z[label+'_positions'][0], x)
                assert np.array_equal(z[label+'_global_positions'][0], x@R)
                assert np.array_equal(z[label+'_local_neighbor_positions'][0], x[:, None]+z['initial_offsets']@R)
                for prefix in ['', 'global_']:
                    assert np.array_equal(z[label+'_'+prefix+'tangent'][0], np.broadcast_to(np.eye(3), (p['markers'], 3, 3)))
            g = tuple(z[k].copy() for k in ['neighbor_indices', 'initial_offsets', 'initial_local_mirror_offsets'])
            if m['seed'] in geometry:
                assert all(np.array_equal(a, b) for a, b in zip(g, geometry[m['seed']]))
            else:
                geometry[m['seed']] = g
        records.append({'name': name, 'trajectory_sha256': m['sha256'], 'spectrum_sha256': m['spectrum_sha256'],
                        'reference_checkpoint_sha256': m['checkpoint_sha256'], 'matched_initial_data': True})
    with np.load(HERE/'targets.npz', allow_pickle=False) as z:
        assert len(z.files) == len(r['runs'])*len(p['horizons'])*3*4
        assert all(z[k].shape == (p['markers'],) and np.isfinite(z[k]).all() for k in z.files)
    report = {'complete': True, 'pairs': len(records), 'forward_evolutions': 2*len(records),
              'identical_centers_neighbors_offsets_and_initial_tangents_across_grids': True,
              'all_seed_grid_step_combinations_match_protocol': True,
              'all_saved_times_and_targets_verified': True, 'records': records,
              'results_sha256': sha(HERE/'results.json'), 'verifier_sha256': sha(Path(__file__))}
    (HERE/'data-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'records'}))


if __name__ == '__main__':
    main()
