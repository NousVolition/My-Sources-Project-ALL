"""Check delivered evidence and exact file coverage without rerunning integrations."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from accuracy_gate import assess

ROOT = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))

def files():
    return {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*')
            if p.is_file() and '__pycache__' not in p.parts
            and p.suffix != '.pyc' and p.name != 'checksums.json'}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    analysis = read('analysis.json')
    assert analysis['analysis_sha256'] == digest(ROOT / 'analyze.py')
    assert analysis['all_geometry_passed'] and analysis['all_energy_passed']
    by = {r['key']: r for r in analysis['runs']}
    for row in analysis['runs']:
        path = ROOT / row['relative_path']
        record = json.loads(path.read_text())
        assert digest(path.with_suffix('.npz')) == row['raw_sha256'], row['key']
        for name, sha in record['source_hashes'].items():
            assert digest(path.parent.parent / name) == sha, (row['key'], name)
    final = [r for r in analysis['runs'] if r['key'].startswith('solver/')]
    protocol = read('data/executed-protocol.json')
    assert {r['name'] for r in final} == {r['name'] for r in protocol['attempts']}
    for record in protocol['attempts']:
        actual = read('solver/data/' + record['name'] + '.json')
        assert all(actual[key] == value for key, value in record.items())
    for name in ['validation', 'jacobian-checks', 'spatial-checks']:
        assert read('solver/data/' + name + '.json')['all_passed'], name
    gates = read('data/final-checks.json')['checks']
    for gate in gates:
        coarse, fine = ROOT / gate['coarse'], ROOT / gate['fine']
        assert digest(coarse) == gate['coarse_sha256']
        assert digest(fine) == gate['fine_sha256']
        assert gate['checker_sha256'] == digest(ROOT / 'finalize_checks.py')
        actual = assess(coarse, fine, gate['kind'])
        for key, value in actual.items():
            if key not in ['coarse', 'fine']:
                assert gate[key] == value, (gate['name'], key)
    # The deliberately inadequate coarse preview must remain visibly rejected.
    assert not next(g for g in gates if g['name'] == 'gate-final-coarse-preview')['passed']
    assert all(g['passed'] for g in gates if g['name'] in
               ['gate-final-one', 'gate-final-two', 'gate-time-one', 'gate-time-two'])
    for comparison in analysis['solver_comparisons']:
        if comparison['a'].startswith('solver/'):
            assert comparison['passed'] and comparison['full_requested_interval']
    fragment = (ROOT / 'molding-fragment.html').read_text(encoding='utf-8')
    payload = json.loads(fragment.split('<script type="application/json" id="mh-data">')[1].split('</script>')[0])
    for key, record in payload['records'].items():
        raw = np.load((ROOT / by[key]['relative_path']).with_suffix('.npz'))
        assert record['triangles'] == raw['triangles'].tolist()
        assert record['initial'] == raw['initial_positions'].round(5).tolist()
        for frame in record['frames']:
            index = int(np.argmin(abs(raw['time'] - frame['time'])))
            assert abs(raw['time'][index] - frame['time']) < 1e-10
            count = len(raw['initial_positions'])
            positions = raw['state'][index, :2 * count].reshape(count, 2)
            assert positions.round(5).tolist() == frame['positions']
            assert raw['state'][index, -7:-4].round(5).tolist() == frame['lorenz']
    browser = read('data/browser-verification.json')
    assert browser['sourceSha256'] == digest(ROOT / 'report.html')
    for name in ['browser-verification.json', 'browser-verification-inline.json']:
        if name.endswith('-inline.json') and not (ROOT / 'data' / name).exists():
            continue  # Standalone reproduction does not require the Codex host.
        result = read('data/' + name)
        assert result['fragmentSha256'] == digest(ROOT / 'molding-fragment.html')
        assert result['testSha256'] == digest(ROOT / 'test_playback.cjs')
        assert not result['errors'] and result['comparisons'] == len(payload['cases'])
        assert result['selections'] == 3 * len(payload['cases'])
        assert result['clockMismatch'] and result['invalidSavedState'] and result['reducedMotion']
        assert all(x['scroll'] <= x['width'] for x in result['layouts'])
    assert 'Results are generated from saved states' not in (ROOT / 'README.md').read_text()
    current = {name: digest(path) for name, path in sorted(files().items())}
    if args.write_manifest:
        (ROOT / 'checksums.json').write_text(json.dumps(current, indent=2) + '\n', encoding='utf-8')
    assert read('checksums.json') == current, 'Manifest must cover the exact current file set'
    print(json.dumps(dict(files=len(current), released_attempts=len(final),
        released_completed=sum(r['failure'] is None for r in final),
        verified_recorded_frames=sum(len(r['frames']) for r in payload['records'].values()),
        checks=[dict(name=g['name'], passed=g['passed']) for g in gates]), indent=2))

if __name__ == '__main__':
    main()
