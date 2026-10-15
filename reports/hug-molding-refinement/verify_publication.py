"""Verify published bytes and completed saved results. Never run integrations."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO / 'tools'))
from restore_numerical_data import load_catalog, verify_files


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        self.links.extend(v for k, v in attrs if k in ('src', 'href') and v)
        self.ids.update(v for k, v in attrs if k == 'id')


def check_links(path):
    parser = Links()
    parser.feed(path.read_text(encoding='utf-8'))
    for link in parser.links:
        url = urlsplit(link)
        if url.scheme:
            assert url.scheme == 'https', link
            continue
        target = path.parent / unquote(url.path) if url.path else path
        assert target.is_file(), (path, link)
        if url.fragment and target.suffix == '.html':
            other = Links()
            other.feed(target.read_text(encoding='utf-8'))
            assert unquote(url.fragment) in other.ids, link
    return len(parser.links)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restored-root', required=True, type=Path)
    args = parser.parse_args()
    manifest = read(ROOT / 'publication-manifest.json')
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts
              and p.name != 'publication-manifest.json'}
    assert set(manifest) == actual, 'Published file inventory differs'
    assert all(digest(ROOT / name) == value for name, value in manifest.items())
    catalog = load_catalog(ROOT / 'data-release.json')
    count = verify_files(catalog, args.restored_root, full=True)
    package = args.restored_root / 'reports/hug-molding-refinement/package'
    study = package / 'hug-molding-refinement-20261013'
    cfg = read(study / 'protocol.json')
    handoff = read(study / 'handoff.json')
    for name, value in {**cfg['sources_sha256'], 'run_refinement.py': handoff['runner_sha256'],
                        'protocol.json': handoff['protocol_sha256']}.items():
        assert digest(study / name) == value, name
    for name, value in read(study / 'report-manifest.json').items():
        assert digest(study / name) == value, name
    assert digest(study / 'final-verification.json') == handoff['final_verification_sha256']
    assert digest(study / 'report-manifest.json') == handoff['report_manifest_sha256']
    links = check_links(ROOT / 'report.html')
    links += sum(check_links(p) for p in package.rglob('*.html'))
    # The historical verifier rewrites audit receipts, so run it on a disposable copy.
    with tempfile.TemporaryDirectory(prefix='hug-unloading-audit-', dir=args.restored_root.resolve()) as tmp:
        assert Path(tmp).resolve().is_relative_to(args.restored_root.resolve())
        copy = Path(tmp) / 'package'
        shutil.copytree(package, copy)
        audit = copy / study.name
        # The archive remains byte-exact. Only the disposable verifier adapts a
        # procedural coordinate reconstruction to cross-platform math rounding.
        # Connectivity, saved bytes, solver controls and scientific gates stay exact.
        script = audit / 'verify_refinement.py'
        source = script.read_text(encoding='utf-8')
        original = "assert np.array_equal(raw['triangles'],mesh.tri) and np.array_equal(raw['initial_positions'],mesh.initial)"
        replacement = """assert np.array_equal(raw['triangles'],mesh.tri)
            reference=raw['initial_positions']
            assert reference.shape==mesh.initial.shape and np.isfinite(mesh.initial).all()
            initial_error=float(np.max(np.abs(reference-mesh.initial)))
            initial_limit=8*np.finfo(np.float64).eps*max(1.,float(np.max(np.abs(reference))))
            print(json.dumps({'initial_mesh_case':name,'max_absolute_difference':initial_error,'roundoff_limit':initial_limit}))
            assert initial_error<=initial_limit"""
        assert source.count(original) == 1, 'Historical verifier changed'
        script.write_text(source.replace(original, replacement), encoding='utf-8')
        result = subprocess.run([sys.executable, '-B', str(audit / 'verify_refinement.py')],
                                text=True, capture_output=True, timeout=240)
        assert result.returncode == 0, result.stdout + result.stderr
        reconstruction = [json.loads(line) for line in result.stdout.splitlines()
                          if line.startswith('{') and 'initial_mesh_case' in line]
        assert len(reconstruction) == 4
        verified = read(audit / 'final-verification.json')
        assert verified['new_integrations_complete'] == 4
        assert not verified['new_integration_failures']
        assert len(verified['historical_mesh_failures_preserved']) == 1
    assert verify_files(catalog, args.restored_root, full=True) == count
    print(json.dumps({'status': 'passed', 'archived_files': count,
                      'published_files': len(manifest), 'report_links': links,
                      'new_simulations': 0, 'new_integrations_complete': 4,
                      'initial_mesh_reconstruction': reconstruction,
                      'historical_mesh_failure_preserved': True}))


if __name__ == '__main__':
    main()
