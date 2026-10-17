"""Prepare a verified local release after the running study finishes.

This prepares files only. GitHub publication is a separately verified step.
"""
import argparse,json,time
from pathlib import Path
from final_experiment import ROOT,dump,digest
from package_extension import package
from roundtrip_extension import roundtrip

def prepare(bundle,restored):
    deadline=time.monotonic()+86400
    while not (ROOT/'results'/'verification.json').exists():
        if (ROOT/'results'/'gpu-pipeline-failure.json').exists():raise RuntimeError('Numerical pipeline failed; no release prepared')
        if time.monotonic()>deadline:raise RuntimeError('Verified completion wait timed out')
        time.sleep(10)
    verification=json.loads((ROOT/'results'/'verification.json').read_text())
    if not verification['completion_and_integrity_passed']:raise RuntimeError('Study integrity did not pass')
    package(bundle);roundtrip(bundle,restored)
    dump(ROOT/'results'/'release-ready.json',dict(status='ready for external publication',controller_sha256=digest(ROOT/'prepare_release.py'),
        raw_roundtrip_passed=True,recording_count=verification['recording_count'],request_id='2300ecd1-d208-4695-9b39-62b93d705892'))
    # Include readiness after its successful proof, without repacking the raw tar.
    from package_extension import compact
    compact(bundle)
    print('Verified release files ready; GitHub upload remains',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--restored',type=Path,required=True);a=ap.parse_args();prepare(a.bundle.resolve(),a.restored.resolve())
