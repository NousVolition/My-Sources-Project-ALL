"""Restore release parts into an isolated directory and repeat the analysis.

No physical simulations are rerun. Exact scientific outputs must match before
the data package is considered ready for publication.
"""
import argparse,json,subprocess,sys,zipfile
from pathlib import Path
from final_experiment import ROOT,dump,digest
from package_extension import compact

SCIENTIFIC=['main-analysis.json','common-analysis.json','unforced-analysis.json',
    'controls-analysis.json','extension-summary.json','spatial-analysis.json']

def roundtrip(bundle,out):
    if out.exists() and any(out.iterdir()):raise RuntimeError('Roundtrip destination must be empty')
    out.mkdir(parents=True,exist_ok=True);target=out.resolve()
    code=bundle/'fluid-memory-extension-code-and-results.zip'
    with zipfile.ZipFile(code) as z:
        if z.testzip():raise RuntimeError('Compact archive CRC failure')
        for name in z.namelist():
            dest=(target/name).resolve()
            if not dest.is_relative_to(target):raise RuntimeError('Unsafe compact archive member')
        z.extractall(target)
    # Force assembly from the release parts rather than reusing the original tar.
    subprocess.run([sys.executable,str(target/'restore_recordings.py'),'--parts',str(bundle),
        '--out',str(target),'--archive',str(target/'joined-recordings.tar.gz')],cwd=target,check=True)
    subprocess.run([sys.executable,str(target/'reproduce_extension.py'),'--analysis-only'],cwd=target,check=True)
    matches={}
    files=SCIENTIFIC+(['baseline-'+name for name in SCIENTIFIC if name!='spatial-analysis.json'] if (ROOT/'fresh-grid-refinement.json').exists() else [])
    for name in files:
        matches['results/'+name]=json.loads((ROOT/'results'/name).read_text())==json.loads((target/'results'/name).read_text())
    for name in ['README.md','report.html']:
        matches[name]=digest(ROOT/name)==digest(target/name)
    for p in (ROOT/'figures').glob('*.png'):
        matches[p.relative_to(ROOT).as_posix()]=digest(p)==digest(target/'figures'/p.name)
    if not all(matches.values()):raise RuntimeError('Restored outputs differ: '+json.dumps(matches))
    manifest=json.loads((ROOT/'raw-recordings.json').read_text())
    proof=dict(passed=True,scope='All release parts rejoined, every original recording checksum verified, all scientific analyses/report/figures regenerated exactly; no NS simulations repeated',
        recording_count=manifest['recording_count'],archive_sha256=manifest['archive_sha256'],matches=matches,
        sources={f:digest(ROOT/f) for f in ['roundtrip_extension.py','restore_recordings.py','reproduce_extension.py']},
        request_id='2300ecd1-d208-4695-9b39-62b93d705892')
    dump(ROOT/'results'/'raw-roundtrip.json',proof);dump(bundle/'raw-roundtrip.json',proof)
    # Include the now-completed proof in the downloadable compact bundle.
    compact(bundle);print('Raw release restoration and scientific-output roundtrip passed',flush=True)
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();roundtrip(a.bundle.resolve(),a.out.resolve())
