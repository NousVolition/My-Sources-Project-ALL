"""Check final file, audit and browser provenance before publication."""
from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    hashes=json.loads((ROOT/'checksums.json').read_text())
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='checksums.json'}
    assert actual==set(hashes),'File manifest is incomplete'
    for name,digest in hashes.items():assert sha(ROOT/name)==digest,name
    a=json.loads((ROOT/'data/analysis.json').read_text());rows=json.loads((ROOT/'data/all-runs.json').read_text())
    assert a['all_declared_attempts_recorded'] and not a['missing_recordings']
    assert a['configurations']==a['declared_configurations']==len(rows)
    assert a['all_completed_audits_passed']
    assert a['analysis_source_sha256']==sha(ROOT/'analyze_results.py')
    for v in a['audits']:assert v['raw_sha256']==sha(ROOT/'data'/(v['name']+'.npz'))
    fragment=ROOT/'molding-fragment.html';assert fragment.stat().st_size<1000000
    payload=json.loads(re.search(r'<script type="application/json" id="mh-data">(.*?)</script>',fragment.read_text(encoding='utf-8'),re.S)[1])
    for name in ['browser-verification-inline.json','browser-verification.json']:
        b=json.loads((ROOT/'data'/name).read_text());assert b['errors']==[]
        assert b['fragmentSha256']==sha(fragment) and b['testSha256']==sha(ROOT/'test_playback.cjs')
        assert b['comparisons']==len(payload['cases']) and b['selections']==3*len(payload['cases'])
        if name=='browser-verification.json':assert b['sourceSha256']==sha(ROOT/'report.html')
    print(json.dumps(dict(files=len(hashes),recorded_attempts=a['configurations'],completed=a['completed'],limited_attempts=len(a['failures']),known_unresolved_state_checks=[c['name'] for c in a['solver_checks'] if not c['passed'] and 'superseded_by' not in c],browser_comparisons=len(payload['cases']),verified=True)))
if __name__=='__main__':main()
