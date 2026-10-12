"""Verify numerical evidence, displayed states, checks and complete file coverage."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from analyze import load,sha
ROOT=Path(__file__).resolve().parent
def read(p):return json.loads((ROOT/p).read_text())
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write-manifest',action='store_true');a=ap.parse_args()
    evidence=read('analysis.json');assert evidence['analysis_sha256']==sha(ROOT/'analyze.py') and evidence['adaptive_sha256']==sha(ROOT/'adaptive_mesh.py')
    assert evidence['all_geometry_passed'] and evidence['all_energy_passed']
    by={r['name']:r for r in evidence['runs']}
    for r in evidence['runs']:
        path=ROOT/r['path'];assert sha(path)==r['metadata_sha256'] and sha(path.with_suffix('.npz'))==r['raw_sha256'];load(path)
    for file,source,key in [('data/adaptive-checks.json','test_adaptive.py','test_sha256'),('data/contact-checks.json','check_contact.py','test_sha256'),('data/stiff-jacobian-checks.json','diagnose_stiff.py','source_sha256')]:
        result=read(file);assert result['all_passed'] and result[key]==sha(ROOT/source)
    assert read('data/adaptive-checks.json')['driver_sha256']==sha(ROOT/'adaptive_mesh.py')
    assert read('solver/data/validation.json')['all_passed']
    if (ROOT/'data/reproduction-check.json').exists():
        replay=read('data/reproduction-check.json')
        assert replay['all_passed'] and replay['source_sha256']==sha(ROOT/'check_reproduction.py')
    if (ROOT/'data/recorded-standard-one-resolution.json').exists():
        assert read('data/recorded-standard-one-resolution.json')['driver_sha256']==sha(ROOT/'history/adaptive_mesh_v1.py')
    protocol=read('data/executed-protocol.json');assert len(protocol['attempts'])==len(by)
    assert protocol['report_sha256']==sha(ROOT/'build_report.py')
    for c in protocol['attempts']:
        saved=read(by[c['name']]['path']);assert all(saved[k]==v for k,v in c.items())
    fragment=(ROOT/'molding-fragment.html').read_text(encoding='utf-8');payload=json.loads(fragment.split('<script type="application/json" id="mh-data">')[1].split('</script>')[0]);frames=0
    for name,record in payload['records'].items():
        path=ROOT/by[name]['path'] if name in by else ROOT/'parents'/(name+'.json');r,z,t,y,d=load(path)
        assert record['triangles']==z['triangles'].tolist() and record['initial']==z['initial_positions'].round(7).tolist()
        n=len(z['initial_positions'])
        for f in record['frames']:
            i=int(np.argmin(abs(t-f['time'])));assert abs(t[i]-f['time'])<1e-10
            assert f['positions']==y[i,:2*n].reshape(n,2).round(7).tolist()
            assert f['lorenz']==y[i,-7:-4].round(5).tolist()
            assert f['gap']==d[i]['gap_arclength_rms'];frames+=1
    for file in ['browser-verification.json','browser-verification-inline.json']:
        if file.endswith('-inline.json') and not (ROOT/'data'/file).exists():continue
        result=read('data/'+file)
        assert result['fragmentSha256']==sha(ROOT/'molding-fragment.html') and result['testSha256']==sha(ROOT/'test_playback.cjs')
        assert not result['errors'] and result['comparisons']==len(payload['cases']) and result['selections']==3*len(payload['cases'])
        assert result['invalidSavedState'] and result['clockMismatch'] and result['reducedMotion']
        assert all(x['scroll']<=x['width'] for x in result['layouts'])
        if file=='browser-verification.json':assert result['sourceSha256']==sha(ROOT/'report.html')
    assert 'Results are generated from the saved attempts' not in (ROOT/'README.md').read_text()
    current={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='checksums.json'}
    if a.write_manifest:(ROOT/'checksums.json').write_text(json.dumps(current,indent=2)+'\n')
    assert read('checksums.json')==current
    print(json.dumps(dict(files=len(current),new_attempts=evidence['attempts'],completed=evidence['completed'],verified_frames=frames,incomplete=[r['name'] for r in evidence['runs'] if r['failure']]),indent=2))
if __name__=='__main__':main()
