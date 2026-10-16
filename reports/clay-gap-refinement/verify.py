"""Verify sources, raw records, parents, audits, protocol, playback and manifest."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,numpy as np
R=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def files():return sorted(p for p in R.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='checksums.json')
def main():
 from audit_one import load
 a=read(R/'analysis.json');assert a['analysis_sha256']==sha(R/'analyze.py') and a['audit_one_sha256']==sha(R/'audit_one.py')
 paths=sorted([*(R/'baseline/data').glob('*.json'),*(R/'numerical/data').glob('*.json'),*(R/'solver/data').glob('*.json')]);records={p.stem:p for p in paths if 'parameters' in read(p)}
 assert len(a['runs'])==len(records)==a['attempts'] and set(records)=={r['name'] for r in a['runs']}
 assert all(r['geometry_passed'] and r['energy_passed'] for r in a['runs'])
 for row in a['runs']:
  p=records[row['name']];r=read(p);assert row['metadata_sha256']==sha(p) and row['raw_sha256']==sha(p.with_suffix('.npz'))
  assert any(all((f/n).exists() and sha(f/n)==h for n,h in r['source_hashes'].items()) for f in [R/v for v in ['baseline','baseline-fresh','numerical','solver']])
  c=read(R/'data/audits'/(r['name']+'.json'));assert c['key']==dict(audit_sha256=sha(R/'audit_one.py'),metadata_sha256=sha(p),raw_sha256=sha(p.with_suffix('.npz'))) and c['result']==row
  load(p)
 for p in (R/'parents').glob('*.json'):
  r=read(p);extract=r.get('checkpoint_extract')
  if not extract:continue
  src=R/extract['source'];assert sha(src)==extract['source_metadata_sha256'] and sha(src.with_suffix('.npz'))==extract['source_raw_sha256'];z=np.load(p.with_suffix('.npz'));v=np.load(src.with_suffix('.npz'));i=extract['index']
  for key in z.files:assert np.array_equal(z[key],v[key][:i+1] if key in ['state','time'] else v[key]),(p,key)
 protocol=read(R/'data/executed-protocol.json');assert len(protocol['attempts'])==len(records) and protocol['thread_profile']==dict(numpy=1,scipy=1)
 for r in protocol['attempts']:
  original=read(records[r['name']]);assert all(original[k]==v for k,v in r.items())
 for name,script in [('regularization-check.json','check_regularization.py'),('derivative-check.json','check_derivatives.py'),('quadrature-check.json','check_quadrature.py')]:
  c=read(R/'data'/name);assert c['passed'] and c['check_sha256']==sha(R/script) and all(sha(R/'solver'/n)==h for n,h in c['source_hashes'].items())
 assert read(R/'data/factor-check.json')['source_sha256']==sha(R/'numerical/work_blocks.py')
 assert read(R/'data/factor-check.json')['check_sha256']==sha(R/'check_factor.py')
 catalog=read(R/'data-release.json');assert catalog['archive_readback_passed']
 for f in catalog['files']:assert (R/f['path']).stat().st_size==f['bytes'] and sha(R/f['path'])==f['sha256']
 restore=read(R/'data/archive-restore-check.json');assert restore['passed'] and restore['fresh_restore'] and restore['matching_files_preserved'] and restore['files']==len(catalog['files']) and restore['archive_sha256']==catalog['archive']['sha256'] and restore['catalog_sha256']==sha(R/'data-release.json') and restore['restorer_sha256']==sha(R/'restore_data.py') and restore['check_sha256']==sha(R/'check_restore.py')
 upload=read(R/'data/remote-upload-check.json');assets={v['name']:v for v in upload['assets']};assert upload['passed'] and upload['tag']==catalog['tag'] and len(assets)==3
 assert assets[catalog['archive']['name']]['sha256']==catalog['archive']['sha256'] and assets[catalog['archive']['name']]['bytes']==catalog['archive']['bytes']
 for name,path in [('data-release.json',R/'data-release.json'),('SHA256SUMS.txt',R/'data/SHA256SUMS.txt')]:assert assets[name]['sha256']==sha(path) and assets[name]['bytes']==path.stat().st_size
 replays=list((R/'data').glob('replay-*.json'));assert replays
 for p in replays:
  q=read(p);src=R/q['record'];assert q['all_arrays_bitwise_equal'] and all(q['arrays'].values()) and q['source_metadata_sha256']==sha(src) and q['source_raw_sha256']==sha(src.with_suffix('.npz'))
 scenes=read(R/'data/playback.json');names=[('length100-168-one','length100-168-two'),('original168-long-one','original168-long-two'),('middle72-one','middle72-two'),('weak72-one','weak72-two-tighter' if 'weak72-two-tighter' in records else 'weak72-two')];assert len(scenes)==4
 for scene,pair in zip(scenes,names):
  A,B=[load(records[n]) for n in pair];assert scene['records']==list(pair) and scene['n']==A[0]['parameters']['angles'] and scene['layers']==A[0]['parameters']['layers']
  assert np.array_equal(scene['triangles'],A[1]['triangles']) and np.array_equal(scene['outline'],np.round(A[1]['outline'],6))
  for f in scene['frames']:
   for s,panel in zip([A,B],f['panels']):
    r,z,t,y,d=s;idx=np.flatnonzero(abs(t-f['time'])<1e-10);assert len(idx)==1;i=int(idx[0]);n=len(z['initial_positions']);assert np.array_equal(panel['xy'],np.round(y[i,:2*n].reshape(n,2),5)) and np.array_equal(panel['lorenz'],np.round(y[:i+1,-7:-4],6))
    assert panel['gap']==d[i]['gap_arclength_rms'] and panel['imprint']==d[i]['imprint']
 browser=read(R/'data/browser-check.json');assert browser['passed'] and browser['cases']==len(scenes) and browser['report_sha256']==sha(R/'report.html') and browser['test_sha256']==sha(R/'test_browser.cjs')
 manifest={p.relative_to(R).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in files()}
 if '--write-manifest' in sys.argv:(R/'checksums.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
 assert read(R/'checksums.json')==manifest
 print(json.dumps(dict(passed=True,files=len(manifest)+1,attempts=len(records),completed=a['completed'],viewer_cases=len(scenes),manifest_sha256=sha(R/'checksums.json'))))
if __name__=='__main__':main()
