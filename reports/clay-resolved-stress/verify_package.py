"""Verify raw evidence, check fingerprints, displayed states and file coverage."""
from pathlib import Path
import argparse,json,numpy as np
from analyze import load,sha
from select_mesh import select
ROOT=Path(__file__).resolve().parent
def read(p):return json.loads((ROOT/p).read_text())
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write-manifest',action='store_true');a=ap.parse_args();e=read('analysis.json');assert e['analysis_sha256']==sha(ROOT/'analyze.py') and e['all_geometry_passed'] and e['all_energy_passed'];by={r['name']:r for r in e['runs']}
 for r in e['runs']:
  p=ROOT/r['path'];assert sha(p)==r['metadata_sha256'] and sha(p.with_suffix('.npz'))==r['raw_sha256'];load(p)
 for f,s in [('sparse-checks.json','check_sparse.py'),('indexed-checks.json','check_indexed.py'),('shared-drive-checks.json','check_shared_drive.py'),('resume-checks.json','check_resume.py'),('acceptance-checks.json','test_acceptance.py'),('fine-derivative-checks.json','check_fine_derivatives.py'),('quadrature-checks.json','check_quadrature.py')]:
  v=read('data/'+f);assert v['all_passed'] and v['source_sha256']==sha(ROOT/s)
 assert read('data/acceptance-checks.json')['analysis_sha256']==sha(ROOT/'analyze.py')
 assert read('solver/data/validation.json')['all_passed']
 assert read('solver-versions/indexed/data/validation.json')['all_passed']
 assert read('solver-versions/shared/data/validation.json')['all_passed']
 assert read('solver-versions/shared-resume/data/validation.json')['all_passed']
 for checkpoint in (ROOT/'parents').glob('*.json'):
  r=json.loads(checkpoint.read_text());c=r.get('checkpoint_extract')
  if not c:continue
  p=ROOT/c['source_record'];original=json.loads(p.read_text());u=np.load(p.with_suffix('.npz'));v=np.load(checkpoint.with_suffix('.npz'));i=c['index']
  assert sha(p)==c['source_metadata_sha256'] and sha(p.with_suffix('.npz'))==c['source_raw_sha256']
  assert all(np.array_equal(v[k],u[k][:i+1] if k in ['state','time'] else u[k]) for k in u.files)
  assert r['diagnostics']==original['diagnostics'][:i+1] and r['audit']==original['audit'] and r['parameters']==original['parameters'] and r['initial']==original['initial'] and r['completed_end']==c['time']==u['time'][i]
 actual=[p for p in [*(ROOT/'solver/data').glob('*.json'),*(ROOT/'solver-versions').glob('*/data/*.json')] if 'parameters' in json.loads(p.read_text())]
 assert {p.relative_to(ROOT).as_posix() for p in actual}=={r['path'] for r in e['runs']},'Analysis omits a completed attempt'
 for p in (ROOT/'data').glob('*-resolution.json'):
  r=json.loads(p.read_text());assert r['driver_sha256']==sha(ROOT/'run_refinement.py') and r['accepted']
 for p in (ROOT/'data').glob('*-matched-decision.json'):
  r=json.loads(p.read_text());assert r['driver_sha256']==sha(ROOT/'run_matched_refinement.py')
  decision=select([ROOT/q for q in r['record_paths']]);assert all(r[k]==v for k,v in decision.items()),'Matched refinement decision differs from raw evidence'
 if (ROOT/'data/reproduction-check.json').exists():
  r=read('data/reproduction-check.json');assert r['all_passed'] and r['source_sha256']==sha(ROOT/'check_reproduction.py')
 if (ROOT/'data/provenance-correction.json').exists():
  c=read('data/provenance-correction.json');p=ROOT/c['record'];r=read(c['record']);archive=ROOT/r['metadata_correction']['original_metadata']
  assert sha(archive)==c['before_metadata_sha256']==r['metadata_correction']['original_metadata_sha256'] and sha(p)==c['after_metadata_sha256'] and sha(p.with_suffix('.npz'))==c['raw_sha256']==r['metadata_correction']['raw_sha256_unchanged']
  old=json.loads(archive.read_text());old['parent'].pop('upstream_commit');old['parent'].update(checkpoint_origin='local-recording',source_record='solver/data/fine120-two.json');assert {k:v for k,v in r.items() if k!='metadata_correction'}==old
 protocol=read('data/executed-protocol.json');assert len(protocol['attempts'])==len(by) and protocol['report_sha256']==sha(ROOT/'build_report.py')
 for r in protocol['attempts']:assert all(read(by[r['name']]['path'])[k]==v for k,v in r.items())
 fragment=(ROOT/'molding-fragment.html').read_text(encoding='utf-8');payload=json.loads(fragment.split('<script type="application/json" id="mh-data">')[1].split('</script>')[0]);frames=0
 assert len(fragment.encode())<1000000
 for name,record in payload['records'].items():
  r,z,t,y,d=load(ROOT/by[name]['path']);assert record['triangles']==z['triangles'].tolist() and record['initial']==z['initial_positions'].round(7).tolist();n=len(z['initial_positions'])
  for f in record['frames']:
   i=int(np.argmin(abs(t-f['time'])));assert abs(t[i]-f['time'])<1e-10 and f['positions']==y[i,:2*n].reshape(n,2).round(7).tolist() and f['lorenz']==y[i,-7:-4].round(5).tolist() and f['gap']==d[i]['gap_arclength_rms'];frames+=1
 for f in ['browser-verification.json','browser-verification-inline.json']:
  v=read('data/'+f);assert v['fragmentSha256']==sha(ROOT/'molding-fragment.html') and v['testSha256']==sha(ROOT/'test_playback.cjs') and not v['errors'] and v['comparisons']==len(payload['cases']) and v['selections']==3*len(payload['cases'])
  assert v['invalidSavedState'] and v['clockMismatch'] and v['reducedMotion'] and all(q['scroll']<=q['width'] for q in v['layouts'])
  if f=='browser-verification.json':assert v['sourceSha256']==sha(ROOT/'report.html')
 assert 'Results are generated after' not in (ROOT/'README.md').read_text()
 current={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='checksums.json'}
 if a.write_manifest:(ROOT/'checksums.json').write_text(json.dumps(current,indent=2)+'\n')
 assert read('checksums.json')==current
 print(json.dumps(dict(files=len(current),attempts=e['attempts'],completed=e['completed'],frames=frames,incomplete=[r['name'] for r in e['runs'] if r['failure']]),indent=2))
if __name__=='__main__':main()
