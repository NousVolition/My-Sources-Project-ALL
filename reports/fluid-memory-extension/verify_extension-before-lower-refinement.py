"""Audit every declared recording and executed source, including failed attempts."""
import json,subprocess,sys,zipfile,hashlib
from final_experiment import ROOT,P,digest,dump
from analysis import read
from production import tasks
FOLDERS=['data','spinup','fast-data','fast-spinup','production-data','production-spinup','final-data','final-spinup','common-probe-data','unforced-readout-data','common-unforced-data','independent-data','independent-full-data','independent-common-data']
def verify():
    done=json.loads((ROOT/'results'/'gpu-pipeline-complete.json').read_text())
    expected={k:len(v) for k,v in tasks().items()}
    if done['counts']!=expected:raise RuntimeError('Production counts differ from frozen matrix')
    for stage,count in expected.items():
        info=json.loads((ROOT/'results'/('execution-'+stage+'.json')).read_text())
        if info['completed']!=count or info['planned']!=count:raise RuntimeError('Incomplete stage '+stage)
    cpu=json.loads((ROOT/'results'/'cpu-pipeline-complete.json').read_text())
    if len(cpu['completed'])!=12:raise RuntimeError('Incomplete independent comparisons')
    extra=0
    if (ROOT/'fresh-grid-refinement.json').exists():
        plan=json.loads((ROOT/'fresh-grid-refinement.json').read_text());refined=json.loads((ROOT/'results'/'refinement-complete.json').read_text())
        if refined['plan_sha256']!=digest(ROOT/'fresh-grid-refinement.json') or refined['planned']!=len(plan['tasks']) or len(refined['completed'])!=len(plan['tasks']):raise RuntimeError('Refinement completion/provenance mismatch')
        if [r['task'] for r in refined['completed']]!=plan['tasks']:raise RuntimeError('Refinement task list mismatch')
        extra=len(plan['tasks'])
    result=subprocess.run([sys.executable,'-m','unittest','test_integer_waves','test_corrected_controls','test_corrected_common','test_inference','test_spatial','-v'],cwd=ROOT,capture_output=True,text=True)
    dump(ROOT/'results'/'final-tests.json',dict(exit_code=result.returncode,output=result.stdout+result.stderr,sources={f.name:digest(f) for f in ROOT.glob('test*.py')}))
    if result.returncode:raise RuntimeError('Final tests failed: '+result.stdout+result.stderr)
    records=[]
    for folder in FOLDERS:
        for f in sorted((ROOT/folder).glob('*.npz')):
            if not f.with_suffix('.json').exists():raise RuntimeError('Orphan raw recording '+str(f))
            _,m=read(f,[])
            with zipfile.ZipFile(f) as z:
                hashes={}
                for name in z.namelist():
                    h=hashlib.sha256()
                    with z.open(name) as src:
                        for chunk in iter(lambda:src.read(8*1024**2),b''):h.update(chunk)
                    hashes[name]=h.hexdigest()
            records.append(dict(file=f.relative_to(ROOT).as_posix(),bytes=f.stat().st_size,sha256=m['raw_sha256'],array_members_sha256=hashes))
    # A numerical failure can be a valid completed outcome. This verification
    # confirms truthful completion/integrity; it does not waive scientific gates.
    summary=json.loads((ROOT/'results'/'extension-summary.json').read_text())
    audit=dict(completion_and_integrity_passed=True,production_counts=expected,refinement_configurations=extra,independent_comparisons=12,recording_count=len(records),records=records,
        all_numerical_checks_passed=summary['all_numerical_checks_passed'],sources={f.name:digest(f) for f in ROOT.glob('*.py')},
        request_id='2300ecd1-d208-4695-9b39-62b93d705892')
    dump(ROOT/'results'/'verification.json',audit);print('Verified '+str(len(records))+' recordings; numerical outcomes remain separately reported',flush=True);return audit
if __name__=='__main__':verify()
