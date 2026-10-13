"""Verify actual hosted acquisition and model-freeze chronology."""
import datetime
import xml.etree.ElementTree as ET
from common import ROOT,read,save,protocol,digest

def test_receipt(name):
    suites=list(ET.parse(ROOT/'results'/name).getroot().iter('testsuite'))
    n=sum(int(s.attrib['tests']) for s in suites)
    assert n==22 and all(int(s.attrib['failures'])==0 and int(s.attrib['errors'])==0 for s in suites)
    return n

def main():
    p=protocol();records=[]
    for path in sorted((ROOT/'data/cloud').glob('*_environment.json')):records.append(read(str(path.relative_to(ROOT))))
    for stage,n in [('pilot',4),('calibration',32),('test',32),('control',8)]:
        selected=[q for q in records if q['stage']==stage]
        assert len(selected)==n and set(q['index'] for q in selected)==set(range(n)),(stage,len(selected))
    freeze=read('results/freeze_receipt.json');frozen_at=datetime.datetime.fromisoformat(freeze['created_utc'])
    for seed in p['splits']['test']:
        meta=read(f'data/runs/base_{seed}.json');assert datetime.datetime.fromisoformat(meta['started_utc'])>frozen_at
    for q in records:
        if q['stage']=='test':assert q['frozen_models_sha256']==freeze['sha256']
    commits={q['github_sha'] for q in records};assert len(commits)==1
    save('results/execution.json',dict(acquisition_commit=next(iter(commits)),workflow_run_ids=sorted({q['github_run_id'] for q in records}),calibration_hosted_seeds=p['splits']['train']+p['splits']['validation'],test_hosted_seeds=p['splits']['test'],matched_control_pairs=[[s,label] for s in p['refinement_seeds'] for label in ['half','fine']],hosted_environments=records,local_test_count=test_receipt('tests.xml'),linux_test_count=test_receipt('linux_tests.xml'),model_freeze_sha256=freeze['sha256'],model_freeze_utc=freeze['created_utc'],all_test_starts_after_freeze=True,selection='Every predeclared run retained. Whole independent baseline fields are sampling units; snapshots, strengths and dependent controls are not extra independent runs.'))
    print('Acquisition and freeze chronology verified')

if __name__=='__main__':main()
