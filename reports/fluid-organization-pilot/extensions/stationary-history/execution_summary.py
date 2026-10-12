"""Describe actual acquisition locations without treating replicas as independent."""
from common import ROOT,protocol,read,save
import xml.etree.ElementTree as ET
p=protocol();cloud=ROOT/'data/cloud';cal=[];test=[];controls=[];records=[]
for path in sorted(cloud.glob('*_environment.json')):
 q=read(str(path.relative_to(ROOT)));records.append(q)
 if path.name.startswith('group_'):
  imported=path.with_name(path.name.replace('_environment','_import'))
  if imported.exists():cal.extend(q['seeds'])
 elif q.get('action')=='base':test.extend(q['seeds'])
 elif q.get('action')=='control':controls.extend((s,q['label']) for s in q['seeds'])
assert len(cal)==len(set(cal)) and len(test)==len(set(test)) and len(controls)==len(set(controls))
expected_cal=p['splits']['train']+p['splits']['validation'];expected_test=p['splits']['test']
def tests_pass(name):
 path=ROOT/'results'/name
 if not path.exists():return False
 suites=list(ET.parse(path).getroot().iter('testsuite'))
 return sum(int(s.attrib['tests']) for s in suites)==20 and all(int(s.attrib['failures'])==0 and int(s.attrib['errors'])==0 for s in suites)
save('results/execution.json',dict(calibration_local_seeds=[s for s in expected_cal if s not in cal],calibration_hosted_seeds=sorted(cal),test_hosted_seeds=sorted(test),test_local_seeds=[s for s in expected_test if s not in test],hosted_control_pairs=sorted(controls),hosted_environments=records,local_environment=read('results/environment.json'),local_20_test_pass=tests_pass('tests.xml'),linux_20_test_pass=tests_pass('linux_tests.xml'),selection='Hosted calibration groups were imported before the corresponding local job started. Local queued jobs then used the hash-verified cache. No duplicated seed is counted as another independent run.',cross_platform_note='Original acquisition bytes are retained. Different platforms can produce small floating-point differences in regenerated trajectories.'))
print('Execution provenance recorded.')
