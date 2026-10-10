"""Verify saved experiment semantics, test status and artifact integrity."""
from pathlib import Path
import csv,hashlib,json,re,xml.etree.ElementTree as ET
import numpy as np
from model import to_state,from_state,voltage,power_residual
ROOT=Path(__file__).resolve().parent


def main():
    c=json.loads((ROOT/'protocol.json').read_text());s=json.loads((ROOT/'results/summary.json').read_text())
    ref=json.loads((ROOT/'results/reference_summary.json').read_text())
    if not ref['refinement']['all_passed'] or not ref['power_gate_passed']:raise AssertionError('Reference checks failed')
    checks=[]
    def check(name,ok):
        checks.append({'check':name,'passed':bool(ok)})
        if not ok:raise AssertionError(name)
    initial=np.load(ROOT/'data/initial_conditions.npz')['phases'];initial_state,sign=to_state(initial)
    for setting in c['settings']:
        with np.load(ROOT/f'data/{setting["name"]}.npz') as z:
            prefix=setting['name'];states=z['states'];starts=z['initial'];ends=z['endpoint'];settled=z['settled']
            check(prefix+' starts match recorded trajectory',np.array_equal(starts,states[:,0]))
            check(prefix+' fresh resets exact state',np.array_equal(starts[:,0],np.broadcast_to(initial_state,starts[:,0].shape)))
            check(prefix+' continuation from S+W endpoints',np.array_equal(starts[1:,1:],ends[:-1,1:]))
            check(prefix+' initial sweeps cloned correctly',np.array_equal(starts[0],np.broadcast_to(initial_state,starts[0].shape)))
            check(prefix+' settled states occur at declared time',np.array_equal(settled,states[:,round(setting['settle']/c['record_stride'])]))
            check(prefix+' end states exact',np.array_equal(ends,states[:,-1]))
            expected=voltage(from_state(settled,z['sign']),from_state(ends,z['sign']),c['base_window'])
            check(prefix+' saved individual voltage recomputes',np.array_equal(expected,z['voltage']))
            check(prefix+' actual current order',np.array_equal(z['bias_order'],np.array([c['biases'],c['biases'],c['biases'][::-1]]).T))
            ib=z['bias_order'][:,None,:,None,None]
            phases=from_state(states,z['sign']);d=z['dissipated']
            residual=power_residual(phases.transpose(1,0,2,3,4,5),d.transpose(1,0,2,3,4),z['bias_order'][:,:,None,None])
            check(prefix+' power residual independently recomputes',np.allclose(residual,z['residual'],rtol=0,atol=1e-13))
    for settle in [100,200]:
        with np.load(ROOT/f'data/reference_S{settle}.npz') as z:
            prefix=f'reference S={settle}';starts=z['initial'];ends=z['endpoint'];states=z['states']
            check(prefix+' initial draws unchanged',np.array_equal(starts[0],np.broadcast_to(initial_state,starts[0].shape)))
            check(prefix+' fresh reset',np.array_equal(starts[:,0],np.broadcast_to(initial_state,starts[:,0].shape)))
            check(prefix+' continuation uses S+200 endpoint',np.array_equal(starts[1:,1:],ends[:-1,1:]))
            check(prefix+' settled states at declared time',np.array_equal(z['settled'],states[:,settle//5]))
            expected=voltage(from_state(z['settled'],z['sign']),from_state(ends,z['sign']),200)
            check(prefix+' measured voltage recomputes',np.array_equal(expected,z['voltage']))
            check(prefix+' finite and dissipative',np.isfinite(states).all() and np.all(np.diff(z['dissipated'],axis=1)>=-1e-10))
    with np.load(ROOT/'data/reference_S200.npz') as b,np.load(ROOT/'data/reference_windows.npz') as w:
        check('reference windows reuse exact settled states',np.array_equal(b['settled'],w['states'][0]))
    with np.load(ROOT/'data/reference_windows.npz') as a,np.load(ROOT/'data/reference_windows_RK4.npz') as b:
        check('cross-solver windows reuse exact state and sign',np.array_equal(a['states'][0],b['states'][0]) and np.array_equal(a['sign'],b['sign']))
    with np.load(ROOT/'data/both.npz') as base,np.load(ROOT/'data/windows.npz') as w:
        check('window replay starts from identical settled states',np.array_equal(base['settled'],w['states'][0]))
        check('window replay preserves hidden log separation',np.array_equal(base['sign'],w['sign']))
        check('window endpoint W=200 matches primary',np.array_equal(base['endpoint'],w['states'][round(200/c['record_stride'])]))
    for path,count in [('primary_voltages.csv',7344),('window_voltages.csv',5508),('disjoint_blocks.csv',7344)]:
        with (ROOT/'results'/path).open() as f:rows=list(csv.DictReader(f))
        check(path+' row count',len(rows)==count)
        if 'v1' in rows[0]:check(path+' total equals individual sum',all(abs(float(r['v1'])+float(r['v2'])-float(r['total_voltage']))<1e-12 for r in rows))
    tests=ET.parse(ROOT/'tests.xml').getroot();suites=tests.findall('testsuite')
    check('18 automated tests passed',sum(int(x.attrib['tests']) for x in suites)==18 and all(int(x.attrib['failures'])==int(x.attrib['errors'])==0 for x in suites))
    for link in re.findall(r'(?:src|href)="([^"]+)"',(ROOT/'report.html').read_text(encoding='utf-8')):
        if '://' not in link and not link.startswith('#'):check('report link '+link,(ROOT/link).exists())
    (ROOT/'results/artifact_checks.json').write_text(json.dumps({'all_passed':True,'checks':checks},indent=2)+'\n')
    final={'status':'Completed with documented initial numerical failures and passing additional refinement checks',
           'original_four_setting_step_gate_passed':s['numerical_gates']['max_absolute_voltage_step_difference']['passed'],
           'additional_refinement_gates_passed':ref['refinement']['all_passed'],
           'reference_power_gate_passed':ref['power_gate_passed'],
           'automated_tests_passed':18,'saved_artifact_checks_passed':len(checks),
           'scope':'Refined full-history references support finite-time statements only. Original fixed-step gate remains failed; no asymptotic hysteresis or new physical memory law is established.'}
    (ROOT/'results/final_verification.json').write_text(json.dumps(final,indent=2)+'\n')
    files=[]
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or '.pytest_cache' in path.parts:continue
        rel=path.relative_to(ROOT)
        if rel.parts[0].startswith('preliminary_direct_') or rel.as_posix()=='manifest.json':continue
        b=path.read_bytes();files.append({'path':rel.as_posix(),'size':len(b),'sha256':hashlib.sha256(b).hexdigest()})
    (ROOT/'manifest.json').write_text(json.dumps({'scope':'Final accepted source, protocol, tests, report and arrays; archived initial failed-check summaries included; superseded direct-phase raw traces excluded.','files':files},indent=2)+'\n')
    print(f'{len(checks)} saved-artifact checks passed; {len(files)} files hashed')


if __name__=='__main__':main()
