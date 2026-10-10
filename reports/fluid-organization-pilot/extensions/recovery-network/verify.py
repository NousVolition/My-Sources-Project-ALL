"""Fail visibly when a declared numerical/physical gate is not satisfied."""
import hashlib,json,platform,sys
import xml.etree.ElementTree as ET
import numpy as np
import scipy,matplotlib
from common import ROOT,P,save

def load(name):return json.loads((ROOT/'results'/name).read_text(encoding='utf-8'))
def main():
    checks={};metas=[]
    for path in (ROOT/'data'/'fluid').glob('*.npz'):
        meta=json.loads(path.with_suffix('.json').read_text(encoding='utf-8'));z=np.load(path);metas.append(meta)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256']:raise AssertionError('Corrupt '+str(path))
        e=z['energy'];e0=e[0,0];t=z['time']
        checks[path.stem+'_unforced_energy_decay']=bool(np.max(np.diff(e[:,1]))/e0<1e-9)
        checks[path.stem+'_after_ramp_energy_decay']=bool(np.max(np.diff(e[t>=2-1e-9,2]))/e0<1e-9)
    c=P['fluid']['gates']
    checks['fluid_energy_budget']=max(m['energy_budget_residual'] for m in metas)<c['normalized_budget']
    checks['fluid_preparation_energy_budget']=max(m['burn_budget_residual'] for m in metas)<c['normalized_budget']
    checks['fluid_divergence']=max(m['divergence_max'] for m in metas)<c['divergence']
    checks['fluid_cfl']=max(m['max_cfl'] for m in metas)<c['cfl']
    checks['fluid_clones']=max(m['clone_error'] for m in metas)==0
    r=load('recovery_summary.json');nv=load('network_validation.json');rec=load('recurrence.json')
    checks['fluid_refinements']=len(r['validation']['refinement'])==8 and all(q['pass'] for q in r['validation']['refinement'])
    checks['network_refinements']=len(nv['refinement'])==12 and all(q['pass'] for q in nv['refinement'])
    checks['recurrence_replay_matches']=all(q['replay_energy_error']<1e-12 for q in rec['records'])
    checks['learning_parameter_control']=load('learned_network_parameters.json')['max_abs_coefficient_error']<.003
    checks['network_oracle_control']=load('network_forecast.json')['oracle']['mean']<1e-7
    checks['network_positive_bounded_states']=True
    for path in (ROOT/'data'/'network').glob('*.npz'):
        z=np.load(path)
        for k in ('base','perturbed'):
            # The kick can briefly exceed one; compare against its own initial maximum.
            yy=z[k];x=np.exp(yy);bound=max(1.,float(x[0].max()))
            checks['network_positive_bounded_states'] &= bool(np.isfinite(yy).all() and (x>0).all() and np.max(x)<=bound+1e-8)
    checks['upstream_operators_unchanged']=hashlib.sha256((ROOT/'vendor'/'numerics.py').read_bytes()).hexdigest()=='d95ed8e4cb9895cf7247f4f30722479ab60055af8625d94e8396a124c5aab502'
    tests=ET.parse(ROOT/'results'/'tests.xml').getroot()
    cases=tests.findall('.//testcase');failures=tests.findall('.//failure')+tests.findall('.//error')
    checks['automated_tests']=len(cases)>=14 and len(failures)==0
    summary={'all_passed':all(checks.values()),'checks':checks,'automated_tests_passed':len(cases),
             'max_fluid_budget':max(m['energy_budget_residual'] for m in metas),
             'environment':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},
             'limits':'Grid checks cover four fluid seeds; no measurement-noise or turbulence generality is claimed.'}
    save(ROOT/'results'/'verification.json',summary)
    entries=[]
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or any(x in path.parts for x in ('__pycache__','.pytest_cache')) or path.name in ('manifest.json','recordings.zip'):continue
        entries.append({'path':path.relative_to(ROOT).as_posix(),'size':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    save(ROOT/'manifest.json',{'protocol_sha256':hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),'files':entries})
    if not summary['all_passed']:raise AssertionError([k for k,v in checks.items() if not v])
    print(f'PASS: {len(checks)} gates, {len(cases)} automated tests, {len(entries)} hashed files',flush=True)

if __name__=='__main__':main()
