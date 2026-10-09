"""Verify saved numerical data and delivery integrity without rerunning the suite."""
import json
import csv
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from core import ROOT,sha,save_json,PilotFlow,energy,wrap


def main():
    records=[];maxerr=0.;maxmarker=0.;clones=0
    for path in (ROOT/'data').rglob('*.json'):
        m=json.loads(path.read_text())
        if 'diagnostics' not in m:continue
        file=path.with_suffix('.npz');assert sha(file)==m['sha256']
        for name,expected in m['sources'].items():
            src=ROOT/('vendor/numerics.py' if name=='numerics.py' else name)
            assert sha(src)==expected,(path,name,'Source changed')
        z=np.load(file);assert all(np.isfinite(z[k]).all() for k in z.files)
        c=m['config'];f=PilotFlow(c['n'],c['nu'],c['nonlinear'],c['forcing']);last=m['diagnostics'][-1]
        for key,h in [('baseline_energy',z['final_hb']),('disturbed_energy',z['final_hp']),('perturbation_energy',z['final_hp']-z['final_hb'])]:
            err=abs(energy(f,h)-last[key]);assert err<1e-12;maxerr=max(maxerr,err)
        marker=float(np.sqrt(np.mean(np.sum(wrap(z['positions_perturbed'][-1]-z['positions_base'][-1])**2,axis=1))))
        assert abs(marker-last['marker_difference_rms'])<1e-14;maxmarker=max(maxmarker,abs(marker-last['marker_difference_rms']))
        ix=np.argmin(abs(z['time']-c['kick_time']))
        assert np.array_equal(z['positions_base'][:ix+1],z['positions_perturbed'][:ix+1])
        if c['amplitude']==0:
            assert np.array_equal(z['final_hb'],z['final_hp']);assert np.array_equal(z['positions_base'],z['positions_perturbed']);clones+=1
        assert m['max_cfl']<.4 and m['diffusion_number']<.5
        maxbudget=max(abs(r[k]) for r in m['diagnostics'] for k in ['baseline_energy_residual','disturbed_energy_residual','perturbation_budget_residual'])
        assert maxbudget<1e-4,(path,'Energy gate failed')
        records.append({'file':str(file.relative_to(ROOT)).replace('\\','/'),'sha256':sha(file),'seed':m['seed'],'max_budget_residual':maxbudget})
    assert len(records)==93,len(records)
    particles=[]
    for path in (ROOT/'data/particles').glob('*.json'):
        m=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));assert all(np.isfinite(z[k]).all() for k in z.files)
        assert m['gas_reproduction_max_error']<1e-10
        assert max(m['max_particle_Re_by_species'])<.1
        source=list((ROOT/'data').rglob(m['source']));assert any(sha(q)==m['source_sha256'] for q in source)
        assert m['source_code_sha256']==sha(ROOT/'particles.py')
        particles.append({'file':path.name,'max_particle_Re':max(m['max_particle_Re_by_species'])})
    assert len(particles)==7
    tests=ET.parse(ROOT/'results/tests.xml').getroot();suites=list(tests.iter('testsuite'))
    ntests=sum(int(x.attrib['tests']) for x in suites);fails=sum(int(x.attrib.get('failures',0))+int(x.attrib.get('errors',0)) for x in suites)
    assert ntests>=22 and fails==0
    with (ROOT/'results/numerical_controls.csv').open() as f:
        sensitivity=list(csv.DictReader(f))
    failures=[{'regime':r['regime'],'seed':int(r['seed']),'comparison':r['comparison'],'marker_rms_change':float(r['marker_rms_change']),'limit':.005}
              for r in sensitivity if float(r['marker_rms_change'])>.005]
    stress=json.loads((ROOT/'results/stress_extension.json').read_text())
    oscillator_checks={}
    for name in ['phase_locking','dynamics','pendulum','relaxation','weak_nonlinear','parametric','bifurcations']:
        z=np.load(ROOT/'results'/(name+'.npz'))
        assert all(np.isfinite(z[k]).all() for k in z.files),name
        oscillator_checks[name]='saved arrays finite'
    pen=json.loads((ROOT/'results/pendulum.json').read_text())
    assert max(r['relative_period_error'] for r in pen['records'])<1e-6
    assert max(r['max_absolute_energy_error'] for r in pen['records'])<1e-8
    weak=json.loads((ROOT/'results/weak_nonlinear.json').read_text())
    assert max(abs(r['measured_period']/r['integral_period']-1) for r in weak['duffing'])<1e-7
    assert max(r['max_energy_drift'] for r in weak['duffing'])<1e-8
    relaxation=json.loads((ROOT/'results/relaxation.json').read_text())
    assert max(r['max_energy_balance_absolute_residual'] for r in relaxation['records'])<1e-5
    parametric=json.loads((ROOT/'results/parametric.json').read_text())
    assert max(r['max_energy_work_residual'] for r in parametric['trajectory_controls'])<1e-8
    assert max(abs(r['determinant']-1) for r in parametric['floquet'])<1e-8
    bifurcations=json.loads((ROOT/'results/bifurcations.json').read_text())
    assert max(r['relative_error'] for r in bifurcations['bottleneck'])<1e-7
    assert max(r['max_radial_solution_error'] for r in bifurcations['hopf'])<1e-7
    output={'status':'passed','fluid_configurations':len(records),'particle_configurations':len(particles),
            'test_count':ntests,'test_failures':fails,'identical_no_kick_configs':clones,
            'max_recomputed_energy_error':maxerr,'max_recomputed_marker_rms_error':maxmarker,
            'fluid_records':records,'particle_records':particles,
            'numerical_sensitivity_failures':failures,
            'stress_followup':stress,
            'oscillator_checks':oscillator_checks,
            'scope':'Passed data integrity, automated tests and energy/CFL gates. Numerical sensitivity failures are retained explicitly; this is not a proof of continuum convergence or unrestricted physical validity.'}
    save_json(ROOT/'results/verification.json',output)
    files={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='manifest.json' and '__pycache__' not in str(p) and '.pytest_cache' not in str(p)}
    save_json(ROOT/'manifest.json',{'algorithm':'SHA-256','files':files})
    print(json.dumps({k:v for k,v in output.items() if k not in ['fluid_records','particle_records']},indent=2))


if __name__=='__main__':main()
