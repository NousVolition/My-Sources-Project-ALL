"""Validate, index and checksum the completed experiment, preserving prior work."""
from pathlib import Path
from html.parser import HTMLParser
import ast,hashlib,json,csv
import numpy as np
from water_battery import cases,save

HERE=Path(__file__).resolve().parent;D=HERE/'data'


def main():
    result=json.loads((D/'results.json').read_text());models=json.loads((D/'model_results.json').read_text());validation=json.loads((D/'integration_validation.json').read_text())
    assert result['complete'] and result['completed_main_runs']==156 and len(result['refinements'])==12
    assert models['checks']['passed'] and validation['passed'] and all(result['analysis_checks'].values())
    registry=[]
    for path in sorted(D.glob('*.npz')):
        if path.name=='model_trajectories.npz':continue
        with np.load(path) as data:
            config=json.loads(str(data['configuration_json']));digest=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
            assert digest==str(data['configuration_sha256']),path
            assert abs(float(data['time_ps'][-1])-config['case']['duration'])<1e-10
            for key in data.files:
                if data[key].dtype.kind in 'fci':assert np.isfinite(data[key]).all(),(path,key)
            registry.append(dict(file=path.name,kind=config['kind'],replica=config['replica'],case=config['case']['name'],duration_ps=float(data['time_ps'][-1]),configuration_sha256=digest,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    assert len(registry)==168 and abs(sum(r['duration_ps'] for r in registry)-6348)<1e-8
    save(D/'run_registry.json',registry)
    protocol=json.loads((D/'protocol.json').read_text());protocol.update(phase_rule='Window begins after max(one period, one-quarter run). |Pxy|>=0.15 on >=95% frames, longest valid segment>=80%, concentration>=0.9, relative drift<=0.05, phase excursion<pi. Threshold sensitivity at 0.1 and 0.2.',
        source_module_sha256={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE.parent/'experiment.py',HERE.parent/'driven-water-extension'/'drive_water.py',HERE.parent/'response-extension'/'impulse_response.py']},
        precision='OpenCL double for every molecular production/refinement run. Separate hardware benchmarks did not change the production engine.',
        nonlinear_fit_note='The unrestricted GLV mapping and a constrained-competition-plus-source variant are exploratory model comparisons; the second was added after observing the first fit. Assessment trajectories were excluded from coefficient fitting, not from all development decisions.',
        main_completed=156,refinements_completed=12,total_recorded_molecular_ps=6348)
    save(D/'protocol.json',protocol)
    rows=result['rows'];ledger=[]
    for case in cases():
        selected=[r for r in rows if r['case']==case['name']]
        ledger.append(dict(family=case['name'],domain='molecular simulation and ghost-dipole control',runs=len(selected),status='measured; all three configurations retained',tracking_screen_passes_water=sum(r.get('phase',{}).get('finite_window_tracking',False) for r in selected if r['kind']=='water') if case['mode']=='rotate' else None))
    rec=result['recurrence']
    ledger.extend([
        dict(family='unrestricted GLV fit to squared polarization',domain='reduced model fitted to water',status='does not outperform persistence in these two excluded-from-fit trajectories' if all(r['GLV_MSE']>r['persistence_MSE'] for r in rec['GLV_forecasts']) else 'mixed or improved relative to persistence',ratios=[r['GLV_MSE']/r['persistence_MSE'] for r in rec['GLV_forecasts']],negative_competition_coefficients=rec['unrestricted_GLV_has_negative_competition_coefficients']),
        dict(family='nonnegative competition plus constant source',domain='exploratory reduced model fitted to water',status='measured; model selection is exploratory',ratios=[r['MSE']/r['persistence_MSE'] if r['MSE'] is not None else None for r in rec['constrained_activity_forecasts']]),
        dict(family='pronounced chunks in chosen nine-state network',domain='engineered reference ODE',status='not observed by >80% late group-dominance criterion',fractions=[r['late_fraction_with_one_chunk_above_80pct'] for r in models['engineered_hierarchy']['rows']]),
        dict(family='heteroclinic network in molecular water',domain='hypothesis',status='not established; clustering and phase following do not identify invariant saddles and connecting manifolds'),
        dict(family='numerical convergence and phase-analysis fixtures',domain='validation',status='passed; finite measured error retained')])
    save(D/'outcome_ledger.json',ledger)
    save(D/'development_notes.json',dict(items=[
       'A protein-model equilibrium exactly on a search-grid point was initially missed; exact zeros and sign-change roots are now both retained. This affected the reference root finder, not water trajectories.',
       'An initial deterministic audit attempted a 0.0025 ps field update with 0.001 ps MD steps and rounded the count. An integer-step guard and a 0.0005 ps audit step corrected it. Production cases already used integer step ratios.',
       'The first serial water worker was interrupted for scheduling and resumed with disjoint replica workers. Completed archives were kept and their configuration hashes checked; incomplete integration was rerun.',
       'The nine-state parameters were not retuned to force strong chunks after the weak-group result. Unrestricted and constrained reduced-model outcomes are both reported.',
       'A SciPy roundoff warning occurred at a stationary reference trajectory; the returned solutions and invariants were checked. No production molecular integration failure was recorded.'
    ]))
    older=['','response-extension','damping-extension','switching-extension','stability-extension','oscillator-extension','driven-water-extension'];old_count=0
    for name in older:
        folder=HERE.parent/name
        for relative,expected in json.loads((folder/'manifest_sha256.json').read_text()).items():
            assert hashlib.sha256((folder/relative).read_bytes()).hexdigest()==expected,str(folder/relative);old_count+=1
    class Check(HTMLParser):
        def __init__(self):super().__init__();self.images=0
        def handle_starttag(self,tag,attrs):
            if tag=='img':
                d=dict(attrs);assert d.get('alt') and d.get('src','').startswith('data:image/png;base64,');self.images+=1
    parser=Check();parser.feed((HERE/'report.html').read_text(encoding='utf-8'));assert parser.images==9
    for p in HERE.glob('*.py'):ast.parse(p.read_text(encoding='utf-8'))
    files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ['manifest_sha256.json','partial_results.json'] and not p.name.endswith('.partial.npz')]
    save(HERE/'manifest_sha256.json',{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)})
    print(json.dumps(dict(main_runs=156,refinement_runs=12,total_ps=6348,reference_model_families=12,figures=9,unchanged_earlier_files=old_count,release_files=len(files)+1,release_bytes=sum(p.stat().st_size for p in files)+(HERE/'manifest_sha256.json').stat().st_size)),flush=True)


if __name__=='__main__':main()
