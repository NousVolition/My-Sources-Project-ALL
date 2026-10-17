"""Expand resolution controls uniformly before inspecting fresh memory effects.

This one-time amendment preserves original files and execution plans. It changes
only additional run scheduling and the final choice of grids for inference.
"""
import json,datetime,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    file=ROOT/'fresh-grid-refinement.json';plan=json.loads(file.read_text())
    if plan.get('additional_independent_comparisons')==2:
        assert len(plan['tasks'])==54
        print('Lower-viscosity refinement amendment already recorded');return
    assert len(plan['tasks'])==25
    failure=json.loads((ROOT/'results'/'lower-viscosity-first-grid.json').read_text())
    for name in ['fresh-grid-refinement.json','refinement-analysis-freeze.json','refinement_gpu_pipeline.py','report_extension.py','spatial_assessment.py','parallel_cpu_pipeline.py','final_analysis.py','verify_extension.py']:
        p=ROOT/name;backup=p.with_name(p.stem+'-before-lower-refinement'+p.suffix)
        if not backup.exists():backup.write_bytes(p.read_bytes())
    # Each evolution is independent. Smaller streaming batches change resource
    # use only. Original helpers and their source fingerprints remain intact.
    for old,new in [('common_probe.py','safe_common_probe.py'),('unforced_readout.py','safe_unforced_readout.py'),('common_unforced.py','safe_common_unforced.py')]:
        text=(ROOT/old).read_text().replace("'"+old+"'","'"+new+"'")
        if old=='common_probe.py':
            text=text.replace('for off in range(0,len(states),8):','batch=4 if n>=128 else 8\n    xp.get_default_memory_pool().free_all_blocks()\n    for off in range(0,len(states),batch):').replace('states[off:off+8]','states[off:off+batch]')
        elif old=='unforced_readout.py':
            text=text.replace("for off in range(0,len(a['start_fields']),8):","batch=4 if n>=128 else 8\n    for off in range(0,len(a['start_fields']),batch):").replace("a['start_fields'][off:off+8]","a['start_fields'][off:off+batch]")
        else:
            text=text.replace('for off in range(0,len(starts),16):','batch=4 if n>=128 else 16\n    for off in range(0,len(starts),batch):').replace('starts[off:off+16]','starts[off:off+batch]')
        (ROOT/new).write_text(text,encoding='utf-8')
    text=(ROOT/'independent_check_full.py').read_text().replace("n=32 if regime=='smooth' else max(64,min(effective_regimes()[regime]['grids']));horizon=.4","n=96;horizon=.4").replace("'independent-full-data'","'independent-refined-data'").replace("'independent_check_full.py'","'refined_independent_check.py'")
    (ROOT/'refined_independent_check.py').write_text(text,encoding='utf-8')
    plan['analysis_overrides']['forced_nu002']={'grids':[96,128]}
    for seed in range(82000,82004):
        for layout in ['coincident','writer_moved','probe_moved']:
            plan['tasks'].append(dict(kind='main',seed=seed,regime='forced_nu002',layout=layout,n=128,dt=.002))
            if seed<82002:plan['tasks'].append(dict(kind='main',seed=seed,regime='forced_nu002',layout=layout,n=128,dt=.001))
    for seed in range(82000,82004):
        plan['tasks'].append(dict(kind='common',seed=seed,regime='forced_nu002',layout='coincident',n=128,dt=.002))
        if seed==82000:plan['tasks'].append(dict(kind='common',seed=seed,regime='forced_nu002',layout='coincident',n=128,dt=.001))
    for kind in ['unforced','common_unforced']:
        for seed in range(82000,82002):
            plan['tasks'].append(dict(kind=kind,seed=seed,regime='forced_nu002',n=128,dt=.002))
            if seed==82000:plan['tasks'].append(dict(kind=kind,seed=seed,regime='forced_nu002',n=128,dt=.001))
    assert len(plan['tasks'])==54
    plan.update(extra_configurations=54,additional_independent_comparisons=2,
        lower_viscosity_failure=failure,lower_viscosity_amendment_time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        lower_viscosity_reason='The first fresh 80/96 pair fails the unchanged FTLE threshold. Refine the entire lower-viscosity regime to 96/128, including designated half steps and all common/off controls; add normal/tight independent comparisons at 96. No fresh memory effects inspected.',
        independent_refinement=dict(seed=82000,regime='forced_nu002',n=96,horizon=.4,tolerances=['normal','tight'],scope='Selected full-window readouts; does not repeat preparation or all trajectories'))
    file.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    print('Recorded 54 additional GPU configurations and 2 additional independent comparisons.')

if __name__=='__main__':main()
