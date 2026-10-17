"""Preserve original outcomes, then apply a numerically motivated grid refinement.

Grid selection was frozen before fresh memory effects were inspected. Physical
source files and the original execution matrix remain unchanged.
"""
import json,copy
import analysis as an
import followup_analysis as follow
import production                 # Bind original matrix settings before overrides.
import verify_extension as verify
from final_experiment import ROOT,dump,digest

def apply_probe_stability(final):
    """A converged response also requires a converged common intervention."""
    safeguard=json.loads((ROOT/'common-probe-grid-safeguard.json').read_text())
    threshold=safeguard['pulse_relative_limit']
    common=json.loads((ROOT/'results'/'common-analysis.json').read_text())
    off=json.loads((ROOT/'results'/'unforced-analysis.json').read_text())
    passed={}
    for regime in an.effective_regimes():
        rows=[r for r in common['grid'] if r['regime']==regime]
        passed[regime]=bool(rows) and all(r['pulse_relative']<threshold for r in rows)
        final['checks']['common/'+regime]['probe_grid']=passed[regime]
    for regime,result in off['regimes'].items():
        rows=[r for r in result['grid'] if r['method']=='common']
        final['checks']['unforced/'+regime]['common_probe_grid']=bool(rows) and all(r['pulse_relative']<threshold for r in rows)
    for key,claim in common['claims'].items():
        pulse_pass=passed[key.split('/')[0]]
        claim['probe_grid_passed']=pulse_pass
        claim['numerically_qualified']=claim['numerically_qualified'] and pulse_pass
        claim['physical_claim_passed']=claim['physical_claim_passed'] and pulse_pass
        for value in final['memory_assessment']['common_probe_cells'][key].values():
            value['probe_grid_passed']=pulse_pass
            value['numerically_qualified']=value['numerically_qualified'] and pulse_pass
            value['qualified_positive']=value['qualified_positive'] and pulse_pass
    strict=final['memory_assessment']['common_probe_cells']['smooth/coincident']
    for answer,metric in [('identical_common_probe_history_response_qualified','retained'),
                          ('identical_common_probe_positive_attenuation_qualified','attenuation'),
                          ('identical_common_probe_scrambling_specificity_qualified','specificity')]:
        final['memory_assessment'][answer]=strict[metric]['qualified_positive']
    final['common_claims']=common['claims']
    final['probe_grid_stability']=dict(safeguard_sha256=digest(ROOT/'common-probe-grid-safeguard.json'),
        pulse_relative_limit=threshold,regime_passed=passed,scope=safeguard['scope'])
    dump(ROOT/'results'/'common-analysis.json',common)

def main():
    freeze=json.loads((ROOT/'refinement-analysis-freeze.json').read_text())
    for name,expected in freeze['sources'].items():
        if digest(ROOT/name)!=expected:raise RuntimeError('Analysis/amendment source changed after freeze: '+name)
    observation=json.loads((ROOT/'results'/'refined-observation-validation.json').read_text())
    if not observation['all_passed']:raise RuntimeError('Refined CUDA/direct observation agreement failed')
    if observation['source_raw_sha256']!=digest(an.filename(82000,'forced_nu002','coincident',128,.002)):raise RuntimeError('Refined observation input fingerprint mismatch')
    if any(digest(ROOT/f)!=h for f,h in observation['sources'].items()):raise RuntimeError('Refined observation source fingerprint mismatch')
    streaming=json.loads((ROOT/'results'/'bounded-streaming-validation.json').read_text())
    if not streaming['passed'] or any(digest(ROOT/f)!=h for f,h in streaming['sources'].items()):raise RuntimeError('Packed streaming validation/provenance failed')
    for name in ['gpu-pipeline-complete.json','cpu-pipeline-complete.json','refinement-complete.json']:
        if not (ROOT/'results'/name).exists():raise RuntimeError('Incomplete study: '+name)
    plan=json.loads((ROOT/'fresh-grid-refinement.json').read_text())
    done=json.loads((ROOT/'results'/'refinement-complete.json').read_text())
    if done['planned']!=len(plan['tasks']) or len(done['completed'])!=len(plan['tasks']):raise RuntimeError('Incomplete refinement matrix')
    original=follow.all_results()
    for name in ['main-analysis.json','common-analysis.json','unforced-analysis.json','controls-analysis.json','extension-summary.json']:
        dump(ROOT/'results'/('baseline-'+name),json.loads((ROOT/'results'/name).read_text()))
    original_controls=json.loads((ROOT/'results'/'controls-analysis.json').read_text())
    final_controls=copy.deepcopy(original_controls)
    refined=[]
    for f in sorted((ROOT/'independent-refined-data').glob('*.npz')):
        _,m=an.read(f,['modes']);refined.append(dict(folder='independent-refined-data',file=f.name,**m))
    if len(refined)!=plan.get('additional_independent_comparisons',0):raise RuntimeError('Incomplete refined independent solver checks')
    final_controls['independent_solver'].extend(refined)
    old=an.effective_regimes;cfg=copy.deepcopy(old())
    for regime,changes in plan['analysis_overrides'].items():cfg[regime].update(changes)
    selected=lambda:copy.deepcopy(cfg)
    an.effective_regimes=selected;follow.effective_regimes=selected
    def selected_controls():
        result=copy.deepcopy(final_controls);dump(ROOT/'results'/'controls-analysis.json',result);return result
    old_controls=follow.controls;follow.controls=selected_controls
    final=follow.all_results()
    apply_probe_stability(final)
    final['checks']['refined_observation']={'direct_cpu_cuda_agreement':observation['all_passed']}
    final['checks']['bounded_streaming']={'input_and_coupled_step_equivalence':streaming['passed']}
    final['all_numerical_checks_passed']=all(all(group.values()) for group in final['checks'].values())
    final['memory_assessment']['completed_followups']={
        'equal_work_probes':'Three methods, including an identical common impulse with matched work and RMS strength within each comparison',
        'preparation_probe_locations':'Three paired layouts; identical common impulse repeated at all three for smooth flows',
        'wider_flow_conditions':'Two finite-time forced viscosities; four paired fresh flow seeds, without a stationarity claim',
        'independent_solver':'Fourteen full-window comparisons at normal/tight tolerances; separate implementation by the same author',
        'structure_intervention':'Four phase scrambles, displacement-matched intact translation and identical-state reset',
        'spatial_extent':'Endpoint radial footprint of velocity-response contrast; descriptive, not persistent collective-memory evidence'}
    final['grid_refinement']=dict(plan_sha256=digest(ROOT/'fresh-grid-refinement.json'),extra_configurations=len(plan['tasks']),
        analysis_controller_sha256=digest(ROOT/'final_analysis.py'),
        analysis_regimes=cfg,original_all_numerical_checks_passed=original['all_numerical_checks_passed'],
        original_failed_checks={group:[key for key,passed in value.items() if not passed] for group,value in original['checks'].items() if not all(value.values())},
        selection_basis='Numerical grid agreement, before inspecting fresh memory effects; all regime flows/layouts refined; thresholds unchanged')
    dump(ROOT/'results'/'extension-summary.json',final)
    import report_extension
    report_extension.report()
    # The original 212-configuration execution plan is audited separately from
    # the additional refinement, not silently replaced with new grid counts.
    an.effective_regimes=old;follow.effective_regimes=old;follow.controls=old_controls
    verify.verify()

if __name__=='__main__':main()
