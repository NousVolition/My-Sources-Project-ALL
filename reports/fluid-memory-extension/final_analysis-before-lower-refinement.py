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

def main():
    for name in ['gpu-pipeline-complete.json','cpu-pipeline-complete.json','refinement-complete.json']:
        if not (ROOT/'results'/name).exists():raise RuntimeError('Incomplete study: '+name)
    plan=json.loads((ROOT/'fresh-grid-refinement.json').read_text())
    done=json.loads((ROOT/'results'/'refinement-complete.json').read_text())
    if done['planned']!=len(plan['tasks']) or len(done['completed'])!=len(plan['tasks']):raise RuntimeError('Incomplete refinement matrix')
    original=follow.all_results()
    for name in ['main-analysis.json','common-analysis.json','unforced-analysis.json','controls-analysis.json','extension-summary.json']:
        dump(ROOT/'results'/('baseline-'+name),json.loads((ROOT/'results'/name).read_text()))
    original_controls=json.loads((ROOT/'results'/'controls-analysis.json').read_text())
    old=an.effective_regimes;cfg=copy.deepcopy(old())
    for regime,changes in plan['analysis_overrides'].items():cfg[regime].update(changes)
    selected=lambda:copy.deepcopy(cfg)
    an.effective_regimes=selected;follow.effective_regimes=selected
    old_controls=follow.controls;follow.controls=lambda:copy.deepcopy(original_controls)
    final=follow.all_results()
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
