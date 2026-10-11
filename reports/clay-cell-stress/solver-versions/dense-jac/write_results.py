"""Populate the human report exclusively from recorded, audited results."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def main():
    a=json.loads((DATA/'analysis.json').read_text());rows=json.loads((DATA/'all-runs.json').read_text());r={v['name']:v for v in rows};aud={v['name']:v for v in a['audits']}
    text=f"{a['completed']} of {a['configurations']} final-protocol attempts completed their requested durations. "
    text+=f"All completed geometry/energy audits passed: **{a['all_completed_audits_passed']}**. Maximum completed scaled energy-budget error is {a['maximum_completed_budget_error']:.7g} (target 0.001). "
    text+='A stopped computation is listed below and is excluded from completed comparisons. Timings depend on the machine.\n\n'
    old=a['prior_baseline'][0];new=r.get('prior-extreme-new-one')
    if new:
        text+='The same previous extreme one-way forcing (alpha 16, rho 160, viscosity 0.15, capacity 0.01, initial Lorenz [-8,8,27], tilted ellipse, time 1) gives:\n\n'
        text+='| Measurement | Earlier material | Reinforced material |\n| --- | ---: | ---: |\n'
        for title,previous,current in [('Smallest swept cell area / reference',old['minimum_swept_area_ratio'],new['audit']['minimum_swept_area_ratio']),('Smallest saved shape quality / reference',old['min_quality_ratio'],new['metrics']['min_quality_ratio']),('Smallest saved altitude / reference',old['min_altitude_ratio'],new['metrics']['min_altitude_ratio'])]:text+=f'| {title} | {previous:.7g} | {current:.7g} |\n'
        text+='\nThe shape and altitude comparison uses saved frames on the same 0.04 recording grid. Swept area is the minimum certified across all accepted straight steps. These are different measures and are not silently interchanged.\n\n'
    text+='![Cell measurements over model time](figures/cell-quality.png)\n\n![Old and revised meshes at the same forcing and saved time](figures/molding-comparison.png)\n\n'
    text+='The combined stress cases double and quadruple the former maximum amplitude, increase rho to 256 and 300, and reduce viscosity to 0.075 and 0.0375. Each pair has matched parameters and starting state.\n\n'
    text+='| Pair | Requested duration | Reached one-way / two-way | Smallest swept area, one-way / two-way |\n| --- | ---: | --- | --- |\n'
    for label in ['prior-extreme-new','harder32','harder64']:
        names=[label+'-'+d for d in ['one','two']]
        if all(n in r for n in names):
            x,y=[r[n] for n in names];text+=f"| {label} | {x['requested_end']:g} | {x['completed_end']:.6g} / {y['completed_end']:.6g} | {x['audit']['minimum_swept_area_ratio']:.6g} / {y['audit']['minimum_swept_area_ratio']:.6g} |\n"
    text+='\nDefault barriers enforce lower domains at area ratio 0.25 and relative shape quality 0.35. Observed margins above those floors are evidence from these finite calculations; the floors themselves are chosen constraints. Compression still occurs.\n\n'
    text+='| Completed pair / second-half window | Material motion change with feedback | Lorenz Z fluctuation change | Final surface gap, one-way / two-way |\n| --- | ---: | ---: | --- |\n'
    for pair in a['pairs']:
        if pair['key'].startswith(('old-method','zero-port')) or '-fine' in pair['key']:continue
        c=pair['changes_percent'];fmt=lambda v:'undefined (zero reference)' if v is None else f'{v:+.4f}%'
        text+=f"| {pair['key']} / {pair['start']:g}–{pair['end']:g} | {fmt(c['material_motion_rms'])} | {fmt(c['z_std'])} | {pair['one_gap']:.5g} / {pair['two_gap']:.5g} |\n"
    text+='\n![Recorded reciprocal effects](figures/reciprocal.png)\n\nThese comparisons measure the change when reaction is enabled. They do not mean that stronger forcing always improves contact or that the two systems received identical force histories after feedback began. Larger surface gaps mean looser surrounding contact.\n\n'
    if a['long_blocks']:
        for key in ['material_motion_rms','z_std']:
            values=[b['changes_percent'][key] for b in a['long_blocks']];text+=f"Across {len(values)} completed six-unit blocks, {key.replace('_',' ')} changes range from {min(values):+.5f}% to {max(values):+.5f}%. "
        text+='These are observed ranges across three starts and time blocks, not confidence intervals or independent identically distributed samples.\n\n![Sustained finite-window effects](figures/long-windows.png)\n\n'
    checks=a['solver_checks'];text+=f"{len(checks)} solver/control comparisons are recorded; all pass their stated limits: **{a['all_solver_checks_passed']}**. "
    text+=f"Largest scaled state discrepancy is {max((c['error'] for c in checks),default=0):.7g}. The usual short-state target is 0.0001; the zero-port control uses 1e-12. See the exact durations and references in `data/analysis.json`.\n\n"
    if a['work_refinement']:
        maxwork=max(v for c in a['work_refinement'] for v in c['scaled_differences'].values());text+=f"The largest scaled discrepancy between base and tighter BDF gross forward/return-work quadrature is {maxwork:.7g}. Both signs of work are integrated by accepted-step trapezoids and are approximate; net work has its own ODE state.\n\n"
    text+='| Steady forming mesh / time 6 | Final imprint | Final surface gap RMS |\n| --- | ---: | ---: |\n'
    for n in [24,48,72]:
        name='steady-mesh'+str(n)
        if name in r and not r[name]['failure']:
            f=r[name]['metrics']['final'];text+=f"| {n} angles | {f['imprint']:.7g} | {f['surface_gap_rms']:.7g} |\n"
    for kind in ['True','False']:
        name='form-'+kind;release='release-'+kind
        if name in r and release in r and not r[name]['failure'] and not r[release]['failure']:
            initial=r[name]['metrics']['final']['imprint'];end=r[release]['metrics']['final']['imprint'];text+=f"\n{'Yielding' if kind=='True' else 'Elastic'} forming/release control: imprint after forming {initial:.7g}, after twelve release units {end:.7g}. "
            text+='Release begins from the complete saved forming state, including its Lorenz coordinates, rather than restarting from the named seed.\n'
    sensitivity=[r[n] for n in ['strength25','strength400','floor0.15','floor0.35'] if n in r]
    if sensitivity:
        text+='\n| Changed material choice / same prior extreme | Smallest swept area ratio | Smallest swept quality ratio |\n| --- | ---: | ---: |\n'
        for v in sensitivity:text+=f"| {v['name']} | {v['audit']['minimum_swept_area_ratio']:.6g} | {v['audit']['minimum_swept_quality_ratio']:.6g} |\n"
    if a['failures']:
        text+='\n**Incomplete attempts:**\n\n'
        for v in a['failures']:text+=f"- {v['name']}: reached {v['completed_end']:.7g} of {v['requested_end']:g}; {v['failure']['message']}.\n"
    else:text+='\nNo final-protocol attempt has an integration failure. Initial pilot budget stops and the corrected control-test setup error are retained under `initial-attempts`; they are not counted as successful final experiments.\n'
    text+='\nThe raw states, all-run ledger, per-run audits and source/file hashes accompany this report. These are numerical tests of an abstract constitutive choice; they do not validate physical clay against published measurements.\n'
    path=ROOT/'README.md';s=path.read_text(encoding='utf-8');s=re.sub(r'<!-- RESULTS -->.*?<!-- /RESULTS -->','<!-- RESULTS -->\n'+text+'<!-- /RESULTS -->',s,flags=re.S);path.write_text(s,encoding='utf-8')
if __name__=='__main__':main()
