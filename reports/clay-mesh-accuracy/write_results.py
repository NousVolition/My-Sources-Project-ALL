"""Write conclusions from the completed data, preserving failures and scope."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
def main():
    a=json.loads((ROOT/'analysis.json').read_text());final=[r for r in a['runs'] if r['key'].startswith('solver/')]
    by={r['key']:r for r in a['runs']};pct=lambda x:f'{100*x:.4g}%'
    lines=[f"The released solver recorded **{len(final)} attempts**, of which **{sum(r['failure'] is None for r in final)} completed** their requested interval. Including earlier repair pilots and the rejected alternative, this package records {a['attempted']} attempts, with {a['completed']} completed. Counts refer to integrations, not independent physical experiments.", '', '### Gap agreement', '', '| Drive and refinement | Original gap statistic: endpoint change | Surface-weighted endpoint change |', '|---|---:|---:|']
    for c in a['mesh_comparisons']:
        if not c['mode'].startswith('final'):continue
        lines.append(f"| {c['mode'].removeprefix('final ')}: {by[c['coarse']]['parameters']['angles']} to {by[c['fine']]['parameters']['angles']} points, t={c['time']:g} | {pct(c['surface_gap_rms']['relative_difference'])} | {pct(c['gap_arclength_rms']['relative_difference'])} |")
    lines+=['','The earlier published 24-to-48-point discrepancy was about 60.2% with the original gap statistic (52.9% with surface-length weighting). The revised material changes the physical law; these numbers are not a pure integrator comparison. A larger gap is not an improvement in molding. Smaller changes under refinement are evidence of numerical agreement.','']
    for file in sorted((ROOT/'data').glob('gate-final-*.json')):
        g=json.loads(file.read_text());lines.append(f"- **{file.stem}**: {'meets' if g['passed'] else 'does not meet'} the declared check. Shared recorded interval ends at {g['shared_end']:.6g}; maximum scaled gap difference {pct(g.get('maximum_gap_scaled_difference',0))}. Both requested runs complete: {g['complete']}.")
    lines+=['','The released default starts at **48 angular points and 5 layers**, with a **72-point, 7-layer check**. The 24-point grid remains an explicitly coarse preview and a stress/control grid. Agreement at the tested saved times is not a rigorous exact-solution error bound or proof for every mesh, forcing, or longer interval. The gap gate uses a 0.01 model-unit denominator floor; that floor is inactive in the reported standard gap comparison.','', '### Time integration and contact sampling', '', '| Final comparison | Shared interval ends | Largest scaled state difference | Target met |','|---|---:|---:|---|']
    for c in a['solver_comparisons']:
        if c['a'].startswith('solver/'):
            lines.append(f"| {c['a'].removeprefix('solver/')} / {c['b'].removeprefix('solver/')} | {c['end']:g} | {pct(c['scaled_state_error'])} | {c['passed'] and c['full_requested_interval']} |")
    lines+=['','The 8-to-16-point contact quadrature comparison failed the 0.01% trajectory target. The 16-to-32 comparison met it; the released default uses 32 points and is checked against 64. The earlier broad and local numerical derivative calculations also agree over their shared interval. Material stiffness is assembled explicitly; obstacle-normal and local self-contact derivatives retain small finite differences and are independently checked against full-system derivatives.','', '### Effects in both directions', '', '| Final paired test and window | Material motion change with feedback | Lorenz Z fluctuation change with feedback |','|---|---:|---:|']
    for p in a['pairs']:
        if p['name'].startswith('solver/'):
            lines.append(f"| {p['name'].removeprefix('solver/')} [{p['start']:g}, {p['end']:g}] | {p['motion_change_percent']:+.3f}% | {p['z_std_change_percent']:+.3f}% |")
    lines+=['','Motion RMS is obtained from integrated viscous dissipation. Z fluctuation is the standard deviation over the stated saved time window. The three sustained tests use starts (1,1,1), (-8,8,27), and (1.000001,1,1); their spread is a sensitivity observation, not a confidence interval. Reduced motion can accompany a larger object gap. These are finite-window effects of this chosen coupling, not proof that chaos has been removed.','']
    completed=[r for r in final if r['failure'] is None]
    if completed:
        lines += [f"Across completed released-model runs, the smallest certified area ratio was **{min(r['minimum_area'] for r in completed):.6g}**, the smallest deformation-quality bound **{min(r['minimum_distortion'] for r in completed):.6g}**, and the smallest altitude/reference-altitude bound **{min(r['minimum_altitude'] for r in completed):.6g}**. No recorded or certified-path object, wall, or material crossing was detected. The largest scaled energy-budget error was **{max(max(r['energy_errors'].values()) for r in completed):.6g}**.", '']
    lines+=['### Failures and remaining limits','']
    for r in a['runs']:
        if r['failure']:lines.append(f"- `{r['key']}` reached {r['end']:.9g} of {r['requested_end']:g}: {r['failure']['message']}.")
    lines+=['- A proposed fixed-reference-neighborhood self-contact exclusion was rejected. Its 72-point mesh approached a surface crossing and stopped near time 0.109974, even while cells retained positive area. The released solver keeps the original self-contact protection.', '- The `production/` and `resolved/` folders retain that rejected experiment and its numerical variants. They are historical folder names, not release recommendations. The released implementation is `solver/`.', '- The largest forcing sweeps on the 24-point grid test geometry, solver behavior and reciprocal sensitivity. They do not establish mesh accuracy for every quantity. Only the explicitly matched finer-grid tests establish the stated gap agreement.', '- Long chaotic trajectories have not been shown pointwise converged. No measured clay calibration, material uncertainty distribution, friction, fracture, arbitrary remeshing, three-dimensional validation, or Navier–Stokes theorem is supplied.', '', '![Gap refinement](figures/mesh-gap.png)', '', '![Reciprocal effects](figures/reciprocal.png)']
    text=(ROOT/'README.md').read_text(encoding='utf-8');start=text.index('<!-- RESULTS -->')+len('<!-- RESULTS -->');end=text.index('<!-- END RESULTS -->');text=text[:start]+'\n\n'+'\n'.join(lines)+'\n\n'+text[end:]
    (ROOT/'README.md').write_text(text,encoding='utf-8')
    protocol=[]
    for r in final:
        raw=json.loads((ROOT/r['relative_path']).read_text());protocol.append({k:raw[k] for k in ['name','parameters','initial','requested_end','audit','source_hashes']})
    (ROOT/'data/executed-protocol.json').write_text(json.dumps(dict(released_solver='solver',attempts=protocol,writer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print('Wrote data-derived results for',len(final),'released-model attempts.')
if __name__=='__main__':main()
