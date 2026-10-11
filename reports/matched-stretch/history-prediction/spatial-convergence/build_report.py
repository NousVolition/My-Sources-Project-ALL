"""Render the completed frozen-screen analysis without tuning its criteria."""
from pathlib import Path
import json
from collections import Counter
HERE = Path(__file__).resolve().parent


def main():
    r = json.loads((HERE/'results.json').read_text()); p = r['protocol']; grids = r['grids']
    verdict = ('All predeclared operational numerical-stability screens pass for these fixed-state, short-interval measurements.'
               if r['all_screens_pass'] else 'The strict spatial-convergence claim remains unestablished: some predeclared numerical-stability screens fail.')
    text = [f'# Spatial refinement of mirrored arrangements and motion direction\n\n**{verdict}**\n',
    'This study tests the same four reused initial conditions at finer grids, with a common saved fluid state and identical marker positions. Arrangement/direction effects and convergence of individual deformation outcomes are evaluated separately. A surviving aggregate effect does not certify local event locations.\n',
    '## Frozen design and scope\n',
    f'The [protocol](protocol.json) was frozen before the grid-80 pilot. Seeds 1024–1027 start from their verified grid-56, step-0.005 full Fourier checkpoints at time 0.3 from the [velocity-reversal study](../velocity-reversal/README.md). Exact Fourier mode transfer, with FFT normalization, supplies the same band-limited present velocity field at every resolution. Central positions, global mirror positions, neighbor identities and anchored reflected offsets are identical across grids. No initial-condition replay was needed.\n',
    f'The main grid ladder is {", ".join(str(n)+"³" for n in grids)}, with time step 0.0025, half the previous study\'s step. All four seeds use each grid. Seed 1024 additionally uses step 0.00125 at grids {", ".join(r["half_step_comparisons"])}. The initial scope is sixteen main pairs and two half-step pairs; the predeclared grid-224 extension adds four main pairs and one half-step pair only if initial screens fail. This completed record contains {len(r["runs"])} pairs / {2*len(r["runs"])} forward fluid evolutions, including the pilot once through cache reuse.\n',
    'Every branch uses the unchanged parent Navier–Stokes solver, Heun method and periodic trilinear velocity sampler with its exact piecewise spatial Jacobian. Domain side 6, viscosity 0.02, zero mean and zero force are preserved. Branches evolve for 0.2 model-time units; saved horizons are 0.05, 0.1 and 0.2, with primary horizon 0.2. Velocity negation is a forward physical intervention with positive viscosity, not exact backward evolution.\n',
    'The four cells are original/global-mirrored marker clouds in normal/reversed fluid motion. A global mirror moves starting locations relative to a fixed-frame fluid field. For anchored geometry, eight neighbor offsets around each unchanged center are reflected. The 4096 virtual neighbor probes are independent overlapping diagnostic clouds, not a physical packing. All 5120 probes remain passive. They cannot move the fluid\'s deformation field or alter its stress.\n',
    'The pointwise target is forward log-largest-tangent-singular-value / horizon (numerical FTLE). High events use the unchanged original training cutoff **1.0106081570657628**. Finite-cloud deformation is a regularized least-squares offset-map log-singular-value rate; the pair target is mean neighbor-center log-distance change / horizon. Those finite targets have no binary cutoff and are not interchangeable with pointwise FTLE. The cloud regularization remains 1e-10 times the Gram trace, with a floor.\n',
    '## Do arrangement and direction effects persist?\n',
    'Means below average whole independent runs, rather than pooling labels or treating grids as new seeds.\n',
    '| Grid | Original-cloud direction RMS difference | Original event disagreement | Original event Jaccard | Mirrored event disagreement | Mirrored event Jaccard |\n|---|---:|---:|---:|---:|---:|']
    for n in grids:
        a = r['summary'][str(n)]['global_comparisons']['direction_original']; b = r['summary'][str(n)]['global_comparisons']['direction_mirror']
        text.append(f'| {n} | {a["same_label_RMSE"]["mean"]:.5f} | {100*a["event_disagreement_fraction"]["mean"]:.2f}% | {100*a["event_Jaccard"]["mean"]:.2f}% | {100*b["event_disagreement_fraction"]["mean"]:.2f}% | {100*b["event_Jaccard"]["mean"]:.2f}% |')
    text += ['\nFor anchored finite clouds let op, om, mp and mm denote original/normal, original/reversed, mirrored/normal and mirrored/reversed. A=[(mp-op)+(mm-om)]/2, D=[(om-op)+(mm-mp)]/2, I=(mm-mp)-(om-op). The table averages each contrast\'s per-run RMS. I has different scaling; these numbers do not partition variance or establish which factor dominates.\n',
             '| Grid | Arrangement A, cloud | Direction D, cloud | Interaction I, cloud | Arrangement A, pair | Direction D, pair | Interaction I, pair |\n|---|---:|---:|---:|---:|---:|---:|']
    keys = ['arrangement_average', 'direction_average', 'interaction_difference_in_differences']
    for n in grids:
        e = r['summary'][str(n)]['effects']
        text.append('| '+str(n)+' | '+' | '.join(f'{e[f][k]["mean"]:.5f}' for f in ['cloud', 'pair'] for k in keys)+' |')
    text += ['\nDescriptive 2000-repeat bootstraps resample the four whole runs. At the finest grid, finite-cloud contrast intervals are:']
    for k in keys:
        e = r['summary'][str(grids[-1])]['effects']['cloud'][k]
        text.append(f'- {k}: [{e["descriptive_ci95"][0]:.5f}, {e["descriptive_ci95"][1]:.5f}].')
    text += ['\nThese are four reused seeds, not a fresh confirmatory population. RMS is nonnegative, so an interval above zero is not itself a null test. Geometry effects reflect how finite offsets sample nonuniform motion; they do not establish arrangement-induced stress, material memory or persistent-identity control. The consistent-label permutation, passive-probe independence, affine-map and whole-state mirror controls from the parent suite remain applicable. No history predictor was refitted or newly scored out of sample.\n',
             '## Numerical refinement and frozen screens\n',
             'The required finest-pair target relative RMS difference is at most 1% in every seed/cell/family, normalized by the finer target\'s RMS. Pointwise high-event Jaccard must be at least 0.95 in every seed/cell. Successive RMS differences must contract. Finest-grid time-step differences on seed 1024 must be at most one tenth of its adjacent spatial differences. Per-run anchored contrast RMS and global direction RMS must drift by at most 5% between the finest pair. These are operational finite-run screens, not a mathematical convergence proof.\n',
             '| Adjacent grids | Pointwise relative RMS range | Cloud relative RMS range | Pair relative RMS range | Minimum pointwise event Jaccard | Maximum shared-band Eulerian relative RMS |\n|---|---:|---:|---:|---:|---:|']
    for a, b in zip(grids[:-1], grids[1:]):
        records = [r['spatial_comparisons'][str(seed)][f'{a}-{b}'] for seed in p['seeds']]
        ranges = []
        for f in ['global', 'cloud', 'pair']:
            v = [d[f][c]['relative_RMSE']*100 for d in records for c in ['original_normal', 'original_reversed', 'mirror_normal', 'mirror_reversed']]
            ranges.append(f'{min(v):.3f}–{max(v):.3f}%')
        j = min(d['global'][c]['event_Jaccard'] for d in records for c in d['global'])
        e = max(v for d in records for v in d['shared_band_Eulerian_relative_RMS'].values())
        text.append(f'| {a}→{b} | '+ ' | '.join(ranges)+f' | {j:.5f} | {e:.3e} |')
    counts = Counter(x['kind'] for x in r['screens']); passed = Counter(x['kind'] for x in r['screens'] if x['pass'])
    text += ['\n| Screen | Passed / evaluated |\n|---|---:|']
    text += [f'| {k} | {passed[k]} / {v} |' for k, v in counts.items()]
    if r['failures']:
        text += ['\nEvery failed screen is retained below; tolerances were not relaxed after seeing results.\n', '| Screen and identifier | Value | Required |\n|---|---:|---:|']
        text += [f'| {x["kind"]}: {x["identifier"]} | {x["value"]:.6g} | {x["comparison"]} {x["limit"]:.6g} |' for x in r['failures']]
    else:
        text += ['\nAll screens pass. This establishes operational numerical stability at the chosen tolerances for the measured finite intervals and these four inputs. It does not certify all local extrema, population generalization, or a continuum theorem.\n']
    text += ['\nThe frozen pass/fail screens apply to the primary horizon 0.2 only. Shorter horizons are descriptive checks, with no newly calibrated event threshold. Their finest-pair differences are:\n',
             '| Horizon | Pointwise relative RMS range | Cloud relative RMS range | Pair relative RMS range |\n|---|---:|---:|---:|']
    for horizon, rec in r['secondary_horizon_spatial_comparisons'].items():
        ranges = []
        for family in ['global', 'cloud', 'pair']:
            values = [v[family][cell]['relative_RMSE']*100 for v in [rec[str(seed)][f'{grids[-2]}-{grids[-1]}'] for seed in p['seeds']] for cell in v[family]]
            ranges.append(f'{min(values):.3f}–{max(values):.3f}%')
        text.append(f'| {horizon} | '+' | '.join(ranges)+' |')
    text += ['\n## Temporal, interpolation and solver diagnostics\n', '| Main grid | Max CFL | Max absolute relative energy residual | Max spectral divergence | Max probe tangent-volume error |\n|---|---:|---:|---:|---:|']
    for n in grids:
        d = [v for rec in r['runs'].values() if rec['n'] == n and rec['dt'] == p['dt'] for v in rec['diagnostics'].values()]
        rows = [row for v in d for row in v['rows']]
        text.append(f'| {n} | {max(v["max_cfl"] for v in d):.5f} | {max(abs(v["energy_budget_error"]) for v in rows):.3e} | {max(v["spectral_divergence_max"] for v in rows):.3e} | {100*max(v["probe_tangent_volume_error_max"] for v in rows):.3f}% |')
    text += ['\nThe Eulerian comparison uses normalized endpoint Fourier coefficients on the shared grid-56 frequency band; it does not measure omitted fine-grid modes. High-band energy fractions are recorded separately. The inherited trilinear interpolant is not exactly divergence-free and its Jacobian jumps at grid-cell faces. Good Eulerian budgets cannot remove that marker error. Successive target-difference alignment is recorded, so shrinking norms are not automatically treated as a smooth signed error expansion. No unsupported Richardson extrapolation or formal order is inferred. This follows the distinction between decreasing grid differences and an established asymptotic error regime in the [NASA grid-convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).\n',
             '| Half-step grid, seed 1024 | Pointwise RMS difference range | Cloud RMS difference range | Pair RMS difference range | Minimum pointwise event Jaccard |\n|---|---:|---:|---:|---:|']
    for n, rec in r['half_step_comparisons'].items():
        cells = list(rec['global']); ranges = []
        for fam in ['global', 'cloud', 'pair']:
            v = [rec[fam][c]['same_label_RMSE'] for c in cells]; ranges.append(f'{min(v):.3e}–{max(v):.3e}')
        j = min(rec['global'][c]['event_Jaccard'] for c in cells)
        text.append(f'| {n} | '+' | '.join(ranges)+f' | {j:.5f} |')
    text += ['\n## Provenance, execution and reproduction\n',
             '[input-verification.json](input-verification.json) records the 27 reused archive/checkpoint/metadata/source files checked against GitHub before launching. New per-job metadata records exact trajectory and endpoint-spectrum hashes, the common reference checkpoint hash, transfer invariants, executed-source hashes and diagnostics. [execution.json](execution.json) records the initial matrix; any required extension has a separate execution record. [extension-decision.json](extension-decision.json) retains the initial screen failures. Publication hashes are in [publication-manifest.json](publication-manifest.json).\n',
             '[scheduling-audit.json](scheduling-audit.json) and per-job claims/logs document disjoint parallel scheduling. The first short-lived-shell launch never initialized numerical work; it created no claims, data or checkpoints. The persistent launcher was verified before retrying. Completed caches are verified by exact source/data hashes. Partial checkpoints are preserved and refused by the sequential controller before any integration; if it arrives before a helper completes, it waits through a documented guarded restart. Numerical sources were not modified while active. Recovery checkpoints are local restart evidence, not published bulk data.\n',
             'The [results](results.json), [all-horizon target arrays](targets.npz), [grid comparison table](grid-comparisons.csv), individual [recorded trajectories and shared-band spectra](recorded-data) and [plot](convergence.png) support the numerical tables. In grid/time comparison records, inherited helper keys named normal/reversed denote the first/second supplied array (coarser/finer grid or main/half step), rather than a change in motion direction. The complete parent and new transfer test suite passed **32 tests**, recorded in [tests.xml](tests.xml). Original long matrices and previously completed studies were not restarted.\n',
             '```text\npython reports/matched-stretch/history-prediction/spatial-convergence/analyze.py\npython reports/matched-stretch/history-prediction/spatial-convergence/build_report.py\npython -m pytest reports/matched-stretch/history-prediction/test_history_prediction.py reports/matched-stretch/history-prediction/resolution-followup/test_analysis.py reports/matched-stretch/history-prediction/mirror-followup/test_chirality.py reports/matched-stretch/history-prediction/velocity-reversal/test_reversal.py reports/matched-stretch/history-prediction/mirror-direction/test_mirror_direction.py reports/matched-stretch/history-prediction/spatial-convergence/test_refinement.py -q -p no:cacheprovider\n```\n',
             'Install the parent requirements and retain the sibling reversal/mirror archives and full reference checkpoints. The sequential runner supports --pilot and --extension; extension execution requires the saved frozen-screen decision. Preserve originals and use an isolated empty output directory for regeneration on a platform with different source bytes. Simulation cache reuse is byte-strict; analysis allows only LF/CRLF-equivalent source verification, not actual source edits.\n',
             '## Scope of any convergence conclusion\n',
             'The common t=0.3 field is a saved band-limited input. Refinement tests evolution and marker interpolation from that identical present state. It does not establish continuum accuracy of the earlier state preparation, the full history-prediction experiment, other Reynolds numbers, or the separate original strain/tube matrix. That original raw strain remains nonsmooth across periodic joins, with initial maxima outside the central tube.\n',
             'Reported events associate initial locations with deformation accumulated over a future interval. They do not establish exact Eulerian peak onset or location. These results cannot show physical material memory, a marker arrangement forcing the flow, leadership, or a Navier–Stokes breakthrough.\n',
             '![Spatial targets, event agreement and arrangement/direction effects](convergence.png)\n']
    (HERE/'README.md').write_text('\n'.join(text)+'\n', encoding='utf-8')


if __name__ == '__main__': main()
