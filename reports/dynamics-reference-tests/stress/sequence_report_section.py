"""Render existing sequence results without rerunning numerical computations."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET


def render():
    root = Path(__file__).resolve().parent
    data = root/"data"/"sequence"
    if not (data/"summary.json").exists():
        return ""
    result = json.loads((data/"summary.json").read_text())
    rows = result["cases"]
    finest = max(r["cycles"] for r in rows)
    shown = [r for r in rows if r["grid"] == 64 and r["cycles"] == finest]
    cells = "".join(f'<tr><td>{r["speed_profile"]}</td><td>{r["schedule"]}</td><td>{r["maximum_sampled_profile_error"]:.4g}</td><td>{r["maximum_sampled_accumulated_phase_error_degrees"]:.4g}°</td><td>{"PASS" if r["accuracy_passed"] else "FAIL"}</td></tr>' for r in shown)
    tests = list(ET.parse(data/"independent-tests.xml").getroot().iter("testcase"))
    failed = sum(t.find("failure") is not None or t.find("error") is not None for t in tests)
    reproduction = "Two selected finest grid-64 pattern trajectories repeated every state exactly in this environment. A repeat of all 42 saved cases was attempted but did not start because automatic approval review hit the account usage limit. It is not a completed verification."
    if (data/"reproduction.json").exists():
        repeat = json.loads((data/"reproduction.json").read_text())
        agreement = "Every saved numerical array matched exactly" if repeat["every_retained_array_identical"] else "The numerical arrays did not all match; inspect the reproduction record"
        reproduction = f'<strong>Full repeat completed:</strong> {repeat["numerical_files"]} case files and {repeat["arrays_compared"]} arrays compared in {repeat["elapsed_seconds"]:.2f} seconds. {agreement}. This covers every saved time and diagnostic value, plus spatial states at cycle boundaries. Two selected cases also repeat every intermediate state. This is same-environment reproducibility; accuracy is assessed separately against the analytic solution. <a href="data/sequence/reproduction.json">Full reproduction record</a>.'
    return f'''<section id="sequence"><h2>Your exact 1 → 2 pattern: tested</h2>
<p>The supplied pairs <strong>12 12 12 22 11 22 11 22 32 23</strong> mean this repeated sequence:</p>
<pre>1, 2, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 3, 2, 2, 3</pre>
<p>Each number multiplies a small base time interval. These are not literal intervals of 1, 2 and 3 model units. The 20-step cycle sums to 35 base intervals. We tested the pattern, its full reverse and a uniform interval of 1.75 base units, giving identical duration and calculation counts within each comparison.</p>
<p><strong>Result:</strong> the pattern converges when refined, but in these two smooth flows uniform steps gave about <strong>2.6 times smaller error</strong> at equal work. Reversing the pattern had very little effect. A repeating pattern does not respond to the current error or flow conditions, so it is different from automatic step selection.</p>
<figure><img src="figures/08_sequence_comparison.png" alt="The exact 20-step pattern and error against calculation count for two speed profiles"><figcaption>RK4, grid 64, identical duration and work at each refinement level. The user and reversed curves nearly coincide. All values come from completed computations. <a href="figures/08_sequence_comparison.svg">SVG</a></figcaption></figure>
<p>Completed <strong>{len(rows)} comparisons</strong>, <strong>{result['checks_passed']}/{result['checks_total']} benchmark checks</strong> and <strong>{len(tests)-failed}/{len(tests)} new independent tests</strong>. The model comparisons, saved outputs and two selected repeats took {result['elapsed_seconds']:.2f} seconds; plotting and the independent test suite are additional. No new long-duration Lorenz or near-boundary run was performed.</p>
<div class="scroll"><table><thead><tr><th>Flow speed order</th><th>Schedule</th><th>Maximum sampled profile error</th><th>Maximum phase error</th><th>All accuracy budgets</th></tr></thead><tbody>{cells}</tbody></table></div>
<p>This table uses 1,280 cycles, or <strong>25,600 RK4 steps and 102,400 field evaluations per case</strong>, over eight model units. The base interval is 1/5600 ≈ 0.000178571, so pattern intervals are approximately 0.000178571, 0.000357143 and 0.000535714. The equal-work uniform interval is 0.0003125.</p>
<p>On grid 64, 80, 160, 320, 640 and 1,280 cycles were tested in both speed orders. Refinement recovered approximately fourth-order time convergence. The 640-cycle pattern already met the selected profile-error target, narrowly; 1,280 cycles provided more margin. At the finest level, grids 32 and 128 were also tested. Grid 128 agreed in accuracy with grid 64, while grid 32 retained aliasing error around 1.99 for all three schedules.</p>
<p>Accuracy, conservation and sampled boundedness remain separate. All cases passed the latter two checks, including inaccurate coarse cases. Budgets are unchanged: profile and amplitude error ≤10⁻⁵, accumulated phase error ≤0.01°, mass drift ≤10⁻¹⁰ and sampled fluctuation magnitude ≤2. These are pilot criteria, not universal standards. Every accepted time and diagnostic value is saved; spatial states are saved at cycle boundaries. Reported maxima are sampled, not rigorous continuous bounds.</p>
<p>{reproduction}</p>
<p class="small">This is the same exact-reference, mode-17 prescribed-advection model as the adaptive study below. Results are deterministic calculations, not observations of particles or bacteria. No random-seed uncertainty applies. They do not establish a universally optimal pattern, a physical speed-control mechanism or a new fluid law.</p>
<p><a href="data/sequence/cases.csv">All 42 comparisons</a> · <a href="data/sequence/summary.json">Results and provenance</a> · <a href="data/sequence/independent-tests.xml">11-test record</a> · <a href="sequence_plan.json">Exact sequence and plan</a> · <a href="sequence_transport.py">Reproducible model</a></p></section>'''
