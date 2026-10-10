"""Render completed adaptive transport results into the existing stress report."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET


def render():
    root = Path(__file__).resolve().parent
    path = root/"data"/"adaptive"/"summary.json"
    if not path.exists():
        return ""
    result = json.loads(path.read_text())
    rows = result["cases"]
    tests = list(ET.parse(path.parent/"independent-tests.xml").getroot().iter("testcase"))
    failures = sum(x.find("failure") is not None or x.find("error") is not None for x in tests)
    manual = [r for r in rows if r["method"] == "RK4"]
    cells = "".join(f'<tr><td>{r["speed_profile"]}</td><td>{r["schedule"]}</td><td>{r["maximum_sampled_profile_error"]:.6g}</td><td>{r["maximum_sampled_accumulated_phase_error_degrees"]:.5g}°</td></tr>' for r in manual)
    fine = [r for r in rows if r["grid"] == 64 and r["method"] == "DOP853" and r["rtol"] == 1e-9]
    return f'''<section id="adaptive"><h2>New test: bigger and smaller intervals, in both orders</h2>
<p><strong>Calculation intervals can change as the motion changes.</strong> This benchmark tests slow → fast → slow and fast → slow → fast. A smaller time step means calculating more often; it does not slow the modeled water. The prescribed flow determines its speed.</p>
<p>Completed: <strong>{len(rows)} comparisons</strong>, <strong>{result['checks_passed']}/{result['checks_total']} benchmark checks</strong> and <strong>{len(tests)-failures}/{len(tests)} additional independent tests</strong>. The comparisons took {result['elapsed_seconds']:.2f} seconds on this machine, excluding plotting and the separate test suite. These are deterministic model results, not physical measurements.</p>
<div class="scroll"><table><thead><tr><th>Flow speed order</th><th>Time-step order</th><th>Maximum sampled profile error</th><th>Accumulated phase error</th></tr></thead><tbody>{cells}</tbody></table></div>
<p>Both manual schedules use exactly <strong>4,500 RK4 steps and 18,000 field evaluations</strong> over eight model units. The small interval is 0.001 and the large interval is 0.008; switches occur at t=2 and t=6. Putting small intervals in the fast portion reduced profile error by about <strong>100-fold</strong>. Reversing the flow reversed which schedule worked better. Neither manual schedule met the strict 10⁻⁵ profile-error target; being better does not make it accurate enough.</p>
<figure><img src="figures/07_adaptive_intervals.png" alt="Two prescribed flow-speed curves and the corresponding automatically selected time intervals"><figcaption>The error-controlled solver increases and decreases its intervals in both directions. Initial-step selection and the final shortened step are visible. Speed alone does not generally determine numerical difficulty. <a href="figures/07_adaptive_intervals.svg">SVG</a></figcaption></figure>
<p><strong>Automatic selection:</strong> on grid 64, DOP853 at relative tolerance 10⁻⁹ used {fine[0]['accepted_steps']} accepted intervals in each flow profile, with worst sampled profile error {max(r['maximum_sampled_profile_error'] for r in fine):.3g}. The median interval was approximately 0.033–0.034 where speed was below 1, and 0.00665–0.00666 where speed exceeded 3. Grid 128 gave comparable accuracy. Different methods have different costs; these counts do not establish a universal speed advantage.</p>
<p><strong>Aliasing remained:</strong> grid 32 retained error around 1.99, despite tighter tolerances and a passing mass budget. On resolved grids, tightening tolerance from 10⁻⁶ to 10⁻⁹ reduced error from about 2×10⁻⁵ to 1.8×10⁻⁸. The loose tolerance missed the chosen global accuracy target. All 16 cases passed the separate mass and finite sampled boundedness checks; only the four tight-tolerance, resolved-grid cases passed every accuracy budget.</p>
<p>The periodic model is q<sub>t</sub> + ω(t)q<sub>θ</sub> = 0, q(θ,0)=cos(17θ), with physical scalar m=60+q. The exact reference is cos(17[θ−∫ω dt]). Spatially uniform angular speed makes this prescribed ring flow incompressible. There is no diffusion, particle interaction, biological force, turbulence model or dynamically solved Navier–Stokes velocity.</p>
<p>The solver evolves q so a large mean cannot hide error under relative tolerance. <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html">SciPy documents DOP853 and its local error controls</a>; solver completion alone is not a global accuracy guarantee. This study separately checks the analytic solution, mass and boundedness. Independent tests repeat both tight grid-64 adaptive integrations exactly in this environment. The full set of 16 cases was not independently repeated.</p>
<p><strong>There is no universally best fixed number pattern.</strong> Define the numbers, their units, an error budget and a cost budget before comparing them. This run tested two manual schedules and error-controlled intervals; it did not optimize arbitrary sequences. Adaptive spatial refinement is a different extension and was not implemented. Adding grid points after discarding a mode cannot uniquely reconstruct the lost information.</p>
<p class="small">Errors are sampled maxima: adaptive dense output every 0.005 model units; manual runs at every accepted step. They are not rigorous continuous bounds. The input is deterministic, so random seeds and statistical confidence intervals do not apply. Accuracy budgets: profile and amplitude errors ≤10⁻⁵; accumulated phase error ≤0.01°. Mass budget: 10⁻¹⁰. A sampled fluctuation bound of 2 is a guard, not proof of stability for all time.</p>
<p><a href="data/adaptive/cases.csv">All 16 comparisons</a> · <a href="data/adaptive/summary.json">Results and source hashes</a> · <a href="data/adaptive/independent-tests.xml">Independent test record</a> · <a href="adaptive_transport_plan.json">Settings and budgets</a> · <a href="adaptive_transport.py">Reproducible model</a></p></section>'''
