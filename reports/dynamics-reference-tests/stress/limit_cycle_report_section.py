"""Describe the completed exact limit-cycle benchmark without rerunning it."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET


def render():
    root = Path(__file__).resolve().parent
    data = root/"data"/"limit-cycle"
    if not (data/"summary.json").exists():
        return ""
    result = json.loads((data/"summary.json").read_text())
    tests = list(ET.parse(data/"independent-tests.xml").getroot().iter("testcase"))
    failures = sum(t.find("failure") is not None or t.find("error") is not None for t in tests)
    rows = [r for r in result["cases"] if r["duration"] == 100 and r["step"] == .05]
    cells = "".join(f'<tr><td>{r["method"]}</td><td>{100*r["maximum_sampled_radial_error"]:.6g}%</td><td>{r["maximum_sampled_phase_error_degrees"]:.6g}°</td><td>{"PASS" if r["radial_accuracy_passed"] else "FAIL"}</td><td>{"PASS" if r["phase_accuracy_passed"] else "FAIL"}</td></tr>' for r in rows)
    return f'''<section id="limit-cycle"><h2>New screenshot: a rotating, attracting circle</h2>
<p>The supplied equations are <strong>θ′ = 1</strong> and <strong>r′ = (1 − r²)r</strong>. Inside the unit circle, positive radius increases; outside it, radius decreases. At radius 1, the trajectory keeps rotating, completing one revolution every 2π model units. This is an attracting periodic orbit, also called a stable limit cycle. <a href="https://uclnatsci.github.io/Mathematics-for-Natural-Sciences-2/DynamicalSystems/limit_cycles.html">UCL's course notes give the same example</a>.</p>
<p><strong>Exception to the screenshot's wording:</strong> r₀=0 stays at the origin. The angle there is undefined as a Cartesian position. Thus “every radius approaches 1” must be restricted to r₀&gt;0. The circle r=1 is invariant; points on it move rather than becoming equilibrium points.</p>
<pre>θ(t) = θ₀ + t
r(t) = r₀ / sqrt(r₀² + (1 − r₀²) exp(−2t)),   t ≥ 0
x(t) = r(t) cos(θ₀+t),   y(t) = r(t) sin(θ₀+t)</pre>
<p>The implemented formula handles r₀=0 separately and otherwise works both inside and outside the circle. In the screenshot's separated integral, the real logarithm requires |1−r²| for r&gt;1. The displayed intermediate expression with a positive exponential covers the inner branch; the final initial-value formula extends to all positive initial radii for forward time.</p>
<p>Completed <strong>{len(result['cases'])} numerical comparisons</strong>, <strong>{result['checks_passed']}/{result['checks_total']} benchmark checks</strong>, and <strong>{len(tests)-failures}/{len(tests)} independent tests</strong>. Six initial radii, including the origin and points outside the circle, were advanced for 20 units with Euler, Heun and RK4 at four steps. Twelve additional unit-circle runs lasted 100 units to expose accumulated phase drift. The computational sweep, saved outputs and one selected repeat took {result['elapsed_seconds']:.2f} seconds; plotting and the independent suite are separate.</p>
<figure><img src="figures/09_limit_cycle.png" alt="Exact spirals approach the unit circle; coarse numerical methods alter the radius and accumulate phase drift"><figcaption>Left: exact trajectories; dots mark their starts. Middle and right: numerical unit-circle starts at h=0.1. Radius and timing must be checked separately. <a href="figures/09_limit_cycle.svg">SVG</a></figcaption></figure>
<p><strong>A close-looking circle can still have wrong timing.</strong> The table uses h=0.05 over 100 units, starting at radius 1. The radial-error budget is 0.1% of the exact unit radius and the phase budget is 1°.</p>
<div class="scroll"><table><thead><tr><th>Method</th><th>Maximum radius error</th><th>Maximum accumulated phase error</th><th>Radius target</th><th>Phase target</th></tr></thead><tbody>{cells}</tbody></table></div>
<p>Heun passes the radius target but fails the phase target. All 84 cases remained within the sampled radius guard of 3, including inaccurate settings. Boundedness is therefore not sufficient validation. RK4 recovered fourth-order convergence against the exact solution. One finest-step 100-unit RK4 trajectory repeated every state exactly in this environment; all 84 cases were not independently repeated.</p>
<p class="small">This is a two-coordinate state-space model. The Cartesian equations are x′=(1−x²−y²)x−y and y′=x+(1−x²−y²)y. Their divergence is 2−4r², so they cannot be treated directly as an incompressible planar water velocity field. Radius is not conserved away from r=1; no artificial mass-conservation test is imposed. The study supplies a solver benchmark, not physical or biological measurements. Errors and boundedness are measured at saved steps, not proven between steps. Initial conditions are deterministic; no stochastic confidence intervals apply.</p>
<p><a href="data/limit-cycle/cases.csv">All 84 comparisons</a> · <a href="data/limit-cycle/summary.json">Results and source hashes</a> · <a href="data/limit-cycle/independent-tests.xml">Independent test record</a> · <a href="limit_cycle_plan.json">Settings and criteria</a> · <a href="limit_cycle.py">Reproducible model</a></p></section>'''
