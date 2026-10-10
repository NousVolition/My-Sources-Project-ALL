"""Build an offline, readable HTML report from saved numerical results."""
from pathlib import Path
import hashlib
import html
import json
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parent


def read_json(name):
    return json.loads((ROOT/"data"/name).read_text(encoding="utf-8"))


def main():
    summary=read_json("summary.json")
    checks=read_json("checks.json")
    environment=read_json("environment.json")
    reproduction=read_json("reproduction.json")
    suites=ET.parse(ROOT/"data"/"independent-tests.xml").getroot()
    cases=list(suites.iter("testcase"))
    failures=sum(case.find("failure") is not None or case.find("error") is not None for case in cases)
    independent_passed=len(cases)-failures
    rates=summary["lorenz_ensemble"]["estimates"][-1]
    wheel=summary["waterwheel"]["finest"]
    f=lambda number: f"{number:.5g}"
    figures=[
        ("01_numerical_methods","Check the method before interpreting the motion",
         "Euler, improved Euler and RK4 recover their expected convergence orders. The positive logistic ODE approaches its equilibrium smoothly. Euler steps that are too large create oscillation and loss of equilibrium stability."),
        ("02_scalar_stability","Four kinds of equilibrium behavior",
         "Representative fields −x³, x³, x² and 0 give attraction, repulsion, attraction from one side, and a line of equilibria. The derivative at zero vanishes in the first three, so the sign of the field supplies information that linearization misses. Four initial values were tested for each field."),
        ("03_nonuniqueness","One initial value can have many solutions here",
         "Six exact departing or waiting solutions and the stationary solution all start at zero for x′=real_cuberoot(x). A solver starting at exactly zero returned zero. That result cannot establish uniqueness. The perturbed examples use a fixed step and do not uniformly resolve the limit as the initial value approaches zero."),
        ("04_overdamped_limit","Neglecting inertia becomes a better approximation",
         "The full damped oscillator agrees with a matrix-exponential benchmark and its energy budget closes. Reducing mass at fixed drag and spring constant brings its displacement closer to exp(−t). The original velocity is zero, so an initial adjustment layer remains."),
        ("05_bifurcation_and_reaction","Equilibria depend on parameters and assumptions",
         "The saddle-node example has two equilibria above r=1, with opposite stability. The cropped reaction example is tested under two explicit assumptions: a closed supply of A+X, and A maintained by an external reservoir. These are illustrative mass-action equations, not measured bacterial kinetics."),
        ("06_lorenz_volume_and_separation","Shrinking state volume and growing separation coexist",
         "The integrated tangent matrix obeys log(V/V₀)=−(σ+1+β)t. A separate pair of initially close trajectories eventually separates. This volume is a volume of possible system states, not the physical volume of a parcel of water."),
        ("07_lorenz_ensemble","Sensitivity survives repeated starts and step refinement",
         "Four independent starting states were each run at two time steps. All eight largest finite-time Lyapunov rates were positive. The plotted orbit uses a densely sampled 30-unit segment; sparse QR samples were not joined to draw the orbit."),
        ("08_waterwheel_budget","The waterwheel closes its mass budget",
         "Numerical total mass follows the exact inflow-minus-leakage law. The first harmonic and rotation agree with Lorenz β=1; an unforced third harmonic decays without feeding the torque. This smooth validation case is not a simulation of a particular apparatus.")]
    chart_html="\n".join(f'<section><h2>{html.escape(title)}</h2><figure><img src="figures/{name}.png" alt="{html.escape(title)}" loading="lazy"><figcaption>{html.escape(caption)} <a href="figures/{name}.svg">Vector figure</a></figcaption></figure></section>' for name,title,caption in figures)
    checks_html="\n".join(f'<tr><td>{html.escape(item["name"])}</td><td class="{"pass" if item["passed"] else "fail"}">{"PASS" if item["passed"] else "FAIL"}</td></tr>' for item in checks)
    report=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Dynamics tests — numerical results</title>
<style>
:root{{--ink:#172c3c;--muted:#4c6578;--line:#d7e1e8;--accent:#176777;--paper:#f3f6f8}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.6 system-ui,Segoe UI,sans-serif}}
main{{max-width:1180px;margin:auto;padding:38px 28px 70px}}header{{padding:10px 0 30px;border-bottom:2px solid var(--accent)}}
h1{{font-size:clamp(30px,4vw,47px);line-height:1.12;letter-spacing:-.025em;margin:10px 0 16px}}h2{{font-size:24px;line-height:1.25;margin:0 0 16px}}
h3{{font-size:18px;margin:25px 0 8px}}p{{margin:10px 0 16px}}a{{color:#075f85;text-underline-offset:3px}}.eyebrow{{color:var(--accent);font-weight:700;letter-spacing:.09em;text-transform:uppercase;font-size:12px}}
.lead{{font-size:20px;max-width:880px}}.status{{display:flex;flex-wrap:wrap;gap:10px;margin:22px 0 5px}}.status span{{padding:6px 13px;border:1px solid #b7d5ca;background:#e8f4ee;border-radius:5px;font-weight:650}}
section{{margin:30px 0;background:white;padding:27px;border:1px solid var(--line);border-radius:8px}}.grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:15px;margin-top:25px}}
.card{{border-top:3px solid var(--accent);padding:15px;background:white}}.value{{font-size:28px;line-height:1.2;font-weight:700;margin:5px 0}}.small,figcaption{{font-size:14px;color:var(--muted)}}
table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:11px 12px;border-bottom:1px solid var(--line);vertical-align:top}}th{{background:#edf3f6}}
.scroll{{overflow-x:auto}}figure{{margin:0}}img{{display:block;width:100%;height:auto}}figcaption{{padding-top:12px}}li{{margin:6px 0}}code,pre{{font-family:Consolas,monospace}}pre{{background:#eef3f6;padding:18px;overflow-x:auto;line-height:1.45;font-size:13px}}
.pass{{color:#12683c;font-weight:700}}.fail{{color:#b52020;font-weight:700}}summary{{cursor:pointer;font-weight:700}}footer{{font-size:13px;color:var(--muted)}}
@media(max-width:700px){{main{{padding:22px 14px}}section{{padding:18px}}.grid{{grid-template-columns:1fr}}}}
@media print{{body{{background:white}}main{{max-width:none;padding:0}}section{{break-inside:avoid;border:0;padding:10px 0}}.grid{{grid-template-columns:repeat(3,1fr)}}details{{display:none}}}}
</style></head><body><main>
<header><div class="eyebrow">Computational study · 9 October 2026 · Local results</div>
<h1>Stability, numerical accuracy<br>and the Lorenz waterwheel</h1>
<p class="lead">The supplied examples pass analytic and independent numerical checks. They also show why solver artifacts, nonuniqueness and chaotic sensitivity must be distinguished.</p>
<div class="status"><span>{summary['checks']['passed']}/{summary['checks']['total']} model checks passed</span><span>{independent_passed}/{len(cases)} independent tests passed</span><span>{reproduction['numerical_files_compared']} numerical files reproduced</span></div>
<p class="small">All results below are mathematical computations. No laboratory or biological tests were performed. Full assumptions and equations: <a href="README.md">research README</a>.</p></header>
<div class="grid"><div class="card"><div class="eyebrow">Method accuracy</div><div class="value">1.00 · 1.98 · 3.98</div><div>Observed logistic convergence orders for Euler, Heun and RK4.</div></div>
<div class="card"><div class="eyebrow">Lorenz sensitivity</div><div class="value">{rates['mean_largest_rate']:.3f}</div><div>Mean largest growth rate per model time unit at the finer step.</div></div>
<div class="card"><div class="eyebrow">Waterwheel conservation</div><div class="value">{wheel['mass_error']:.2e}</div><div>Maximum absolute mass-budget error in the finest run.</div></div></div>
<section><h2>What the tests establish</h2><div class="scroll"><table><thead><tr><th>Test</th><th>Measured result</th><th>Meaning</th></tr></thead><tbody>
<tr><td>Logistic, h=0.025</td><td>Euler error 3.20e−3; Heun 2.74e−5; RK4 4.81e−10</td><td>Expected error orders recovered against an exact solution over t=0…4.</td></tr>
<tr><td>Coarse Euler, h=2.8</td><td>25 equilibrium crossings in the final 30 samples</td><td>The method introduced oscillation absent from the positive continuous logistic solution.</td></tr>
<tr><td>Overdamped limit</td><td>Error falls 0.07093 → 0.009313 → 0.0009831</td><td>Reducing mass 0.1 → 0.01 → 0.001 improves the approximation at fixed drag and spring constant.</td></tr>
<tr><td>Cube-root equation</td><td>7 branches; maximum residual {summary['nonuniqueness']['max_residual']:.2e}</td><td>Different exact solutions have the same initial value. A zero solver result does not imply uniqueness.</td></tr>
<tr><td>Lorenz Hopf threshold</td><td>{summary['lorenz_local']['hopf']:.10f}</td><td>A local eigenvalue crossing was verified; global bifurcations were not computed.</td></tr>
<tr><td>Lorenz volume contraction</td><td>Rate −13.66667; log-volume error {summary['lorenz_local']['phase_volume_log_error']:.2e}</td><td>Contraction of state-space volume coexists with a positive largest growth rate.</td></tr>
<tr><td>Waterwheel → Lorenz</td><td>β=1; projection error {wheel['lorenz_projection_error']:.2e}</td><td>The idealized first-mode reduction differs from the classic β=8/3 Lorenz example.</td></tr>
</tbody></table></div></section>
{chart_html}
<section><h2>Uncertainty, controls and limits</h2>
<p>Four independently seeded initial states were each tested with Lorenz steps 0.005 and 0.0025. After 30 time units of transient, the growth spectrum was measured for 180 time units. At the finer step the approximate 95% interval across starts is <strong>[{rates['approximate_95_t_interval'][0]:.5f}, {rates['approximate_95_t_interval'][1]:.5f}]</strong>. It is a finite-run summary, not physical measurement uncertainty or a bound on every source of bias. The two time steps use paired starts.</p>
<p>Exact benchmarks, a separate high-order solver, matrix exponentials and conservation identities test numerical accuracy. Waterwheel ablations remove torque, inflow asymmetry or the extra harmonic. Three grids and three steps are compared. Only low Fourier modes are populated, so the grid comparison demonstrates that this smooth case is resolved.</p>
<p>Unstable periodic orbits, subcritical Hopf criticality, the homoclinic threshold, global attractor coexistence and basin boundaries were not independently computed. These would require additional continuation and basin tests. Fixed-point exercise answers that follow from sign arguments are documented separately in the README and are not counted as numerical tests.</p>
<p><strong>For the fluid research:</strong> phase-space contraction is not physical compression of water. The wheel is an open system with inflow and leakage. These results do not establish a new Navier–Stokes law, prove uniqueness for the full fluid equations, or measure mixing, biological stress response or freezing. The original <a href="../water-biology-study/report.html">three-question water/biology study</a> remains a separate package.</p></section>
<section><h2>Reproduce and inspect</h2><p>The package includes raw CSV and NPZ outputs, all equations and solvers, SVG/PNG figures, dependency versions and test records. A full repeated run returned identical numerical contents for <strong>{reproduction['numerical_files_compared']} files</strong> on this environment. Exact chaotic trajectories need not be identical on different hardware or library builds.</p>
<pre>python -m pip install -r requirements-lock.txt
python run_tests.py
python -m pytest test_independent.py -q --junitxml=data/independent-tests.xml
python verify_reproduction.py
python make_report.py
python verify_package.py</pre>
<p class="small">Use a virtual environment; Windows and macOS/Linux instructions are in the <a href="README.md">README</a>. Recorded versions: Python {environment['python'].split()[0]}, NumPy {environment['numpy']}, SciPy {environment['scipy']}, Matplotlib {environment['matplotlib']}.</p>
<p><a href="data/summary.json">Numerical summary</a> · <a href="data/checks.json">All check measurements</a> · <a href="data/reproduction.json">Reproduction record</a> · <a href="data/lorenz_ensemble.csv">Seed-level growth rates</a> · <a href="data/waterwheel_convergence.csv">Grid/time comparisons</a></p>
<details><summary>View all {len(checks)} model checks</summary><table><thead><tr><th>Check</th><th>Result</th></tr></thead><tbody>{checks_html}</tbody></table></details></section>
<section><h2>Source material and evidence</h2><p>The supplied textbook screenshots specify the examples. Partial crops are explicitly marked where assumptions were needed. The original page images are not republished.</p>
<p><a href="https://doi.org/10.1175/1520-0469(1963)020%3C0130:DNF%3E2.0.CO;2">Lorenz (1963), Deterministic Nonperiodic Flow</a> supplies the primary context for the classic model. <a href="https://arxiv.org/abs/1202.5508">Illing and colleagues (2012), Experiments with a Malkus-Lorenz water wheel</a> is external experimental evidence concerning a physical apparatus. No measurements from those papers are presented here as new data or reproduced experimental observations.</p></section>
<footer>Generated entirely from this local study package. No external scripts, fonts or analytics. Mathematical model results; all simulation parameters and limitations are documented.</footer>
</main></body></html>'''
    stress_summary=ROOT/"stress"/"data"/"summary.json"
    if stress_summary.exists() and (ROOT/"stress"/"report.html").exists():
        expanded=json.loads(stress_summary.read_text())
        extension=f'''<section><h2>Expanded stress results are available</h2>
<p>The next pass took <strong>{expanded['elapsed_seconds']/60:.2f} minutes</strong>: 24 long Lorenz runs, 120 parameter-scan trajectories, tougher spatial resolution checks and separate amplitude/phase accuracy rules.</p>
<p><strong>{expanded['checks_passed']}/{expanded['checks_total']} declared stress checks passed.</strong> The oscillation study additionally evaluates 144 method/step/duration combinations and explicitly rejects stable-but-inaccurate settings.</p>
<p><a href="stress/report.html">Open the expanded illustrated report</a> · <a href="stress/README.md">Stress methods and reproduction</a></p></section>'''
        report=report.replace('</header>','</header>'+extension,1)
    (ROOT/"report.html").write_text(report,encoding="utf-8")
    excluded={"manifest.json","package-verification.json"}
    paths=[p for p in ROOT.rglob("*") if p.is_file() and p.name not in excluded and not any(part in ("__pycache__",".pytest_cache",".venv") for part in p.parts)]
    manifest={p.relative_to(ROOT).as_posix():dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(paths)}
    (ROOT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(f"Report generated; {len(manifest)} files recorded in manifest.")


if __name__=="__main__": main()
