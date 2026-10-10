"""Plot and describe completed stress outputs without rerunning any model."""
from pathlib import Path
import csv
import html
import json
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import NullLocator, ScalarFormatter
from adaptive_report_section import render as adaptive_section
from sequence_report_section import render as sequence_section
from limit_cycle_report_section import render as limit_cycle_section
from pond_report_section import render as pond_section

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"
FIG=ROOT/"figures"
FIG.mkdir(exist_ok=True)


def load(name): return json.loads((DATA/name).read_text())


def table(name):
    with (DATA/name).open(newline="") as file: return list(csv.DictReader(file))


def save_figure(name,fig):
    fig.tight_layout()
    for extension in ("png","svg"): fig.savefig(FIG/f"{name}.{extension}",dpi=160,bbox_inches="tight")
    plt.close(fig)


def figures(result):
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    groups=result["long_run"]["groups"]
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for i,h in enumerate((.005,.0025,.00125)):
        group=[r for r in groups if r["step"]==h]
        duration=np.array([r["duration"] for r in group])
        rates=np.array([r["mean_largest_rate"] for r in group])
        half=np.array([r["approximate_95_t_interval"][1]-r["mean_largest_rate"] for r in group])
        axes[0].errorbar(duration*(1+(i-1)*.025),rates,yerr=half,marker="o",capsize=3,label=f"h={h}")
        axes[1].plot(duration,[r["mean_z"] for r in group],"o-",label=f"h={h}")
    axes[0].set(xlabel="Measurement duration (model units)",ylabel="Largest growth rate",title="Longer windows; eight independent starts")
    axes[1].set(xlabel="Measurement duration (model units)",ylabel="Mean z",title="Check a second observable too")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks([100,250,500,1000])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(NullLocator())
        ax.legend(fontsize=8)
    save_figure("01_duration_and_step",fig)
    rows=table("transport_convergence.csv")
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for n in (16,32,64,128):
        own=[r for r in rows if int(r["grid"])==n]
        steps=[float(r["step"]) for r in own]
        axes[0].loglog(steps,[float(r["maximum_endpoint_error"]) for r in own],"o-",label=f"grid {n}")
        axes[1].loglog(steps,[float(r["maximum_mass_error"]) for r in own],"o-",label=f"grid {n}")
    for ax in axes:
        ax.set_xticks([.00025,.0005,.001,.002])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(axis="x",labelrotation=25)
        ax.legend(fontsize=8)
    axes[0].set(xlabel="Time step",ylabel="Maximum endpoint profile error",title="Coarse grids retain a large spatial error")
    axes[1].set(xlabel="Time step",ylabel="Maximum mass-budget error",title="Mass conservation alone misses that error")
    save_figure("02_resolution_vs_conservation",fig)
    rows=table("parameter_scan.csv")
    rhos=result["parameter_scan"]["settings"]["rhos"]
    fig,axes=plt.subplots(1,2,figsize=(11,5.8))
    for ax,h in zip(axes,(.005,.0025)):
        matrix=np.zeros((len(rhos),6))
        for i,rho in enumerate(rhos):
            own=[r for r in rows if float(r["rho"])==rho and float(r["step"])==h]
            for row in own: matrix[i,int(row["start_index"])]=row["near_equilibrium_in_sampled_terminal_window"]=="True"
        ax.imshow(matrix,aspect="auto",cmap=ListedColormap(["#eadac6","#167768"]),vmin=0,vmax=1)
        ax.set(xticks=np.arange(6),xticklabels=["Near +","Near -","Start 1","Start 2","Start 3","Start 4"],
               yticks=np.arange(len(rhos)),yticklabels=[f"{r:.5g}" for r in rhos],ylabel="rho",title=f"Finite-window classification, h={h}")
        ax.tick_params(axis="x",labelrotation=35)
        for i in range(len(rhos)):
            for j in range(6): ax.text(j,i,"near" if matrix[i,j] else "open",ha="center",va="center",fontsize=8,color="white" if matrix[i,j] else "#423b32")
    fig.suptitle("Near: within 0.001 of one equilibrium throughout sampled t=180…200\nOpen: criterion not met; this does not identify chaos",fontsize=11)
    save_figure("03_parameter_scan",fig)
    rows=result["parameter_scan"]["local_times"]
    fig,ax=plt.subplots(figsize=(10,4.5))
    efold=np.array([r["linear_efold_time"] for r in rows])
    ax.bar(np.arange(len(rows)),efold,color=["#b76b20" if value>1000 else "#237d8d" for value in efold])
    ax.set(yscale="log",xticks=np.arange(len(rows)),xticklabels=[f"{r['rho']:.5g}" for r in rows],
           xlabel="rho",ylabel="Local linear e-folding time",title="Near the stability boundary, 1,000 units is a short run")
    ax.tick_params(axis="x",labelrotation=30)
    ax.axhline(200,color=".4",ls=":",label="parameter-scan duration: 200")
    ax.axhline(1000,color="#ad3e32",ls="--",label="long-run duration at rho=28: 1,000")
    ax.legend(fontsize=8)
    save_figure("04_critical_timescale",fig)
    rows=table("rk4_frequency_stress.csv")
    fig,ax=plt.subplots(figsize=(8,4))
    positions=np.arange(len(rows))
    values=[float(r["observed_final_amplitude"]) for r in rows]
    ax.scatter(positions,values,color=["#227c91","#b56b1b","#ae3230"],s=80,zorder=3)
    ax.axhline(1,color=".25",ls="--",label="exact amplitude = 1")
    ax.set(yscale="log",xticks=positions,xticklabels=[r["step"] for r in rows],xlim=(-.5,2.5),
           xlabel="RK4 time step",ylabel="Amplitude after one model time unit",
           title="Stable can be inaccurate: a pure oscillation is damped away")
    for x,value in zip(positions,values): ax.annotate(f"{value:.3g}",(x,value),xytext=(0,10),textcoords="offset points",ha="center")
    ax.set_ylim(1e-9,1e26)
    ax.legend(fontsize=8)
    save_figure("05_stability_is_not_accuracy",fig)


def main():
    result=load("summary.json")
    checks=load("checks.json")
    audit=load("selected-reproduction.json")
    diagnostic=load("reference-derivative-diagnostic.json")
    oscillation=load("oscillation-summary.json")
    tests=list(ET.parse(DATA/"independent-tests.xml").getroot().iter("testcase"))
    failures=sum(x.find("failure") is not None or x.find("error") is not None for x in tests)
    figures(result)
    fine=next(r for r in result["long_run"]["groups"] if r["step"]==.00125 and r["duration"]==1000)
    comparisons=[r for r in checks if r["kind"]=="predeclared robustness criterion"]
    cells="".join(f'<tr><td>{html.escape(r["name"])}</td><td>{r["difference"]:.5g}</td><td>{r["tolerance"]:.5g}</td><td>{"PASS" if r["passed"] else "FAIL"}</td></tr>' for r in comparisons)
    charts=[
        ("01_duration_and_step","Duration and resolution","Error bars are approximate 95% t intervals across eight independent starts. The same starts and overlapping duration prefixes are paired; blocks and prefixes are not independent experiments."),
        ("02_resolution_vs_conservation","A conservation check can miss a wrong solution","The source contains angular mode 17. Grids 16 and 32 alias that structure. Their mass budgets remain close to roundoff, yet profile errors persist as the time step shrinks. Grids 64 and 128 resolve it and recover fourth-order time convergence."),
        ("03_parameter_scan","120 finite-window trajectories","The near-equilibrium criterion uses sampled distances during the final 20 model time units, with tolerance 0.001. ‘Open’ means the criterion was not met by t=200. It does not distinguish long transients, periodic motion or chaos."),
        ("04_critical_timescale","How long is long enough?","At rho_H ± 0.01, a small linear perturbation has an e-folding time near 3,300 units. Five e-folds take about 16,500 units. These are local predictions from eigenvalues, not nonlinear simulations of that duration. At rho=28, the ensemble was tested for 1,000 units after a 50-unit transient."),
        ("05_stability_is_not_accuracy","Accuracy requires more than avoiding blow-up","For y′=−376 i y, the exact amplitude stays one. RK4 h=0.004 erases almost all amplitude while remaining stable. At h=0.008 it amplifies the mode to roughly 10²³. These deliberately inadequate steps are negative controls, not physical instability.")]
    images="".join(f'<section><h2>{title}</h2><figure><img src="figures/{name}.png" alt="{title}"><figcaption>{caption} <a href="figures/{name}.svg">SVG</a></figcaption></figure></section>' for name,title,caption in charts)
    method_labels={"rk4":"RK4","implicit_midpoint":"Implicit midpoint","gauss_legendre4":"Gauss–Legendre order 4"}
    chosen=[r for r in oscillation["choices"] if r["duration"]==100 and r["method"] in method_labels]
    step_rows="".join(f'<tr><td>{method_labels[r["method"]]}</td><td>{r["largest_tested_accepted_step"]:.9g}</td><td>{100*r["amplitude_error"]:.5g}%</td><td>{r["phase_error_degrees"]:.5g}°</td></tr>' for r in chosen)
    accuracy_section=f'''<section><h2>Changed: accuracy now has its own acceptance rules</h2>
<p>For the oscillation benchmark, a setting must satisfy <strong>0.1% relative amplitude error</strong>, <strong>1° accumulated phase error</strong>, numerical stability and at least eight output samples per true period. These are explicit pilot limits, not universal physical requirements. Phase is accumulated without wrapping away missed cycles.</p>
<p>The expanded analysis evaluated <strong>{oscillation['combinations']} method/step/duration combinations</strong>. <strong>{oscillation['accepted_combinations']}</strong> met all requirements; <strong>{oscillation['combinations']-oscillation['accepted_combinations']}</strong> were rejected. Four representative 40,000-step direct integrations independently agreed with the amplification-factor predictions. The 144-case sweep uses exact discrete-map analysis, not 144 separately time-stepped trajectories.</p>
<p><strong>Why amplitude alone is insufficient:</strong> implicit midpoint preserves amplitude for this linear mode but at h=0.00025 accumulates about <strong>158.42°</strong> of phase error over 10 units. Fourth-order Gauss–Legendre reduces that error to about <strong>0.02335°</strong> at the same step and duration. Neither result justifies automatically replacing the solver in a different model.</p>
<div class="scroll"><table><thead><tr><th>Method</th><th>Largest tested accepted step at T=100</th><th>Amplitude error</th><th>Accumulated phase error</th></tr></thead><tbody>{step_rows}</tbody></table></div>
<p class="small">Zero amplitude error is the exact-arithmetic norm identity for these implicit methods on this constant-frequency linear mode; roundoff was measured separately. Implicit stages have a different cost from explicit RK4, so step counts alone do not establish speed. The exact-rotation oracle applies only to the analytically solvable mode.</p>
<figure><img src="figures/06_amplitude_and_phase.png" alt="Amplitude and accumulated phase error across time steps after 100 model time units"><figcaption>Both error budgets must pass. Exact zero norm error is displayed at 10⁻¹⁴ for visibility. Coarse settings with insufficient phase sampling are rejected and omitted from these curves. <a href="figures/06_amplitude_and_phase.svg">SVG</a></figcaption></figure>
<p><a href="data/oscillation_accuracy.csv">All 144 accept/reject decisions</a> · <a href="data/oscillation_accepted_steps.csv">Accepted steps by duration</a> · <a href="data/oscillation_direct_audit.csv">Direct integration audits</a> · <a href="oscillation_accuracy_plan.json">Accuracy plan</a></p></section>'''
    failed="" if not result["failed_checks"] else '<p class="alert">Failed criteria: '+html.escape(", ".join(r["name"] for r in result["failed_checks"]))+'</p>'
    run_minutes=result["elapsed_seconds"]/60
    body=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Expanded dynamics stress tests</title>
<style>*{{box-sizing:border-box}}body{{margin:0;color:#193143;background:#f3f6f8;font:16px/1.6 system-ui,Segoe UI,sans-serif}}main{{max-width:1140px;margin:auto;padding:32px 24px 65px}}h1{{font-size:clamp(30px,4vw,46px);line-height:1.12;margin:10px 0 20px}}h2{{line-height:1.25;font-size:25px}}section{{background:white;border:1px solid #d5e0e8;border-radius:8px;padding:25px;margin:27px 0}}a{{color:#096781}}.eyebrow{{text-transform:uppercase;letter-spacing:.1em;color:#237487;font-size:12px;font-weight:bold}}.lead{{font-size:20px;max-width:930px}}.badges{{display:flex;flex-wrap:wrap;gap:10px}}.badges span{{background:#e4f0ed;padding:6px 12px;border:1px solid #accfc5;border-radius:5px;font-weight:650}}figure{{margin:0}}img{{width:100%;height:auto}}figcaption,.small{{font-size:14px;color:#496477}}table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{padding:10px;border-bottom:1px solid #dce5eb;text-align:left}}th{{background:#edf3f6}}.scroll{{overflow-x:auto}}pre{{background:#edf3f6;padding:17px;overflow-x:auto;font-size:13px;line-height:1.45}}.alert{{color:#a92e21;font-weight:bold}}li{{margin:8px 0}}@media print{{body{{background:white}}section{{break-inside:avoid;border:0}}main{{padding:0}}}}@media(max-width:650px){{main{{padding:20px 12px}}section{{padding:16px}}}}</style></head><body><main>
<header><div class="eyebrow">Expanded local study · 9 October 2026</div><h1>Stress duration is a convergence question</h1>
<p class="lead">The expanded run took <strong>{run_minutes:.2f} minutes of computation</strong>. It tested 24 long Lorenz trajectories, 120 parameter-scan trajectories and 16 transport grid/time combinations. The computations expose failure modes that a short, smooth benchmark can miss.</p>
<div class="badges"><span>{result['checks_passed']}/{result['checks_total']} declared checks</span><span>{len(tests)-failures}/{len(tests)} independent tests</span><span>1 finest-step long run reproduced</span></div>{failed}
<p class="small">Model-time units are dimensionless and are not seconds of a real apparatus. <a href="../report.html">Original baseline report</a> · <a href="README.md">Methods, plan and limits</a></p></header>
<section><h2>New: can a physical change turn decay into growth?</h2><p>The completed <a href="convection-transition/report.html">heated-fluid transition study</a> compares a thermally driven Boussinesq layer below and above its classical onset threshold. It includes full inertia, separate spatial/time refinement, energy budgets, negative controls, raw records and an independent inertia-removal comparison. <a href="convection-transition/README.md">Read the methods, results and limitations</a>.</p></section>
{pond_section()}
{limit_cycle_section()}
{sequence_section()}
{adaptive_section()}
<section><h2>What changed with the longer test?</h2><ul>
<li>Eight independent starting states, each at steps 0.005, 0.0025 and 0.00125; 50 units of state transient followed by 1,000 units of measurement.</li>
<li>At the finest step, the largest growth rate averaged <strong>{fine['mean_largest_rate']:.5f}</strong> per model unit, with an approximate 95% interval <strong>[{fine['approximate_95_t_interval'][0]:.5f}, {fine['approximate_95_t_interval'][1]:.5f}]</strong>.</li>
<li>Spatial underresolution produced an endpoint error of about <strong>3.80</strong> on grid 16 and <strong>0.434</strong> on grid 32, despite mass-budget errors around 10⁻¹³. At the finest step, resolved grids reduced the profile error to about <strong>1.11×10⁻¹⁰</strong>.</li>
<li>Near the stability boundary, the eigenvalues predict a roughly <strong>3,300-unit</strong> local e-folding time. The finite parameter scan does not establish long-term settling there.</li></ul></section>
<section><h2>Predeclared stopping criteria</h2><p>These pilot tolerances were written before execution. They compare the two finest steps and the 500-unit and 1,000-unit duration prefixes. Passing supports these observables at rho=28; it does not settle every parameter value.</p><div class="scroll"><table><thead><tr><th>Comparison</th><th>Absolute difference</th><th>Tolerance</th><th>Result</th></tr></thead><tbody>{cells}</tbody></table></div></section>
{accuracy_section}
{images}
<section><h2>How long should the next stress run be?</h2><p><strong>For this rho=28 pilot:</strong> assess the comparisons above before extending duration. A longer run is useful when growth-rate or state statistics still drift; refining space or time is necessary when numerical error persists. Repeating an underresolved solution for longer does not repair it.</p>
<p><strong>For parameters very close to a stability boundary:</strong> use the local timescale to design a separate run. Around rho_H ± 0.01, even five local e-folds require roughly 16,500 model units. Such a run was not completed here, and five e-folds would still not prove a global basin or bifurcation claim.</p>
<p><strong>Wall-clock planning:</strong> this machine took {run_minutes:.2f} minutes for the bounded expansion and {audit['elapsed_seconds']/60:.2f} additional minutes to repeat one selected finest-step trajectory. At fixed method, ensemble size and step, cost grows roughly with simulated duration. Doubling the number of starts or halving the step also roughly doubles integration work. These are planning estimates, not guarantees.</p></section>
<section><h2>Failures and verification are retained</h2><p>The first independent derivative check had residual <strong>{diagnostic['values'][0]['max_residual']:.3g}</strong>, marginally exceeding its 10⁻⁷ tolerance. Halving the derivative interval twice reduced the residual to <strong>{diagnostic['values'][-1]['max_residual']:.3g}</strong>. The tolerance and simulated model were unchanged. The <a href="data/independent-tests-initial.xml">initial failure record</a> and <a href="data/reference-derivative-diagnostic.json">refinement diagnostic</a> remain in the package.</p>
<p>The finest-step run for seed {audit['seed']} was repeated: all sampled states, tangent-growth blocks and final spectrum matched exactly in the same environment. This is one reproduced extended trajectory, not a repeat of all 24. The earlier baseline's complete 40-file reproduction remains separately documented.</p>
<p><a href="data/checks.json">Every declared check</a> · <a href="data/long_runs.csv">Long-run results</a> · <a href="data/duration_prefixes.csv">Duration prefixes</a> · <a href="data/parameter_scan.csv">Parameter classifications</a> · <a href="data/selected-reproduction.json">Selected repeat</a> · <a href="stress_plan.json">Original plan</a></p></section>
<section><h2>Reproduce</h2><p>Use the pinned dependencies in the parent package. From that package directory:</p><pre>python stress/run_stress.py
python stress/oscillation_accuracy.py
python stress/adaptive_transport.py
python stress/sequence_transport.py
python stress/limit_cycle.py
python stress/pond_vibration.py
python -m pytest stress/test_pond_vibration.py -q --junitxml=stress/data/pond-vibration/independent-tests.xml
python -m pytest stress/test_limit_cycle.py -q --junitxml=stress/data/limit-cycle/independent-tests.xml
python stress/verify_sequence_reproduction.py
python -m pytest stress/test_sequence_transport.py -q --junitxml=stress/data/sequence/independent-tests.xml
python -m pytest stress/test_adaptive_transport.py -q --junitxml=stress/data/adaptive/independent-tests.xml
python -m pytest stress/test_stress.py stress/test_oscillation_accuracy.py -q --junitxml=stress/data/independent-tests.xml
python stress/audit_stress.py
python stress/make_stress_report.py
python make_report.py
python verify_package.py</pre>
<p class="small">If interrupted, use <code>python stress/run_stress.py --resume</code>. Only completed long-run checkpoints with identical model code and plan are reused. Cheap diagnostics run again. Preserve an earlier copy before regenerating outputs.</p></section>
<footer class="small">No new physical or biological measurements. No proof of a new law or global Navier–Stokes result. External scientific context and primary sources are documented in the <a href="../README.md">baseline README</a>.</footer></main></body></html>'''
    (ROOT/"report.html").write_text(body,encoding="utf-8")
    print(f"Stress report generated: {result['checks_passed']}/{result['checks_total']} declared checks; {len(tests)-failures}/{len(tests)} independent checks.")


if __name__=="__main__": main()
