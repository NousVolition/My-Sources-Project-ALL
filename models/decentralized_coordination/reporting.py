"""Static, exportable scientific plots and a self-contained local report page."""
import html
import os
from pathlib import Path


def make_report(out, summary, rounds, traces, exact):
    # Keep font caches within the selected output, never in the user's home.
    os.environ.setdefault("MPLCONFIGDIR", str(out / ".plot-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import numpy as np

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "#fcfcfa", "axes.facecolor": "#fcfcfa",
                         "savefig.facecolor": "#fcfcfa", "svg.hashsalt": "coordination-study"})
    colors = ["#8056a4", "#cf7832", "#187e8b"]
    labels = ["A · hierarchy", "B · local rule", "C · voluntary model"]

    def save(fig, name):
        fig.savefig(out / (name + ".png"), dpi=160, bbox_inches="tight")
        fig.savefig(out / (name + ".svg"), bbox_inches="tight", metadata={"Date": None})
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    fig.suptitle("SIMULATION ONLY · Nine agents, 120-second limit", fontsize=17, weight="bold")
    stats = summary["main"]
    probabilities = np.array([s["consensus_proportion"] for s in stats])
    errors = np.array([[p - s["consensus_ci95"][0], s["consensus_ci95"][1] - p]
                       for p, s in zip(probabilities, stats)]).T
    axes[0, 0].bar(labels, probabilities, color=colors, yerr=errors, capsize=4)
    axes[0, 0].set(ylim=(0, 1.13), ylabel="Proportion", title="Consensus by 120 s · 95% Monte Carlo intervals")
    for i, p in enumerate(probabilities):
        axes[0, 0].text(i, p + 0.035, f"{p:.1%}", ha="center", weight="bold")
    for c, label, color in zip("ABC", labels, colors):
        data = [r for r in rounds if r["condition"] == c]
        times = [r["consensus_time"] for r in data if r["consensus"]]
        t = np.arange(121)
        axes[0, 1].step(t, [sum(v <= x for v in times) / len(data) for x in t],
                        where="post", label=label, color=color, linewidth=2)
    axes[0, 1].set(xlabel="Seconds", ylabel="Fraction of all starts", ylim=(0, 1.05),
                   title="Cumulative first consensus · failures stay in denominator")
    axes[0, 1].legend(frameon=False, loc="lower right")
    axes[1, 0].bar(labels, [s["mean_capped_seconds"] for s in stats], color=colors)
    axes[1, 0].set(ylabel="Seconds", title="Mean time capped at 120 s · failures count as 120")
    x = np.arange(3)
    axes[1, 1].bar(x - 0.18, [s["mean_switches"] for s in stats], width=0.36, color="#187e8b", label="Switches")
    axes[1, 1].bar(x + 0.18, [s["mean_refusals"] for s in stats], width=0.36, color="#af5152", label="Modeled refusals")
    axes[1, 1].set(xticks=x, xticklabels=labels, ylabel="Mean per round", title="Actions before stopping")
    axes[1, 1].legend(frameon=False)
    save(fig, "condition-comparison")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    grid = np.zeros((4, 3))
    for i, refusal in enumerate((0, 0.15, 0.5, 0.9)):
        for j, innovation in enumerate((0, 0.02, 0.1)):
            name = f"r={refusal:g};i={innovation:g}"
            s = next(s for s in summary["sensitivity"] if s["label"] == name)
            grid[i, j] = s["consensus_proportion"]
    axes[0].imshow(grid, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
    for i in range(4):
        for j in range(3):
            axes[0].text(j, i, f"{grid[i,j]:.1%}", ha="center", va="center",
                         color="white" if grid[i,j] > .55 else "black")
    axes[0].set(xticks=range(3), xticklabels=["0", ".02", ".10"], yticks=range(4),
                yticklabels=["0", ".15", ".50", ".90"], xlabel="Initiative probability per opportunity",
                ylabel="Refusal probability given differing copy proposal", title="C: consensus fraction depends on assumptions")
    times = [0, 10, 20, 30, 40]
    counts = [sum(r["consensus_time"] == t for r in exact) for t in times]
    counts.append(sum(r["consensus_time"] is None for r in exact))
    axes[1].bar(["0", "10", "20", "30", "40", "Never"], counts, color=[colors[1]]*5 + ["#596575"])
    for i, value in enumerate(counts):
        axes[1].text(i, value + 5, str(value), ha="center")
    axes[1].set(xlabel="Seconds to consensus", ylabel="Number of initial states", ylim=(0, 400),
                title="B: exact enumeration of all 512 starts")
    fig.suptitle("SIMULATION / EXACT RULE ANALYSIS · No empirical human observations", fontsize=14)
    save(fig, "sensitivity-and-exact-rule")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained", sharey=True)
    for ax, (key, data) in zip(axes.flat, summary["influence"].items()):
        x = np.arange(9)
        ax.bar(x - .19, data["mean_by_seat"], width=.38, label="Position", color="#187e8b")
        ax.bar(x + .19, data["mean_by_identity"], width=.38, label="Identity", color="#bd7751")
        gains = data["heldout_mse_improvement"]
        ax.set(xticks=x, xticklabels=[str(i) for i in range(9)], xlabel="Seat number / identity number",
               ylabel="Mean switches copied from this source per round",
               title=f"{key}\nHeld-out MSE gain: seat {gains['position']:.1%}, identity {gains['identity']:.1%}")
        ax.legend(frameon=False)
    fig.suptitle("MODEL ATTRIBUTION · Seat 0 is the star hub; identity 0 has traits only in heterogeneous runs", fontsize=13)
    save(fig, "identity-versus-position")

    fig, axes = plt.subplots(4, 3, figsize=(12, 10), layout="constrained")
    for trace in traces:
        ax = axes[trace["block"], "ABC".index(trace["condition"])]
        timeline = np.full((9, 121), np.nan)
        for tick, point in enumerate(trace["trajectory"]):
            start = point["second"]
            end = trace["trajectory"][tick+1]["second"] if tick+1 < len(trace["trajectory"]) else start+1
            timeline[:, start:end] = np.array(point["state"])[:, None]
        cmap = ListedColormap(["#bb4c53", "#377dab"])
        cmap.set_bad("#e9e9e3")
        ax.imshow(timeline, origin="lower", aspect="auto", extent=(-.5,120.5,-.5,8.5), vmin=0, vmax=1, cmap=cmap)
        result = str(trace["consensus_time"]) + " s" if trace["consensus"] else "censored at 120 s"
        ax.set(title=f"Block {trace['block']+1} · {trace['condition']} · {result}",
               xticks=[0,40,80,120], yticks=[0,4,8], xlabel="Seconds", ylabel="Seat")
    fig.suptitle("FOUR SIMULATED PILOT BLOCKS · matched starts · gray means round has stopped", fontsize=14)
    save(fig, "pilot-trajectories")

    def table(headers, rows):
        return "<table><thead><tr>" + "".join(f"<th>{html.escape(str(v))}</th>" for v in headers) + "</tr></thead><tbody>" + "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row) + "</tr>" for row in rows) + "</tbody></table>"

    result_table = table(["Condition", "Consensus", "95% MC interval", "Mean capped seconds", "Mean switches", "Mean refusals"],
                         [[s["label"], f"{s['consensus_proportion']:.2%}",
                           "–".join(f"{v:.2%}" for v in s["consensus_ci95"]),
                           f"{s['mean_capped_seconds']:.2f}", f"{s['mean_switches']:.2f}", f"{s['mean_refusals']:.2f}"] for s in stats])
    influence_table = table(["Network / traits", "Position prediction gain", "Identity prediction gain", "Both", "Top identity persists", "Top seat persists"],
                            [[key, *[f"{s['heldout_mse_improvement'][v]:.2%}" for v in ("position", "identity", "both")],
                              *[f"{s['top_copy_source_persistence'][v]:.2%}" for v in ("identity", "seat")]]
                             for key, s in summary["influence"].items()])
    sensitive_table = table(["Scenario", "Consensus", "Capped seconds", "Mean refusals", "Former leader's copy share after A"],
                            [[s["label"], f"{s['consensus_proportion']:.2%}", f"{s['mean_capped_seconds']:.2f}",
                              f"{s['mean_refusals']:.2f}", f"{s['mean_former_leader_share_after_A']:.2%}"
                              if s["mean_former_leader_share_after_A"] is not None else "—"] for s in summary["sensitivity"]])
    exact_s = summary["B_exact"]
    content = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Coordination with the freedom to refuse — simulation study</title><style>
body{{margin:0;background:#f5f5ef;color:#202c38;font:17px/1.6 system-ui,sans-serif}}main{{max-width:1120px;margin:40px auto;padding:32px;background:white;border-radius:14px}}h1{{font-size:42px;line-height:1.15}}h2{{margin-top:42px}}.badge{{color:#9d3c26;font-weight:750;letter-spacing:.08em}}.callout{{padding:20px;background:#eaf5f3;border-left:5px solid #187e8b}}img{{width:100%;height:auto}}table{{border-collapse:collapse;font-size:14px;width:100%}}th,td{{padding:10px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#ecf0f2}}.scroll{{overflow-x:auto}}a{{color:#146b80}}small{{color:#5a6772}}@media(max-width:650px){{main{{margin:0;padding:20px}}h1{{font-size:32px}}}}</style>
<main><div class="badge">SIMULATION ONLY · NO HUMAN DATA</div><h1>Coordination with the freedom to refuse</h1>
<p>Nine participants in three specified conditions. Reproducible model experiments, exact rule checks, and a protocol for future volunteer testing.</p>
<div class="callout"><b>What the calculations support:</b> local coordination can reach agreement under the stated voluntary-choice model, but the result changes with refusal, initiative, and update timing. A structurally privileged seat can look like a leader; persistent individual traits can also produce influence. These simulations do not decide which explains leadership in real groups.</div>
<h2>The three conditions</h2><p><b>A:</b> a designated leader broadcasts their initial color, followers act every 10 s, with perfect compliance in the primary benchmark. <b>B:</b> a simultaneous two-neighbor rule on a ring every 10 s. <b>C:</b> each second, agents in shuffled order independently have a 25% opportunity to act; they initiate with probability 2%, otherwise select a neighbor and may decline a differing color with probability 15%.</p>
<p>B is a water-inspired analogy. It contains no molecules, forces, thermodynamics, fluid equations, or molecular dynamics. C contains modeled decisions, not a measurement of free will. A, B and C also differ in information and timing, so their difference cannot isolate freedom alone.</p>
<h2>Matched main experiment</h2><p>{stats[0]['rounds']:,} randomized starting-state blocks × three conditions; condition order balanced and shuffled; identities permuted among positions. Initial unanimity is retained and counted at zero seconds. All rounds stop at first unanimity or 120 seconds.</p>
<div class="scroll">{result_table}</div><p><small>Intervals quantify Monte Carlo uncertainty under this model, not uncertainty about people. Timeouts are censored; “capped seconds” includes them as 120 s. First consensus is a stopping event, not evidence of lasting agreement.</small></p>
<img src="condition-comparison.png" alt="Consensus probability, cumulative time, capped time, and action counts in A B C">
<h2>Rule B can freeze in disagreement</h2><p>Exact enumeration: {exact_s['successes_120']}/512 ({exact_s['successes_120']/512:.2%}) reach consensus. All 512 starts reach a fixed point within {exact_s['maximum_transient_ticks']*10} s; there are {exact_s['cycle_basins']} persistent-cycle basins. A neighboring same-color pair cannot be broken by B. Starts containing a red pair and a blue pair therefore cannot reach consensus. Temporary reversals are distinct from a lasting cycle.</p>
<img src="sensitivity-and-exact-rule.png" alt="Sensitivity heatmap for refusal and initiative alongside exact consensus times for rule B">
<h2>Position versus identity</h2><p>The original ring has equivalent seats. The separate star-network extension gives seat 0 eight neighbors and every other seat one. The heterogeneous sensitivity assigns identity 0 a donor-selection weight of 4, refusal of .65, and initiative of .03. These advantages are stipulated. Each identity occupies every seat once per nine-mapping replicate.</p>
<img src="identity-versus-position.png" alt="Mean copy-source counts by identity and seat for ring and star networks with homogeneous and heterogeneous agents">
<div class="scroll">{influence_table}</div><p>In the homogeneous star, position improves held-out prediction by {summary['influence']['star/homogeneous']['heldout_mse_improvement']['position']:.2%}; identity provides no improvement. In the homogeneous ring, neither helps. The inserted identity traits improve ring prediction by only {summary['influence']['ring/heterogeneous']['heldout_mse_improvement']['identity']:.2%}; most round-level variation remains unexplained by fixed identity or seat. Thus network position explains the dominant apparent leadership in this star model, while individual differences have a small, assumption-dependent effect here.</p><p>Prediction gain = 1 − held-out model MSE / held-out intercept MSE. Entire crossed replicates are held out; negative values mean worse prediction. “Influence” means another agent actually switched to a source's color in the model, not authority, intent, or experimentally established causation. Top-source persistence splits ties, skips rounds without copies, and compares adjacent seat permutations. Its 1/9 exchangeable reference is not a significance threshold because starts are matched. Cyclic seat permutations also preserve neighbor identity pairs on the ring; they do not test changing social relationships.</p>
<h2>Assumptions and persistence of hierarchy</h2><div class="scroll">{sensitive_table}</div><p>Former-leader scenarios increase local attention to that identity only when A occurs before C. The baseline has no status memory. A rise in persistence when memory is inserted demonstrates a mechanism, not discovery of a human tendency. All-refuse/no-initiative is a counterexample to guaranteed consensus. Initiative can escape frozen patterns but also disrupt near-consensus.</p>
<h2>Four pilot rounds per condition</h2><p>These twelve rounds are simulated feasibility examples, excluded from the main summaries. They are not a sufficient human sample. Each row uses the same colors and seating across A/B/C; execution order is recorded separately.</p><img src="pilot-trajectories.png" alt="Four matched pilot blocks with seat colors over time for all three conditions">
<h2>What remains untested</h2><p>No volunteer data were collected. There is no empirical calibration, learning in the baseline, deception, enforcement cost, personality assessment, or proof concerning Freud. Consensus is an efficiency outcome; it is not a welfare score or a test of consent. A refusal is a modeled rejection of a specific differing proposal; a human holding their card is not automatically a refusal. Optional volunteer testing needs separate consent and private intention reports.</p>
<p>Related empirical work found network-dependent convention formation in a different task: <a href="https://arxiv.org/abs/1502.06910">Centola &amp; Baronchelli (2015)</a>. It motivates the question but does not validate this simulator.</p>
<h2>Reproduce and inspect</h2><p><a href="../README.md">Run instructions</a> · <a href="../PROTOCOL.md">Volunteer protocol</a> · <a href="../METHODS.md">Definitions and assumptions</a> · <a href="../SOURCES.md">Sources</a> · <a href="summary.json">Full summary</a> · <a href="condition-summary.csv">Condition table</a> · <a href="paired-contrasts.csv">Matched contrasts</a> · <a href="rounds.csv">Round-level data</a> · <a href="manifest.json">Seeds and file hashes</a></p></main></html>"""
    (out / "report.html").write_text(content, encoding="utf-8")
