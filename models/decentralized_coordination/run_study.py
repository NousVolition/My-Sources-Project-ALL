"""Run from this directory: python run_study.py --output results"""
import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import random
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from decentralized_coordination.simulation import (N, Parameters, crossed_mappings, exact_b,
                                                  matched_schedule, seed_for, simulate)
from decentralized_coordination.analysis import influence_analysis, paired_contrasts, summarize


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def write_csv(path, rows):
    if not rows:
        raise ValueError("No rows for " + str(path))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: json.dumps(v, separators=(",", ":")) if isinstance(v, (list, dict))
                         else v for k, v in row.items()} for row in rows)


def compact(result):
    return {k: v for k, v in result.items() if k not in ("trajectory", "events")}


def run_blocks(count, seed, phase, params=None, trace=False):
    schedule = matched_schedule(count, seed_for(seed, phase, "schedule"))
    rounds, traces = [], []
    for spec in schedule:
        for order_index, condition in enumerate(spec["order"]):
            prior_a = spec["order"].index("A") < order_index
            result = simulate(condition, spec["initial"], spec["seat_to_id"],
                              seed_for(seed, phase, spec["block"], condition),
                              params=params, leader_id=spec["leader_id"],
                              former_leader=spec["leader_id"] if condition == "C" and prior_a else None,
                              keep_trace=trace)
            meta = {"phase": phase, "block": spec["block"], "condition_order": spec["order"],
                    "order_index": order_index, "A_before_this_round": prior_a}
            rounds.append({**meta, **compact(result)})
            if trace:
                traces.append({**meta, **result})
    return schedule, rounds, traces


def sensitivity_scenarios():
    cases = {}
    for refusal in (0.0, 0.15, 0.5, 0.9):
        for innovation in (0.0, 0.02, 0.1):
            cases[f"r={refusal:g};i={innovation:g}"] = Parameters(refusal=refusal, innovation=innovation)
    cases.update({"all_hold": Parameters(refusal=1, innovation=0),
                  "synchronous_10s": Parameters(update="sync10"),
                  "slow_opportunities": Parameters(activation=0.08),
                  "fast_opportunities": Parameters(activation=0.6),
                  "heterogeneous_identity_0": Parameters(heterogeneous=True),
                  "former_leader_weight_3": Parameters(former_leader_weight=3),
                  "former_leader_weight_6": Parameters(former_leader_weight=6)})
    return cases


def run_sensitivity(count, seed):
    schedule = matched_schedule(count, seed_for(seed, "sensitivity", "schedule"))
    rows, summaries = [], []
    for name, params in sensitivity_scenarios().items():
        current = []
        for spec in schedule:
            prior_a = spec["order"].index("A") < spec["order"].index("C")
            result = simulate("C", spec["initial"], spec["seat_to_id"],
                              seed_for(seed, "sensitivity", spec["block"]), params=params,
                              leader_id=spec["leader_id"],
                              former_leader=spec["leader_id"] if prior_a else None)
            current.append({"scenario": name, "block": spec["block"], "order": spec["order"],
                            "prior_A": prior_a, **compact(result)})
        summary = summarize(current, name)
        summary["parameters"] = asdict(params)
        shares = [r["former_leader_share"] for r in current if r["former_leader_share"] is not None]
        summary["mean_former_leader_share_after_A"] = sum(shares) / len(shares) if shares else None
        summary["former_leader_share_defined_rounds"] = len(shares)
        rows.extend(current)
        summaries.append(summary)
    # A's perfect obedience is a benchmark assumption, varied explicitly.
    for compliance in (0.7, 0.3):
        name = f"A_compliance={compliance:g}"
        params = Parameters(compliance=compliance)
        current = []
        for spec in schedule:
            result = simulate("A", spec["initial"], spec["seat_to_id"],
                              seed_for(seed, "sensitivity", spec["block"]), params=params,
                              leader_id=spec["leader_id"])
            current.append({"scenario": name, "block": spec["block"], "order": spec["order"],
                            "prior_A": False, **compact(result)})
        summary = summarize(current, name)
        summary.update(parameters=asdict(params), mean_former_leader_share_after_A=None,
                       former_leader_share_defined_rounds=0)
        rows.extend(current)
        summaries.append(summary)
    return rows, summaries


def run_influence(replicates, seed):
    rows, summaries = [], {}
    for topology in ("ring", "star"):
        for profile in ("homogeneous", "heterogeneous"):
            group = []
            for rep in range(replicates):
                rng = random.Random(seed_for(seed, "influence", rep, "start"))
                initial = [rng.randrange(2) for _ in range(N)]
                mappings = crossed_mappings(seed_for(seed, "influence", rep, "mapping"))
                for m, mapping in enumerate(mappings):
                    result = simulate("C", initial, mapping,
                                      seed_for(seed, "influence", rep, m, "choices"),
                                      params=Parameters(heterogeneous=profile == "heterogeneous"),
                                      topology=topology)
                    for seat, identity in enumerate(mapping):
                        group.append({"topology": topology, "profile": profile,
                                      "replicate": rep, "mapping": m, "seat": seat,
                                      "identity": identity, "initial": initial,
                                      "consensus": result["consensus"],
                                      "copied_switches": result["credits_by_seat"][seat]})
            rows.extend(group)
            summaries[topology + "/" + profile] = influence_analysis(group)
    return rows, summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--blocks", type=int, default=4000)
    parser.add_argument("--sensitivity-blocks", type=int, default=1200)
    parser.add_argument("--influence-replicates", type=int, default=240)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    if args.blocks < 2 or args.sensitivity_blocks < 2 or args.influence_replicates < 5:
        parser.error("Need >=2 blocks and >=5 influence replicates")
    out = args.output
    if out.exists() and any(out.iterdir()):
        parser.error("Output directory must be empty; choose a new path to preserve prior results")
    out.mkdir(parents=True, exist_ok=True)
    pilot_schedule, pilot, pilot_traces = run_blocks(4, args.seed, "pilot", trace=True)
    schedule, main_rows, _ = run_blocks(args.blocks, args.seed, "main")
    exact = exact_b()
    write_json(out / "pilot-traces.json", pilot_traces)
    write_csv(out / "pilot-schedule.csv", pilot_schedule)
    write_csv(out / "main-schedule.csv", schedule)
    write_csv(out / "rounds.csv", pilot + main_rows)
    write_csv(out / "exact-B.csv", exact)
    print(f"Pilot: 12 rounds. Main: {len(main_rows)} rounds. Exact B: 512 starts.", flush=True)
    sensitivity, sensitivity_summary = run_sensitivity(args.sensitivity_blocks, args.seed)
    write_csv(out / "sensitivity-rounds.csv", sensitivity)
    write_csv(out / "sensitivity-summary.csv", sensitivity_summary)
    print(f"Sensitivity: {len(sensitivity)} rounds.", flush=True)
    influence, influence_summary = run_influence(args.influence_replicates, args.seed)
    write_csv(out / "influence-observations.csv", influence)
    write_json(out / "influence-summary.json", influence_summary)
    print(f"Influence: {len(influence) // 9} rounds, fully crossed identities and seats.", flush=True)
    summary = {"data_kind": "SIMULATED; no human participants", "seed": args.seed,
               "main": [summarize([r for r in main_rows if r["condition"] == c], c) for c in "ABC"],
               "pilot": [summarize([r for r in pilot if r["condition"] == c], c) for c in "ABC"],
               "paired_contrasts": paired_contrasts(main_rows, seed_for(args.seed, "bootstrap")),
               "B_exact": {"states": len(exact), "successes_120": sum(r["consensus_120"] for r in exact),
                           "fixed_point_basins": sum(r["period"] == 1 for r in exact),
                           "cycle_basins": sum(r["period"] > 1 for r in exact),
                           "maximum_transient_ticks": max(r["transient_ticks"] for r in exact)},
               "sensitivity": sensitivity_summary, "influence": influence_summary}
    write_json(out / "summary.json", summary)
    write_csv(out / "condition-summary.csv", summary["main"])
    write_csv(out / "paired-contrasts.csv", summary["paired_contrasts"])
    import numpy
    manifest = {"data_kind": summary["data_kind"], "seed": args.seed,
                "main_blocks": args.blocks, "pilot_blocks": 4,
                "sensitivity_blocks": args.sensitivity_blocks,
                "influence_replicates": args.influence_replicates,
                "total_simulated_rounds": 12 + len(main_rows) + len(sensitivity) + len(influence) // 9,
                "exact_states_additional": 512, "parameters": asdict(Parameters()),
                "python": platform.python_version(), "numpy": numpy.__version__,
                "randomness": "Python random.Random (MT19937); SHA256 streams; NumPy PCG64 bootstrap",
                "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(Path(__file__).parent.glob("*.py"))},
                "repository_base": "17c14e89940d4c96ec9df315b1a689430d8c5ddf"}
    if not args.no_plots:
        from decentralized_coordination.reporting import make_report
        make_report(out, summary, main_rows, pilot_traces, exact)
        import matplotlib
        manifest["matplotlib"] = matplotlib.__version__
    manifest["data_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(out.iterdir()) if p.suffix in (".csv", ".json")}
    write_json(out / "manifest.json", manifest)
    print(f"Done: {out.resolve()} ({manifest['total_simulated_rounds']} simulated rounds).", flush=True)
    return summary


if __name__ == "__main__":
    main()
