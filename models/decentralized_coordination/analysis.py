"""Descriptive Monte Carlo summaries, matched contrasts, and influence checks."""
import math
import numpy as np


def wilson(successes, total):
    if total <= 0:
        return [None, None]
    p, z = successes / total, 1.959963984540054
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total**2)) / den
    return [max(0.0, center - half), min(1.0, center + half)]


def summarize(rows, label):
    n = len(rows)
    successes = sum(r["consensus"] for r in rows)
    times = [r["consensus_time"] for r in rows if r["consensus"]]
    hhis = [r["copy_hhi"] for r in rows if r["copy_hhi"] is not None]
    proposals = sum(r["differing_proposals"] for r in rows)
    return {"label": label, "rounds": n, "successes": successes,
            "consensus_proportion": successes / n,
            "consensus_ci95": wilson(successes, n),
            "mean_capped_seconds": float(np.mean([r["capped_time"] for r in rows])),
            "median_seconds_successes_only": float(np.median(times)) if times else None,
            "mean_switches": float(np.mean([r["switches"] for r in rows])),
            "mean_refusals": float(np.mean([r["refusals"] for r in rows])),
            "refusal_per_differing_proposal": sum(r["refusals"] for r in rows) / proposals
            if proposals else None,
            "mean_initiatives": float(np.mean([r["initiatives"] for r in rows])),
            "mean_adjacent_tick_reversals": float(np.mean([r["adjacent_tick_reversals"] for r in rows])),
            "proportion_state_revisit": float(np.mean([r["state_revisits"] > 0 for r in rows])),
            "mean_copy_hhi_when_defined": float(np.mean(hhis)) if hhis else None,
            "copy_hhi_defined_rounds": len(hhis)}


def paired_contrasts(rows, seed, draws=1500):
    """Resample whole matched starting-state blocks, not individual actors."""
    groups = {c: sorted([r for r in rows if r["condition"] == c], key=lambda r: r["block"])
              for c in "ABC"}
    n = len(groups["A"])
    if any([r["block"] for r in groups[c]] != [r["block"] for r in groups["A"]] for c in "BC"):
        raise ValueError("Paired contrasts require the same blocks")
    rng = np.random.default_rng(seed)
    out = []
    for left, right in (("A", "B"), ("C", "B"), ("A", "C")):
        for measure in ("consensus", "capped_time", "switches"):
            delta = np.array([float(a[measure]) - float(b[measure])
                              for a, b in zip(groups[left], groups[right])])
            boots = np.array([delta[rng.integers(0, n, n)].mean() for _ in range(draws)])
            out.append({"contrast": left + " minus " + right, "measure": measure,
                        "mean_difference": float(delta.mean()),
                        "ci95": np.quantile(boots, [0.025, 0.975]).tolist(),
                        "resampling_unit": "matched block", "bootstrap_draws": draws})
    return out


def influence_analysis(rows):
    """Balanced 9x9 seat/identity crossed descriptive effects and held-out MSE.

    Repetitions, including their nine mappings, are kept together in train/test.
    These are model diagnostics, not human leadership effect estimates.
    """
    reps = sorted({r["replicate"] for r in rows})
    cut = max(1, int(len(reps) * 0.8))
    train_reps = set(reps[:cut])
    y = np.array([r["copied_switches"] for r in rows])
    ids = np.array([r["identity"] for r in rows])
    seats = np.array([r["seat"] for r in rows])
    train = np.array([r["replicate"] in train_reps for r in rows])
    grand = y.mean()
    by_id = np.array([y[ids == i].mean() for i in range(9)])
    by_seat = np.array([y[seats == i].mean() for i in range(9)])
    total_variance = y.var()
    fractions = {"identity": float(((by_id - grand)**2).mean() / total_variance) if total_variance else 0,
                 "position": float(((by_seat - grand)**2).mean() / total_variance) if total_variance else 0}
    fractions["unexplained"] = max(0, 1 - sum(fractions.values()))
    mean_train = y[train].mean()
    id_effect = np.array([y[train & (ids == i)].mean() - mean_train for i in range(9)])
    seat_effect = np.array([y[train & (seats == i)].mean() - mean_train for i in range(9)])
    baseline_mse = ((y[~train] - mean_train)**2).mean() if (~train).any() else 0
    gains = {}
    for name, effect in (("identity", id_effect[ids]), ("position", seat_effect[seats]),
                         ("both", id_effect[ids] + seat_effect[seats])):
        mse = ((y[~train] - mean_train - effect[~train])**2).mean() if (~train).any() else 0
        gains[name] = float(1 - mse / baseline_mse) if baseline_mse else None
    # Ties receive equal mass; never award a winner to the lowest seat number.
    rounds = {}
    for row in rows:
        rounds.setdefault((row["replicate"], row["mapping"]), []).append(row)
    distributions = {}
    for key, group in rounds.items():
        best = max(r["copied_switches"] for r in group)
        if best == 0:
            continue
        winners = [r for r in group if r["copied_switches"] == best]
        distributions[key] = {axis: {r[axis]: 1 / len(winners) for r in winners}
                              for axis in ("identity", "seat")}
    persistence = {"identity": [], "seat": []}
    for rep in reps:
        for mapping in range(1, 9):
            previous, current = distributions.get((rep, mapping - 1)), distributions.get((rep, mapping))
            if previous is None or current is None:
                continue
            for axis in persistence:
                persistence[axis].append(sum(v * current[axis].get(i, 0) for i, v in previous[axis].items()))
    return {"mean_by_identity": by_id.tolist(), "mean_by_seat": by_seat.tolist(),
            "variance_fraction": fractions, "heldout_mse_improvement": gains,
            "train_replicates": cut, "test_replicates": len(reps) - cut,
            "top_copy_source_persistence": {k: float(np.mean(v)) if v else None for k, v in persistence.items()},
            "persistence_pairs": len(persistence["identity"]),
            "persistence_reference": "1/9 is only an exchangeable, independent reference; matched rounds are dependent"}
