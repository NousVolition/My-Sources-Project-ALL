"""Audit how reply-length adjustment changes clustering of saved phrase counts.

Run from the repository root:
    python language/clustering/study.py --output scratch/clustering
The existing CSV is the input artifact; its missing generator is not recreated.
"""
import argparse
from collections import Counter
import csv
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy
import sklearn
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_samples
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits


SOURCE = Path(__file__).resolve().parents[1] / "data/ai_language_turns_3.csv"
FEATURES = (
    "validation", "apology", "institutional", "clinical_pr", "machine_talk",
    "supremacy_term", "exit_menu", "safety_story",
)
REPRESENTATIONS = ("counts", "per_100_words")
KS = (2, 3, 4, 5, 6)
SEEDS = tuple(range(10))
PRIMARY_K = 3


def load_rows(path):
    """Require valid saved counts; never silently impute malformed observations."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"turn", "primary_mode", "ai_words", *FEATURES}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("The source table is empty.")
    seen = set()
    for row in rows:
        try:
            turn = int(row["turn"])
            values = [float(row[name]) for name in ("ai_words", *FEATURES)]
        except (TypeError, ValueError) as exc:
            raise ValueError("Turn IDs and counts must be numeric.") from exc
        if turn < 1 or turn in seen:
            raise ValueError("Turn IDs must be positive and unique.")
        seen.add(turn)
        if not all(math.isfinite(v) and v >= 0 and v.is_integer() for v in values):
            raise ValueError(f"Turn {turn}: counts must be finite nonnegative integers.")
        if values[0] <= 0:
            raise ValueError(f"Turn {turn}: ai_words must be positive for rate normalization.")
        if not row["primary_mode"] or not row["primary_mode"].strip():
            raise ValueError(f"Turn {turn}: missing descriptive label.")
    return rows


def feature_matrix(rows, representation):
    """Select only the eight count columns; labels and prose never enter X."""
    if representation not in REPRESENTATIONS:
        raise ValueError(f"Unknown representation: {representation}")
    counts = np.asarray([[float(row[name]) for name in FEATURES] for row in rows])
    words = np.asarray([float(row["ai_words"]) for row in rows])
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Counts must be finite and nonnegative.")
    if not np.isfinite(words).all() or (words <= 0).any():
        raise ValueError("Word counts must be finite and positive.")
    values = counts if representation == "counts" else counts * 100.0 / words[:, None]
    scaler = StandardScaler()
    scaled = scaler.fit_transform(values)
    return values, scaled, scaler


def fit_partition(x, k, seed):
    if not 2 <= k < len(x) or len(np.unique(x, axis=0)) < k:
        raise ValueError("Require 2 <= K < row count and at least K distinct vectors.")
    # One thread reduces platform-dependent reduction differences on this tiny table.
    with threadpool_limits(limits=1):
        model = KMeans(n_clusters=k, init="k-means++", n_init=20, random_state=seed,
                       algorithm="lloyd", max_iter=300, tol=1e-4).fit(x)
        sil = silhouette_samples(x, model.labels_, metric="euclidean")
    sizes = np.bincount(model.labels_, minlength=k)
    stats = dict(inertia=float(model.inertia_), silhouette=float(sil.mean()),
                 min_cluster_size=int(sizes.min()), max_cluster_size=int(sizes.max()),
                 singleton_clusters=int((sizes == 1).sum()),
                 negative_silhouette_rows=int((sil < 0).sum()))
    return model, sil, stats


def partition_change(a, b):
    """Compare same-group relations, so arbitrary cluster numbers cannot mislead."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.ndim != 1 or len(a) < 2:
        raise ValueError("Two equally sized label vectors with at least two rows are needed.")
    changed = (a[:, None] == a[None, :]) != (b[:, None] == b[None, :])
    per_row = changed.sum(axis=1)
    pairs = int(np.triu(changed, k=1).sum())
    return dict(adjusted_rand_index=float(adjusted_rand_score(a, b)),
                changed_pairs=pairs, total_pairs=len(a) * (len(a) - 1) // 2,
                per_row_changed=per_row.tolist())


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_study(source=SOURCE, output=Path("scratch/clustering"), make_plot=True):
    source, output = Path(source), Path(output)
    rows = load_rows(source)
    output.mkdir(parents=True, exist_ok=True)
    raw = np.asarray([[float(row[name]) for name in FEATURES] for row in rows])
    label_counts = Counter(row["primary_mode"] for row in rows)
    summary = {
        "question": "How does adjusting phrase counts for reply length change grouping?",
        "source_file": source.name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "rows": len(rows), "features": list(FEATURES),
        "descriptive_labels": len(label_counts),
        "singleton_descriptive_labels": sum(n == 1 for n in label_counts.values()),
        "all_zero_feature_rows": int((raw == 0).all(axis=1).sum()),
        "feature_nonzero_rows": {name: int((raw[:, i] > 0).sum()) for i, name in enumerate(FEATURES)},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scikit_learn": sklearn.__version__, "scipy": scipy.__version__},
        "design": {"K_values": list(KS), "seeds": list(SEEDS), "n_init": 20,
                   "algorithm": "lloyd", "initialization": "k-means++", "max_iter": 300,
                   "tol": 1e-4, "primary_K": PRIMARY_K, "primary_seed": 0,
                   "scaler": "StandardScaler fitted separately for each representation",
                   "selection": "K=3 fixed as an illustrative comparison before seeing results",
                   "labels_used_for_fitting": False,
                   "scope": "Exploratory full-table analysis; no held-out performance claim"},
        "representations": {},
    }
    metrics, stability, profiles = [], [], []
    primary = {}
    for representation in REPRESENTATIONS:
        values, x, scaler = feature_matrix(rows, representation)
        summary["representations"][representation] = {
            "unique_vectors": len(np.unique(values, axis=0)),
            "scaler_mean": scaler.mean_.tolist(), "scaler_scale": scaler.scale_.tolist(),
        }
        for k in KS:
            partitions = []
            for seed in SEEDS:
                model, sil, stats = fit_partition(x, k, seed)
                metrics.append(dict(representation=representation, k=k, seed=seed, **stats))
                partitions.append(model.labels_)
                if k == PRIMARY_K and seed == 0:
                    primary[representation] = (model.labels_, sil)
                    summary["representations"][representation]["primary"] = dict(
                        **stats, sizes=np.bincount(model.labels_, minlength=k).tolist())
                    for cluster in range(k):
                        members = values[model.labels_ == cluster]
                        profiles.append(dict(representation=representation, cluster=cluster,
                                             size=len(members), **dict(zip(FEATURES, members.mean(axis=0)))))
            scores = [adjusted_rand_score(a, b) for a, b in combinations(partitions, 2)]
            group = [m for m in metrics if m["representation"] == representation and m["k"] == k]
            stability.append(dict(representation=representation, k=k,
                                  seed_ari_mean=float(np.mean(scores)), seed_ari_min=float(min(scores)),
                                  silhouette_mean=float(np.mean([m["silhouette"] for m in group])),
                                  silhouette_min=float(min(m["silhouette"] for m in group)),
                                  silhouette_max=float(max(m["silhouette"] for m in group)),
                                  min_cluster_size=min(m["min_cluster_size"] for m in group)))
    a, a_sil = primary["counts"]
    b, b_sil = primary["per_100_words"]
    change = partition_change(a, b)
    summary["primary_comparison"] = {k: v for k, v in change.items() if k != "per_row_changed"}
    summary["primary_seed_stability"] = [s for s in stability if s["k"] == PRIMARY_K]
    assignments = []
    for i, row in enumerate(rows):
        assignments.append(dict(turn=int(row["turn"]), ai_words=int(row["ai_words"]),
                                counts_cluster=int(a[i]), rates_cluster=int(b[i]),
                                counts_silhouette=float(a_sil[i]), rates_silhouette=float(b_sil[i]),
                                changed_group_relations=change["per_row_changed"][i],
                                compared_rows=len(rows) - 1))
    # Existing text labels are used only in this post-fit descriptive cross-tabulation.
    cross_tab = Counter((rep, int(labels[i]), row["primary_mode"])
                        for rep, (labels, _) in primary.items() for i, row in enumerate(rows))
    label_table = [dict(representation=rep, cluster=cluster, existing_label=label, count=count)
                   for (rep, cluster, label), count in sorted(cross_tab.items())]
    review = sorted(assignments, key=lambda row: (-row["changed_group_relations"], row["turn"]))[:8]
    summary["review_queue_turns"] = [row["turn"] for row in review]
    summary["review_queue_status"] = "Automatically prioritized examples; no human adjudication recorded"
    for name, data in (("metrics.csv", metrics), ("stability.csv", stability),
                       ("assignments.csv", assignments), ("cluster-profiles.csv", profiles),
                       ("label-comparison.csv", label_table), ("review-queue.csv", review)):
        save_csv(output / name, data)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if make_plot:
        plot_results(output, stability, a, b, summary)
    return summary


def plot_results(output, stability, a, b, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                                 "svg.hashsalt": "nous-feature-study-v1"})
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.5), layout="constrained")
    colors = {"counts": "#156b64", "per_100_words": "#aa553c"}
    titles = {"counts": "Phrase counts", "per_100_words": "Counts per 100 words"}
    for rep in REPRESENTATIONS:
        group = [row for row in stability if row["representation"] == rep]
        means = np.array([row["silhouette_mean"] for row in group])
        low = means - np.array([row["silhouette_min"] for row in group])
        high = np.array([row["silhouette_max"] for row in group]) - means
        axes[0].errorbar(KS, means, yerr=[low, high], marker="o", capsize=4,
                         label=titles[rep], color=colors[rep], linewidth=2)
    axes[0].set(title="Group separation across K", xlabel="Number of groups (K)",
                ylabel="Mean silhouette; bars show seed range", xticks=KS, ylim=(-0.05, 1))
    axes[0].legend(loc="lower right", frameon=False)
    axes[0].grid(axis="y", alpha=0.2)
    matrix = np.zeros((PRIMARY_K, PRIMARY_K), dtype=int)
    np.add.at(matrix, (a, b), 1)
    axes[1].imshow(matrix, cmap="BuGn", vmin=0, vmax=matrix.max())
    for i in range(PRIMARY_K):
        for j in range(PRIMARY_K):
            axes[1].text(j, i, str(matrix[i, j]), ha="center", va="center",
                         color="white" if matrix[i, j] > matrix.max() * 0.6 else "#172b2a", fontsize=15)
    axes[1].set(title=f"Same {summary['rows']} replies, different representation\nK={PRIMARY_K}; adjusted Rand index {summary['primary_comparison']['adjusted_rand_index']:.3f}",
                xlabel="Per-100-word group (arbitrary ID)", ylabel="Count group (arbitrary ID)",
                xticks=range(PRIMARY_K), yticks=range(PRIMARY_K))
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Representation changes the grouping", fontsize=17, fontweight="bold")
    fig.savefig(output / "study-overview.svg", metadata={"Date": None})
    fig.savefig(output / "study-overview.png", dpi=160, metadata={"Software": "Nous Volition study"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=Path("scratch/clustering"))
    args = parser.parse_args()
    summary = run_study(args.input, args.output)
    print(json.dumps({"rows": summary["rows"], "labels": summary["descriptive_labels"],
                      **summary["primary_comparison"], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
