"""Map saved annotations in turn order, without exporting conversation prose.

This reads the existing third-version table; it does not re-label the replies.
Run: python language/conversation_map.py --output scratch/conversation-map
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

from clustering.study import FEATURES, SOURCE, load_rows


LABELS = ("Validation", "Apology", "Institutional", "Clinical / PR",
          "Machine talk", "Supremacy term", "Exit menu", "Safety story")


def map_data(source):
    """Explicit allowlist excludes user text and assistant openings."""
    source = Path(source)
    rows = sorted(load_rows(source), key=lambda row: int(row["turn"]))
    turns = [dict(turn=int(row["turn"]), words=int(row["ai_words"]),
                  label=row["primary_mode"],
                  counts=[int(row[key]) for key in FEATURES]) for row in rows]
    labels = Counter(row["label"] for row in turns)
    co_occurrence = []
    for i, first in enumerate(FEATURES):
        for j in range(i + 1, len(FEATURES)):
            both = [row["turn"] for row in turns if row["counts"][i] and row["counts"][j]]
            if both:
                co_occurrence.append(dict(first=first, second=FEATURES[j], turns=both))
    return dict(
        schema_version=1, conversation_id="saved-language-map-3",
        annotation_version="saved-v3-unverified-generator",
        source_file=source.name, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        features=[dict(key=k, label=label) for k, label in zip(FEATURES, LABELS)],
        turns=turns,
        summary=dict(replies=len(turns), descriptive_labels=len(labels),
                     no_recorded_matches=[row["turn"] for row in turns if not any(row["counts"])],
                     feature_turns={key: [row["turn"] for row in turns if row["counts"][i]]
                                    for i, key in enumerate(FEATURES)},
                     co_occurrence=co_occurrence),
    )


def plot_map(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.hashsalt": "nous-conversation-map-v1"})
    rows = data["turns"]
    turns = np.array([r["turn"] for r in rows])
    values = np.array([r["counts"] for r in rows]).T
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                     gridspec_kw={"height_ratios": [1, 2.4]}, layout="constrained")
    top.bar(turns, [r["words"] for r in rows], width=.75, color="#26776e")
    top.set_ylabel("Reply length\n(words)")
    top.spines[["top", "right"]].set_visible(False)
    top.grid(axis="y", alpha=.15)
    # A pcolormesh preserves missing turn IDs as gaps rather than squeezing time.
    full = np.zeros((len(FEATURES), int(turns.max() - turns.min() + 1)))
    full[:, turns - turns.min()] = values > 0
    bottom.pcolormesh(np.arange(turns.min() - .5, turns.max() + 1.5),
                      np.arange(len(FEATURES) + 1) - .5, full,
                      cmap=ListedColormap(["#f2f4f3", "#26776e"]), vmin=0, vmax=1,
                      edgecolors="white", linewidth=.8)
    bottom.invert_yaxis()
    bottom.set(yticks=range(len(FEATURES)), yticklabels=LABELS,
               xlabel="Reply number in the saved conversation", ylabel="Saved category")
    bottom.set_xlim(turns.min() - .5, turns.max() + .5)
    top.set_title("Conversation without the reply text", loc="left", fontsize=18, pad=14)
    fig.supxlabel("Filled cell = recorded count above zero; pale cell = no recorded match.\n"
                  "Existing annotations are preserved; overlapping categories are allowed.", fontsize=10)
    for suffix in ("svg", "png"):
        options = {"metadata": {"Date": None}} if suffix == "svg" else {"dpi": 160}
        fig.savefig(output / f"conversation-map.{suffix}", **options)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=Path("scratch/conversation-map"))
    args = parser.parse_args()
    data = map_data(args.input)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "map-data.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    plot_map(data, args.output)
    print(json.dumps({"replies": data["summary"]["replies"],
                      "no_recorded_matches": len(data["summary"]["no_recorded_matches"]),
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
