"""Compare a fresh run with the committed numerical record, within tolerances.

Versions and arbitrary cluster IDs may differ. Membership relations, input
identity, and reported metrics must reproduce. A failure asks for review; it
does not silently refresh the reference.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from study import partition_change


def verify(actual, reference):
    actual, reference = Path(actual), Path(reference)
    a = json.loads((actual / "summary.json").read_text(encoding="utf-8"))
    b = json.loads((reference / "summary.json").read_text(encoding="utf-8"))
    for key in ("source_sha256", "rows", "features", "design", "all_zero_feature_rows"):
        if a[key] != b[key]:
            raise ValueError(f"Changed study input or design: {key}")
    for key, value in b["primary_comparison"].items():
        np.testing.assert_allclose(a["primary_comparison"][key], value, rtol=1e-6, atol=1e-8)
    tables = []
    for directory in (actual, reference):
        with (directory / "assignments.csv").open(encoding="utf-8", newline="") as handle:
            tables.append(list(csv.DictReader(handle)))
    if [r["turn"] for r in tables[0]] != [r["turn"] for r in tables[1]]:
        raise ValueError("Turn order changed.")
    for column in ("counts_cluster", "rates_cluster"):
        labels = [[r[column] for r in table] for table in tables]
        if partition_change(*labels)["changed_pairs"]:
            raise ValueError(f"Primary group membership changed: {column}")
    for filename, id_columns in (("metrics.csv", ("representation", "k", "seed")),
                                  ("stability.csv", ("representation", "k"))):
        with (actual / filename).open(encoding="utf-8", newline="") as handle:
            first = list(csv.DictReader(handle))
        with (reference / filename).open(encoding="utf-8", newline="") as handle:
            second = list(csv.DictReader(handle))
        if len(first) != len(second):
            raise ValueError(f"Row count changed: {filename}")
        for row_a, row_b in zip(first, second):
            for key in row_b:
                if key in id_columns:
                    if row_a[key] != row_b[key]:
                        raise ValueError(f"Changed key: {filename}/{key}")
                else:
                    np.testing.assert_allclose(float(row_a[key]), float(row_b[key]),
                                               rtol=1e-6, atol=1e-8,
                                               err_msg=f"{filename}/{key}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--actual", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    verify(args.actual, args.reference)
    print("Saved numerical record reproduced (tolerance 1e-6 relative / 1e-8 absolute).")
