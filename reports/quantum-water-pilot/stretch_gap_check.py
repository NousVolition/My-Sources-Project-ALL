"""Compare stretch-band positions and separations in the saved vibration table.

This adds a gap diagnostic; it does not change the historical position criterion.
No gap acceptance threshold is assumed. This script runs no molecular simulation.
"""
import argparse
import csv
import hashlib
import io
import json
from decimal import Decimal
from pathlib import Path

def analyze(path):
    raw = path.read_bytes()
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    groups = {}
    for row in rows:
        group = groups.setdefault(row["isotope"], {})
        if row["mode"] in group:
            raise ValueError("Duplicate mode: " + row["isotope"] + "/" + row["mode"])
        group[row["mode"]] = row
    results = []
    h2o_harmonic_comparison = None
    for isotope, modes in sorted(groups.items()):
        sym = modes["symmetric stretch"]
        anti = modes["antisymmetric stretch"]
        p_sym, p_anti = (Decimal(row["predicted_cm1"]) for row in (sym, anti))
        o_sym, o_anti = (Decimal(row["observed_IR_cm1"]) for row in (sym, anti))
        if not all(x.is_finite() for x in (p_sym, p_anti, o_sym, o_anti)):
            raise ValueError("Nonfinite frequency")
        predicted_gap = p_anti - p_sym
        observed_gap = o_anti - o_sym
        gap_error = predicted_gap - observed_gap
        relative_gap_error = (
            abs(gap_error) / abs(observed_gap) * 100 if observed_gap else None
        )
        position_errors = [
            abs(Decimal(row["predicted_cm1"]) - Decimal(row["observed_IR_cm1"]))
            / abs(Decimal(row["observed_IR_cm1"])) * 100
            for row in modes.values()
        ]
        # The error in a separation is the difference between the two band errors.
        assert gap_error == (p_anti - o_anti) - (p_sym - o_sym)
        if isotope == "H2O":
            harmonic_sym = Decimal("3832")
            harmonic_anti = Decimal("3943")
            harmonic_gap = harmonic_anti - harmonic_sym
            harmonic_gap_error = predicted_gap - harmonic_gap
            h2o_harmonic_comparison = {
                "comparison": "predicted harmonic gap against harmonic reference gap",
                "predicted_gap_cm1": str(predicted_gap),
                "reference_symmetric_stretch_cm1": str(harmonic_sym),
                "reference_antisymmetric_stretch_cm1": str(harmonic_anti),
                "reference_gap_cm1": str(harmonic_gap),
                "signed_gap_error_cm1": str(harmonic_gap_error),
                "absolute_relative_gap_error_percent": str(
                    abs(harmonic_gap_error) / abs(harmonic_gap) * 100
                ),
                "reference": "https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0",
            }
        results.append({
            "isotope": isotope,
            "predicted_gap_cm1": str(predicted_gap),
            "observed_gap_cm1": str(observed_gap),
            "signed_gap_error_cm1": str(gap_error),
            "absolute_relative_gap_error_percent": (
                str(relative_gap_error) if relative_gap_error is not None else None
            ),
            "largest_band_position_error_percent": str(max(position_errors)),
            "historical_5_percent_position_check_passes": max(position_errors) <= 5,
        })
    if h2o_harmonic_comparison is None:
        raise ValueError("The saved comparison must include H2O stretch bands.")
    return {
        "source_csv": path.name,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "gap_definition": "antisymmetric stretch minus symmetric stretch",
        "units": "cm^-1",
        "scope": "New diagnostic of saved predictions; no new physical or Lean run.",
        "gap_acceptance_threshold_percent": None,
        "gap_threshold_note": "No gap threshold was specified in the historical position check.",
        "results": results,
        "h2o_harmonic_comparison": h2o_harmonic_comparison,
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv_path", nargs="?", type=Path,
        default=Path(__file__).with_name("vibration_comparison.csv"),
    )
    args = parser.parse_args()
    print(json.dumps(analyze(args.csv_path), indent=2))
