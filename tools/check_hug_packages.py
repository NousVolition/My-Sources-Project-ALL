"""Read-only integrity and saved-control gates for the published HUG packages.

No integration, report rebuilding, or numerical-source changes occur here.
Historical stress failures retain the interpretation in verify_package.py.
"""
from pathlib import Path
import hashlib
import json
import math
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_inventory(root):
    root = root.resolve()
    manifest = read(root / "checksums.json")
    require(isinstance(manifest, dict) and bool(manifest), "Empty or invalid package manifest")
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*")
              if p.is_file() and not any(part in {"__pycache__", ".pytest_cache"}
                                        for part in p.relative_to(root).parts)
              and p.relative_to(root).as_posix() != "checksums.json"}
    require(set(manifest) == actual, "Manifest coverage differs: "
            f"unlisted={sorted(actual-set(manifest))}, missing={sorted(set(manifest)-actual)}")
    for name, digest in manifest.items():
        path = (root / name).resolve()
        require(path.is_relative_to(root), "Manifest path escapes its package")
        require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest),
                "Malformed SHA-256 digest")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"Hash mismatch: {name}")
    return len(manifest)


def finite_record(value):
    if isinstance(value, dict):
        return all(finite_record(v) for v in value.values())
    if isinstance(value, list):
        return all(finite_record(v) for v in value)
    if isinstance(value, float):
        return math.isfinite(value)
    return True


def check_two_way_records(protocol, paired, stress, controls, opening, summary):
    require(all(finite_record(x) for x in [protocol, paired, stress, controls, opening, summary]),
            "Nonfinite saved control record")
    target = protocol["control_target"]
    require(0 < target <= 1e-5, "Invalid saved control target")
    identities = {(r["seed"], r["strength"]) for r in paired}
    require(len(paired) == 12 and identities == {(s, g) for s in range(3)
                                                for g in protocol["feedback_strengths"]},
            "Paired configuration identities differ")
    require(len(stress) == 12 and len({r["name"] for r in stress}) == 12,
            "Stress configuration count differs")
    require(len(controls) == 7 and len({r["name"] for r in controls}) == 7,
            "Independent control count differs")
    for row in controls:
        require(row["passed"] is True and 0 <= row["independent_error"] < target
                and 0 <= row["refinement_error"] < target, "Independent control failed")
    metrics = [r["metrics"] for r in paired + stress] + [r["refined"] for r in paired]
    for row in metrics:
        require(all(row[k] is True for k in ["memory_bounds", "pressure_bounds", "passive_feedback"]),
                "Saved bounds or passivity gate failed")
        require(0 <= row["scaled_activity_budget_error"] < 1e-7, "Saved activity budget gate failed")
    maximum = max(r["scaled_activity_budget_error"] for r in metrics)
    require(summary["controls_passed"] is True and summary["all_bounds_and_passivity"] is True,
            "Summary gate failed")
    require(summary["paired_configurations"] == 12 and summary["paired_solver_runs"] == 24
            and summary["stress_configurations"] == 12 and summary["short_three_solver_controls"] == 7,
            "Summary configuration counts differ")
    require(abs(summary["maximum_scaled_budget_residual"] - maximum) < 1e-15,
            "Stale summary budget maximum")
    cases = opening["cases"]
    require(len(cases) == 12 and len({r["name"] for r in cases}) == 12
            and opening["all_passed"] is True, "Opening control count or gate failed")
    for row in cases:
        require(row["passed"] is True and 0 <= row["scaled_error"] < target
                and 0 <= row["event_error"] < 1e-6 and row["resistance_after_release"] == 0,
                "Opening control failed")
    require(abs(opening["maximum_error"] - max(r["scaled_error"] for r in cases)) < 1e-15
            and abs(opening["max_event_error"] - max(r["event_error"] for r in cases)) < 1e-15,
            "Stale opening error maximum")
    return {"paired": 12, "stress": 12, "independent_controls": 7, "opening_controls": 12}


def check_two_way(root):
    data = root / "two-way-data"
    protocol = read(data / "protocol.json")
    for name, digest in protocol["sources"].items():
        require(hashlib.sha256((root / name).read_bytes()).hexdigest() == digest,
                f"Two-way source drift: {name}")
    result = check_two_way_records(protocol, read(data / "paired.json"), read(data / "stress.json"),
                                   read(data / "controls.json"), read(data / "opening-tests.json"),
                                   read(data / "summary.json"))
    opening = read(data / "opening-tests.json")
    require(hashlib.sha256((root / "test_two_way_opening.py").read_bytes()).hexdigest()
            == opening["source_sha256"], "Opening-control source drift")
    paths = sorted(data.glob("*.npz"))
    require(bool(paths), "No saved two-way arrays")
    for path in paths:
        with np.load(path, allow_pickle=False) as arrays:
            require(bool(arrays.files) and all(np.isfinite(arrays[k]).all() for k in arrays.files),
                    f"Nonfinite or empty saved array: {path.name}")
            if "time" in arrays:
                t = arrays["time"]
                require(t.ndim == 1 and t.size > 1 and np.all(np.diff(t) > 0),
                        f"Invalid saved time axis: {path.name}")
    return dict(result, array_archives=len(paths))


def main():
    reports = ROOT / "reports"
    result = {name: {"hashed_files": check_inventory(reports / name)}
              for name in ["hug-dynamics", "hug-stress"]}
    result["hug-stress"]["two_way"] = check_two_way(reports / "hug-stress")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
