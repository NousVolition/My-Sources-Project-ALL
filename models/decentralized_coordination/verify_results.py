"""Check manifest integrity and exact data reproduction (plots are excluded)."""
import argparse
import hashlib
import json
from pathlib import Path


def integrity(folder):
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["data_sha256"].items():
        actual = hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"Data hash mismatch: {folder / name}")
    return manifest


def verify(actual, expected):
    actual_manifest, expected_manifest = integrity(actual), integrity(expected)
    for field in ("seed", "main_blocks", "pilot_blocks", "sensitivity_blocks", "influence_replicates", "parameters"):
        if actual_manifest[field] != expected_manifest[field]:
            raise AssertionError(f"Different design: {field}")
    for name in expected_manifest["data_sha256"]:
        if actual_manifest["data_sha256"].get(name) != expected_manifest["data_sha256"][name]:
            raise AssertionError(f"Reproduction differs: {name}")
    return len(expected_manifest["data_sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--actual", type=Path, required=True)
    parser.add_argument("--expected", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    print(f"PASS: {verify(args.actual, args.expected)} generated data files match exactly; plot bytes excluded.")
