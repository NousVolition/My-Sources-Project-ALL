"""The CI gate must reject missing files and unsafe or stale saved controls."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("hug_package_gate", ROOT / "tools/check_hug_packages.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def inventory(tmp_path):
    (tmp_path / "source.py").write_bytes(b"x = 1\n")
    manifest = {"source.py": hashlib.sha256(b"x = 1\n").hexdigest()}
    (tmp_path / "checksums.json").write_text(json.dumps(manifest))
    return manifest


@pytest.mark.parametrize("change", ["missing_manifest", "empty_manifest", "missing_file", "unlisted", "changed"])
def test_inventory_rejects_incomplete_or_changed_package(tmp_path, change):
    inventory(tmp_path)
    if change == "missing_manifest": (tmp_path / "checksums.json").unlink()
    elif change == "empty_manifest": (tmp_path / "checksums.json").write_text("{}")
    elif change == "missing_file": (tmp_path / "source.py").unlink()
    elif change == "unlisted": (tmp_path / "nested").mkdir(); (tmp_path / "nested/manifest.json").write_text("{}")
    else: (tmp_path / "source.py").write_text("x = 2\n")
    with pytest.raises((ValueError, FileNotFoundError)): gate.check_inventory(tmp_path)


def test_complete_inventory_ignores_only_runtime_cache(tmp_path):
    inventory(tmp_path)
    (tmp_path / "__pycache__").mkdir()
    (tmp_path / "__pycache__/generated.pyc").write_bytes(b"cache")
    assert gate.check_inventory(tmp_path) == 1


@pytest.fixture
def records():
    folder = ROOT / "reports/hug-stress/two-way-data"
    return [gate.read(folder / name) for name in ["protocol.json", "paired.json", "stress.json",
                                                "controls.json", "opening-tests.json", "summary.json"]]


@pytest.mark.parametrize("change", ["missing_pair", "duplicate_pair", "failed_control", "unsafe_bounds",
                                    "budget", "stale_summary", "missing_opening", "event", "released_force", "nonfinite"])
def test_control_gate_rejects_incomplete_unsafe_and_stale_records(records, change):
    protocol, paired, stress, controls, opening, summary = copy.deepcopy(records)
    if change == "missing_pair": paired.pop()
    elif change == "duplicate_pair": paired[-1] = copy.deepcopy(paired[0])
    elif change == "failed_control": controls[0]["independent_error"] = protocol["control_target"]
    elif change == "unsafe_bounds": paired[0]["metrics"]["memory_bounds"] = False
    elif change == "budget": paired[0]["metrics"]["scaled_activity_budget_error"] = 1e-7
    elif change == "stale_summary": summary["maximum_scaled_budget_residual"] += 1e-6
    elif change == "missing_opening": opening["cases"].pop()
    elif change == "event": opening["cases"][0]["event_error"] = 1e-6
    elif change == "released_force": opening["cases"][0]["resistance_after_release"] = 1e-10
    else: opening["cases"][0]["scaled_error"] = float("nan")
    with pytest.raises(ValueError): gate.check_two_way_records(protocol, paired, stress, controls, opening, summary)


def test_current_saved_control_records_pass(records):
    assert gate.check_two_way_records(*records)["opening_controls"] == 12
