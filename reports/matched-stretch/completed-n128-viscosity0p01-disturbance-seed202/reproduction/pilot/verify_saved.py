"""Independent verification from immutable saved fields, without rerunning a trajectory."""
import sys
from pathlib import Path

import numpy as np

from diagnostics import ROOT, read, save, sha

sys.path.insert(0, str(ROOT / "source"))
from numerics import Flow


def main():
    p, r = read(ROOT / "protocol.json"), read(ROOT / "result.json")
    f = Flow(128, .01, workers=1)
    folder = Path(p["reference_directory"])
    checks, maximum_error, metric_error = [], 0., 0.
    initial = np.load(ROOT / "initial-disturbed.npy", allow_pickle=False)
    h0 = np.load(folder / "field-000.npy", allow_pickle=False)
    assert np.array_equal(initial, f.perturb(h0, .001, 202))
    assert r["protocol_sha256"] == sha(ROOT / "protocol.json")
    for name, digest in p["pilot_source_hashes"].items():
        assert sha(ROOT / name) == digest
    for name, digest in p["original_source_hashes"].items():
        assert sha(Path(p["original_study"]) / name) == digest
    for i, row in enumerate(r["series"]):
        hp = np.load(ROOT / f"fields/field-{i:03d}.npy", allow_pickle=False)
        hb = np.load(folder / f"field-{i:03d}.npy", allow_pickle=False)
        assert np.isfinite(hp).all() and np.isfinite(hb).all()
        delta = np.fft.irfftn(hp - hb, s=(128,) * 3, axes=(-3, -2, -1))
        physical = np.sqrt(np.sum(delta * delta) * f.dx ** 3) / p["reference_initial_l2"]
        error = abs(physical - row["D_relative_initial_reference_l2"]) / max(physical, 1e-30)
        assert error < 1e-10
        maximum_error = max(maximum_error, error)
        u = np.fft.irfftn(hp, s=(128,) * 3, axes=(-3, -2, -1))
        w = np.fft.irfftn(f.curl(hp), s=(128,) * 3, axes=(-3, -2, -1))
        measured = np.array([np.sqrt(np.sum(w * w, axis=0)).max(),
                             .5 * np.sum(u * u) * f.dx ** 3,
                             .5 * np.sum(w * w) * f.dx ** 3])
        recorded = np.array([row[k] for k in ["Wmax", "energy", "enstrophy"]])
        discrepancy = float(np.max(abs(measured - recorded) / np.maximum(abs(recorded), 1e-30)))
        assert discrepancy < 1e-10
        metric_error = max(metric_error, discrepancy)
        assert np.max(abs(u.mean(axis=(1, 2, 3)) - p["mean_velocity"])) < 1e-10
        assert row["divergence_max"] * f.dx / max(row["peak_speed"], 1e-12) < 1e-10
        assert abs(row["energy_budget_relative_residual"]) < .005
        assert sha(folder / f"field-{i:03d}.npy") == row["reference_field_sha256"]
        checks.append({"t": row["t"], "disturbed_sha256": sha(ROOT / f"fields/field-{i:03d}.npy"),
                       "reference_sha256": row["reference_field_sha256"]})
    save(ROOT / "saved-field-verification.json", {"status": "passed", "through": r["series"][-1]["t"],
         "field_pairs": len(checks), "physical_norm_max_relative_error": maximum_error,
         "physical_W_energy_enstrophy_max_relative_error": metric_error,
         "initial_field_prescribed": True, "fields": checks, "new_simulations": 0})
    print(f"Verified {len(checks)} saved field pairs through {r['series'][-1]['t']:.2f}.")


if __name__ == "__main__":
    main()
