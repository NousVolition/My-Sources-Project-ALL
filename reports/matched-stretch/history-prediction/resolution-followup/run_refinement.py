"""Run the bounded follow-up without changing the original solver or pilot."""
from pathlib import Path
import argparse
import datetime
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from simulate import run, sha, save_json
import numpy as np


def validate(path, expected):
    meta = json.loads(path.with_suffix(".json").read_text())
    assert meta["config"] == expected and sha(path) == meta["sha256"]
    assert meta["max_cfl"] < 0.5
    assert max(abs(r["energy_budget_residual"]) for r in meta["diagnostics"]) < .005
    assert max(r["divergence_max"] for r in meta["diagnostics"]) < 1e-10
    with np.load(path, allow_pickle=False) as arrays:
        assert all(np.isfinite(arrays[k]).all() for k in arrays.files)
        assert np.allclose(arrays["time"], np.arange(49)*.025)
        assert arrays["positions"].shape == (49, 512, 3)
    return {"file": path.name, "sha256": meta["sha256"], "seed": meta["seed"],
            "config": meta["config"], "seconds": meta["seconds"], "checked": True}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true")
    args = ap.parse_args()
    p = json.loads((HERE/"protocol.json").read_text())
    manifest = {"protocol_sha256": sha(HERE/"protocol.json"),
                "runner_sha256": sha(Path(__file__)), "records": []}
    jobs = [(s, {**p["base"], "n": p["new_grid"]}) for s in p["seeds"]]
    if args.pilot:
        jobs = jobs[:1]
    else:
        jobs += [(s, {**p["base"], "n": n, "dt": p["half_step_dt"]})
                 for s in p["half_step_seeds"] for n in p["half_step_grids"]]
    for seed, config in jobs:
        path = run(seed, config, HERE/"recorded-data")
        manifest["records"].append(validate(path, config))
        manifest["updated_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest["complete"] = len(manifest["records"]) == len(jobs)
        save_json(HERE/("pilot.json" if args.pilot else "execution.json"), manifest)
    print(json.dumps({"verified_simulations": len(jobs), "pilot": args.pilot}), flush=True)


if __name__ == "__main__":
    main()
