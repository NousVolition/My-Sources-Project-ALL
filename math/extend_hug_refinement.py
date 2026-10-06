"""Extend the saved hug study to 256^3 without repeating completed runs.

--half-step adds the 256^3 time-step control while retaining previous results.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from hug_refinement import compare_fields, run
from navier import PurePythonNavierStokes3D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--half-step", action="store_true")
    args = parser.parse_args()
    path = Path("math/results/hug-refinement.json")
    study = json.loads(path.read_text())
    # Progress reporting does not modify the update or diagnostic calculation.
    original = PurePythonNavierStokes3D.step_with_pressure
    counters = {}
    def progress(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        key = id(self)
        counters[key] = counters.get(key, 0)+1
        if counters[key] % 5 == 0:
            dt = args[0] if args else kwargs['dt']
            print(f"N={self.N} dt={dt:g}: finished step {counters[key]}/{round(.08/dt)}", flush=True)
        return result
    PurePythonNavierStokes3D.step_with_pressure = progress
    try:
        low, low_path = run(128, .001, .08, "scratch/hug-refinement")
        high, high_path = run(256, .001, .08, "scratch/hug-refinement")
        if args.half_step:
            half, half_path = run(256, .0005, .08, "scratch/hug-refinement")
    finally:
        PurePythonNavierStokes3D.step_with_pressure = original
    # Retain the original reviewed results and time-step control. Matching
    # source fingerprints are checked by run() before any cache is reused.
    existing = next(r for r in study["runs"] if r["N"] == 128 and r["dt"] == .001)
    if existing["source_sha256"] != low["source_sha256"]:
        raise RuntimeError("Saved study and cached comparison use different source versions")
    print("Comparing complete fields at matching physical positions", flush=True)
    with np.load(low_path) as a, np.load(high_path) as b:
        comparison = dict(kind="grid", coarse=128, fine=256,
                          initial=compare_fields(a["initial"], b["initial"]),
                          final=compare_fields(a["final"], b["final"]))
    study["runs"] = [r for r in study["runs"] if (r["N"], r["dt"]) != (256, .001)]+[high]
    study["runs"].sort(key=lambda r: (r["N"], -r["dt"]))
    study["spatial_comparisons"] = [p for p in study["spatial_comparisons"]
                                    if (p["coarse"], p["fine"]) != (128, 256)]+[comparison]
    study["spatial_comparisons"].sort(key=lambda p: p["coarse"])
    time_checks = study.get("time_comparisons", [dict(
        grid=study.get("time_comparison_grid",128), steps=[.001,.0005],
        **study["time_comparison"])])
    if args.half_step:
        if half['source_sha256'] != high['source_sha256']:
            raise RuntimeError('Time comparison uses different source versions')
        with np.load(high_path) as a, np.load(half_path) as b:
            time_error = compare_fields(a['final'], b['final'])
            if not np.array_equal(a['initial'], b['initial']):
                raise RuntimeError('Time-step comparison must have identical initial arrays')
        study['runs'] = [r for r in study['runs'] if (r['N'],r['dt']) != (256,.0005)]+[half]
        study['runs'].sort(key=lambda r: (r['N'],-r['dt']))
        time_checks = [c for c in time_checks if c['grid'] != 256]+[
            dict(grid=256,steps=[.001,.0005],**time_error)]
        study['time_comparison'] = time_error
        study['time_comparison_grid'] = 256
        study['extension'] = '256^3 half-step completed; all earlier completed runs reused'
        print(json.dumps(dict(time_comparison_256=time_error),indent=2),flush=True)
    else:
        study.setdefault('time_comparison_grid',128)
        study.setdefault('extension','256^3 spatial run added; completed runs reused')
    study['time_comparisons'] = time_checks
    study['time_comparison_steps'] = [.001,.0005]
    path.write_text(json.dumps(study, indent=2, allow_nan=False)+"\n")
    rows = [row for r in study["runs"] for row in r["rows"]]
    with path.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(comparison, indent=2), flush=True)


if __name__ == "__main__":
    main()
