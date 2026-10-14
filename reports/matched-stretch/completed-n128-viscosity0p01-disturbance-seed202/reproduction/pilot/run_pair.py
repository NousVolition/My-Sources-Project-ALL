"""One newly authorized disturbance; reuse the original viscosity-0.01 reference."""
import argparse
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from diagnostics import ROOT, compare, read, save, sha

sys.path.insert(0, str(ROOT / "source"))
from numerics import Flow


def verify_sources(protocol):
    study = Path(protocol["original_study"])
    for name, digest in protocol["original_source_hashes"].items():
        assert sha(study / name) == digest, "Original active source changed"
    assert sha(ROOT / "source/numerics.py") == protocol["original_source_hashes"]["numerics.py"]
    for name, digest in protocol["pilot_source_hashes"].items():
        assert sha(ROOT / name) == digest, "Pilot source changed"
    assert sha(ROOT / "initial-disturbed.npy") == protocol["initial_disturbed_sha256"]
    assert sha(study / "runs/viscosity-n128-nu0.01/field-000.npy") == protocol["reference_initial_sha256"]


def reference_at(protocol, index):
    folder = Path(protocol["reference_directory"])
    while True:
        r = read(folder / "result.json")
        assert r["source_hashes"] == protocol["original_source_hashes"]
        assert r["job"]["nu"] == .01 and r["job"]["n"] == 128
        assert r["dt"] == protocol["dt"] and r["steps_per_output"] == protocol["steps_per_output"]
        if len(r["series"]) > index:
            row = r["series"][index]
            assert abs(row["t"] - index * .01) < 1e-12
            field = folder / f"field-{index:03d}.npy"
            return np.load(field, allow_pickle=False), row, sha(field)
        if r["status"] not in ("running", "complete"):
            raise RuntimeError("Reference failed; preserve evidence and do not substitute a different flow")
        save(ROOT / "progress.json", {"status": "waiting_for_reference", "pid": os.getpid(),
             "requested_t": index * .01, "reference_saved_t": r["series"][-1]["t"], "updated_epoch": time.time()})
        time.sleep(15)


def run():
    p = read(ROOT / "protocol.json")
    verify_sources(p)
    path = ROOT / "result.json"
    if path.exists() and read(path)["status"] == "complete":
        print("Already complete; no rerun.", flush=True)
        return
    lock = ROOT / "controller.lock"
    with lock.open("x", encoding="utf-8") as handle:
        handle.write(str(os.getpid()))
    try:
        f = Flow(128, .01, workers=1)
        h = np.load(ROOT / "initial-disturbed.npy", allow_pickle=False)
        initial = f.observe(h, 0)
        accum = np.zeros(5)
        series, start, max_cfl = [], 0, 0.
        if path.exists():
            old = read(path)
            assert old["protocol_sha256"] == sha(ROOT / "protocol.json")
            assert old["status"] in ("running", "waiting_for_reference"), "A failure must be investigated before resume"
            assert sha(ROOT / old["field"]) == old["field_sha256"]
            h = np.load(ROOT / old["field"], allow_pickle=False)
            accum = np.array(old["accum"])
            series, start, max_cfl = old["series"], old["last_output"] + 1, old["max_speed_cfl"]
        cached = f.rhs(h)
        begin = time.perf_counter()
        for index in range(start, 41):
            if index:
                for sub in range(p["steps_per_output"]):
                    before_W = cached[1][0]
                    h, integral, speed = f.step(h, p["dt"], cached)
                    if not np.isfinite(h).all():
                        np.save(ROOT / "failed-nonfinite-field.npy", h)
                        raise FloatingPointError("Nonfinite field; preserve last saved checkpoint")
                    cached = f.rhs(h)
                    integral[0] = .5 * p["dt"] * (before_W + cached[1][0])
                    accum += integral
                    max_cfl = max(max_cfl, speed * p["dt"] / f.dx)
                    if max_cfl > .5:
                        np.save(ROOT / "failed-cfl-field.npy", h)
                        raise FloatingPointError("Speed CFL exceeded 0.5")
                    if sub % 50 == 0:
                        save(ROOT / "progress.json", {"pid": os.getpid(), "status": "running",
                             "t": (index - 1) * .01 + (sub + 1) * p["dt"], "end": .4,
                             "step": (index - 1) * p["steps_per_output"] + sub + 1,
                             "total_steps": 40 * p["steps_per_output"], "updated_epoch": time.time()})
            row = f.observe(h, index * .01, accum, initial)
            reference, reference_row, reference_hash = reference_at(p, index)
            assert np.isfinite(reference).all()
            row.update(compare(f, h, reference, p["reference_initial_l2"], p["amplitude"]))
            row["reference_field_sha256"] = reference_hash
            row["reference_Wmax"] = reference_row["Wmax"]
            row["reference_energy"] = reference_row["energy"]
            row["reference_resolved_screen"] = reference_row["resolved_screen"]
            row["resolved_screen"] = bool(row["high_band_energy_fraction"] < .001 and
                row["high_band_enstrophy_fraction"] < .01 and
                (row["width_at_global_peak"]["minimum_chord_cells"] or 0) >= 6)
            scaled_div = row["divergence_max"] * f.dx / max(row["peak_speed"], 1e-12)
            mean_error = max(abs(a - b) for a, b in zip(row["mean_velocity"], p["mean_velocity"]))
            failures = []
            if scaled_div > 1e-10: failures.append("Scaled divergence >1e-10")
            if mean_error > 1e-10: failures.append("Mean velocity drift >1e-10")
            if abs(row["energy_budget_relative_residual"]) > .005: failures.append("Energy budget residual >0.5%")
            if row["difference_mean_max"] > 1e-10: failures.append("Difference mean is not preserved")
            if failures: row["hard_failure"] = failures
            filename = f"fields/field-{index:03d}.npy"
            target = ROOT / filename
            tmp = target.with_suffix(".tmp")
            with tmp.open("wb") as handle: np.save(handle, h)
            os.replace(tmp, target)
            series.append(row)
            result = {"status": "failed_check" if failures else ("complete" if index == 40 else "running"),
                "pid": os.getpid(), "job": p["job"], "protocol_sha256": sha(ROOT / "protocol.json"),
                "last_output": index, "field": filename, "field_sha256": sha(target),
                "dt": p["dt"], "accum": accum.tolist(), "series": series,
                "max_speed_cfl": max_cfl, "updated_epoch": time.time(),
                "elapsed_this_process_seconds": time.perf_counter() - begin}
            save(path, result)
            from build_report import build
            build()
            print(f"t={index*.01:.2f} D={row['D_relative_initial_reference_l2']:.9g} "
                  f"amplification={row['amplification']:.7g} status={result['status']}", flush=True)
            if failures:
                raise RuntimeError("; ".join(failures))
            del reference
        verify_sources(p)
        save(ROOT / "progress.json", {"pid": os.getpid(), "status": "complete", "t": .4, "updated_epoch": time.time()})
    except BaseException:
        save(ROOT / "failure.json", {"error": traceback.format_exc(), "pid": os.getpid(), "updated_epoch": time.time()})
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()
