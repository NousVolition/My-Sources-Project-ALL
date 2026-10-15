"""One authorized viscosity-0.1 pair; original queues and sources stay intact."""
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from diagnostics import ROOT, compare, read, save, sha

sys.path.insert(0, str(ROOT / "source"))
from numerics import Flow


def verify_sources(p):
    for name, digest in p["original_source_hashes"].items():
        assert sha(Path(p["original_study"]) / name) == digest, "Original source changed"
    assert sha(ROOT / "source/numerics.py") == p["original_source_hashes"]["numerics.py"]
    for name, digest in p["source_hashes"].items():
        assert sha(ROOT / name) == digest, "Pair source changed"
    for name, digest in p["initial_hashes"].items():
        assert sha(ROOT / name) == digest, "Prescribed initial field changed"


def screen(row):
    return bool(row["high_band_energy_fraction"] < .001 and
                row["high_band_enstrophy_fraction"] < .01 and
                (row["width_at_global_peak"]["minimum_chord_cells"] or 0) >= 6)


def check(row, f, p):
    failures = []
    if row["divergence_max"] * f.dx / max(row["peak_speed"], 1e-12) > 1e-10:
        failures.append("scaled divergence exceeds 1e-10")
    if max(abs(a-b) for a, b in zip(row["mean_velocity"], p["mean_velocity"])) > 1e-10:
        failures.append("mean velocity drift exceeds 1e-10")
    if abs(row["energy_budget_relative_residual"]) > .005:
        failures.append("energy budget residual exceeds 0.5%")
    return failures


def save_field(path, value):
    tmp = path.with_suffix(".tmp")
    with tmp.open("wb") as handle:
        np.save(handle, value)
    os.replace(tmp, path)


def run():
    p = read(ROOT / "protocol.json")
    verify_sources(p)
    result_path = ROOT / "result.json"
    if result_path.exists() and read(result_path)["status"] == "complete":
        print("Already complete; no rerun.", flush=True)
        return
    lock = ROOT / "controller.lock"
    with lock.open("x", encoding="utf-8") as handle:
        handle.write(str(os.getpid()))
    try:
        f = Flow(p["grid"], p["viscosity"], workers=1)
        fields = {name: np.load(ROOT / f"initial-{name}.npy", allow_pickle=False)
                  for name in ("reference", "disturbed")}
        initial = {name: f.observe(h, 0) for name, h in fields.items()}
        initial_difference_energy = .5 * f.inner(fields["disturbed"] - fields["reference"],
                                                fields["disturbed"] - fields["reference"])
        accum = {name: np.zeros(5) for name in fields}
        max_cfl = {name: 0. for name in fields}
        integrated_difference_net = 0.
        rows, start = [], 0
        if result_path.exists():
            old = read(result_path)
            assert old["status"] == "running", "Investigate a recorded failure before resume"
            assert old["protocol_sha256"] == sha(ROOT / "protocol.json")
            for name in fields:
                assert sha(ROOT / old["fields"][name]) == old["field_hashes"][name]
                fields[name] = np.load(ROOT / old["fields"][name], allow_pickle=False)
                accum[name] = np.array(old["accum"][name])
            rows, start = old["series"], old["last_output"] + 1
            max_cfl = old["max_speed_cfl"]
            integrated_difference_net = old["integrated_difference_net"]
        cached = {name: f.rhs(h) for name, h in fields.items()}
        begin = time.perf_counter()
        for index in range(start, 41):
            if index:
                for sub in range(p["steps_per_output"]):
                    delta = fields["disturbed"] - fields["reference"]
                    before_rate = f.inner(delta, cached["disturbed"][0] - cached["reference"][0])
                    for name in fields:
                        before_W = cached[name][1][0]
                        fields[name], integral, speed = f.step(fields[name], p["dt"], cached[name])
                        if not np.isfinite(fields[name]).all():
                            np.save(ROOT / f"failed-nonfinite-{name}.npy", fields[name])
                            raise FloatingPointError(f"Nonfinite {name} field")
                        cached[name] = f.rhs(fields[name])
                        integral[0] = .5 * p["dt"] * (before_W + cached[name][1][0])
                        accum[name] += integral
                        max_cfl[name] = max(max_cfl[name], speed * p["dt"] / f.dx)
                        if max_cfl[name] > .5:
                            np.save(ROOT / f"failed-cfl-{name}.npy", fields[name])
                            raise FloatingPointError(f"{name} speed CFL exceeds 0.5")
                    delta = fields["disturbed"] - fields["reference"]
                    after_rate = f.inner(delta, cached["disturbed"][0] - cached["reference"][0])
                    integrated_difference_net += .5 * p["dt"] * (before_rate + after_rate)
                    if sub % 10 == 0:
                        save(ROOT / "progress.json", {"pid": os.getpid(), "status": "running",
                             "t": (index-1)*.01 + (sub+1)*p["dt"], "end": .4,
                             "step": (index-1)*p["steps_per_output"] + sub+1,
                             "total_steps": 40*p["steps_per_output"], "updated_epoch": time.time()})
            observations = {name: f.observe(h, index*.01, accum[name], initial[name])
                            for name, h in fields.items()}
            failures = []
            for name, row in observations.items():
                row["resolved_screen"] = screen(row)
                failures += [f"{name}: {message}" for message in check(row, f, p)]
            row = observations["disturbed"]
            row["reference_diagnostics"] = observations["reference"]
            row.update(compare(f, fields["disturbed"], fields["reference"],
                               p["reference_initial_l2"], p["amplitude"]))
            row["integrated_difference_net"] = integrated_difference_net
            residual = row["difference_energy"] - initial_difference_energy - integrated_difference_net
            row["difference_energy_budget_relative_residual"] = residual / initial_difference_energy
            if abs(row["difference_energy_budget_relative_residual"]) > .005:
                failures.append("difference-energy budget residual exceeds 0.5% of initial difference energy")
            if row["difference_mean_max"] > 1e-10:
                failures.append("difference mean drift exceeds 1e-10")
            paths, hashes = {}, {}
            for name, h in fields.items():
                paths[name] = f"fields/{name}/field-{index:03d}.npy"
                save_field(ROOT / paths[name], h)
                hashes[name] = sha(ROOT / paths[name])
            row["field_hashes"] = hashes
            if failures:
                row["hard_failure"] = failures
            rows.append(row)
            status = "failed_check" if failures else ("complete" if index == 40 else "running")
            result = {"status": status, "pid": os.getpid(), "job": p["job"],
                      "protocol_sha256": sha(ROOT / "protocol.json"), "last_output": index,
                      "fields": paths, "field_hashes": hashes, "series": rows,
                      "dt": p["dt"], "accum": {name: x.tolist() for name, x in accum.items()},
                      "max_speed_cfl": max_cfl, "integrated_difference_net": integrated_difference_net,
                      "updated_epoch": time.time(), "elapsed_this_process_seconds": time.perf_counter()-begin}
            save(result_path, result)
            from build_report import build
            build()
            print(f"t={index*.01:.2f} amplification={row['amplification']:.9g} "
                  f"net_rate={row['difference_energy_net_rate']:.9g} status={status}", flush=True)
            if failures:
                raise RuntimeError("; ".join(failures))
        verify_sources(p)
        save(ROOT / "progress.json", {"pid": os.getpid(), "status": "complete", "t": .4,
                                       "updated_epoch": time.time()})
    except BaseException:
        save(ROOT / "failure.json", {"pid": os.getpid(), "error": traceback.format_exc(),
                                      "updated_epoch": time.time()})
        raise
    finally:
        if lock.exists() and lock.read_text() == str(os.getpid()):
            lock.unlink()


if __name__ == "__main__":
    run()
