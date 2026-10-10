"""Repeat all sequence cases and compare every retained numerical array."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from sequence_transport import run

ROOT = Path(__file__).resolve().parent
DATA = ROOT/"data"/"sequence"


def fingerprints():
    result = {}
    for path in sorted(DATA.glob("*.npz")):
        with np.load(path) as arrays:
            result[path.name] = {name: dict(shape=list(arrays[name].shape), dtype=str(arrays[name].dtype),
                sha256=hashlib.sha256(arrays[name].tobytes()).hexdigest()) for name in sorted(arrays.files)}
    return result


def main():
    before = fingerprints()
    if not before:
        raise RuntimeError("Run sequence_transport.py before requesting reproduction")
    previous_source = json.loads((DATA/"summary.json").read_text())["source_sha256"]
    started = time.perf_counter()
    run()
    after = fingerprints()
    current_source = json.loads((DATA/"summary.json").read_text())["source_sha256"]
    result = dict(numerical_files=len(before), arrays_compared=sum(len(v) for v in before.values()),
        every_retained_array_identical=before == after, elapsed_seconds=time.perf_counter()-started,
        previous_run_source_sha256=previous_source, repeated_run_source_sha256=current_source,
        scope="All 42 retained numerical cases; every accepted time and diagnostic, plus each cycle-boundary spatial state. Unsaved intermediate spatial arrays are not part of this full-package comparison; two selected cases separately repeat every state.")
    (DATA/"reproduction.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result))
    if before != after:
        raise AssertionError("Numerical reproduction differs; inspect the retained record")


if __name__ == "__main__":
    main()
