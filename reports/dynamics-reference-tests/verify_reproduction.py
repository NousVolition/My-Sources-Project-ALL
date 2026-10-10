"""Repeat the study and compare numerical outputs with the preceding run.

Exact repeatability is checked on this machine/environment. It is not promised
across platforms for long chaotic trajectories. NPZ timestamps are ignored;
decompressed array values, names, shapes and dtypes are hashed instead.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parent


def snapshot():
    hashes={}
    for path in sorted((ROOT/"data").glob("*")):
        if path.suffix==".csv":
            hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        elif path.suffix==".npz":
            digest=hashlib.sha256()
            with np.load(path,allow_pickle=False) as arrays:
                for key in sorted(arrays.files):
                    a=arrays[key]
                    digest.update(key.encode())
                    digest.update(str(a.shape).encode())
                    digest.update(a.dtype.str.encode())
                    digest.update(a.tobytes())
            hashes[path.name]=digest.hexdigest()
        elif path.name in ("summary.json","checks.json"):
            hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def main():
    before=snapshot()
    if not before: raise RuntimeError("Run run_tests.py first")
    process=subprocess.run([sys.executable,str(ROOT/"run_tests.py")],capture_output=True,text=True)
    (ROOT/"data"/"repeat-run.log").write_text(process.stdout+process.stderr,encoding="utf-8")
    after=snapshot()
    changed=[name for name in sorted(set(before)|set(after)) if before.get(name)!=after.get(name)]
    result=dict(repeated_run_exit_code=process.returncode,numerical_files_compared=len(before),
                identical=not changed and process.returncode==0,changed_files=changed,
                previous_hashes=before,repeated_hashes=after,
                scope="Same machine, package versions, initial conditions and parameters. Image metadata and environment timing excluded.")
    (ROOT/"data"/"reproduction.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("hashes")},indent=2))
    if not result["identical"]: sys.exit(1)


if __name__=="__main__": main()
