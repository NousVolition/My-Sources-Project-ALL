"""Restore byte-identical final arrays from API-sized publication parts; no simulation or fitting."""
from pathlib import Path
import argparse, hashlib, json, os, shutil
HERE = Path(__file__).resolve().parent

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def restore(directory=HERE):
    directory = Path(directory)
    spec = json.loads((directory / "array-parts" / "delivery.json").read_text(encoding="utf-8"))
    primary = directory / spec["outputs"][0]
    for name in spec["outputs"]:
        existing = directory / name
        if existing.exists() and digest(existing) != spec["original_sha256"]:
            raise ValueError("Refusing to overwrite different existing array file: " + str(existing))
    parts = []
    for record in spec["parts"]:
        path = directory / "array-parts" / record["name"]
        if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
            raise ValueError("Part checksum mismatch: " + str(path))
        parts.append(path)
    if not primary.exists():
        temporary = primary.with_name(primary.name + ".restore-tmp")
        if temporary.exists():
            raise FileExistsError("Existing restore temporary file: " + str(temporary))
        try:
            with temporary.open("xb") as output:
                for path in parts:
                    with path.open("rb") as stream:
                        shutil.copyfileobj(stream, output, 1024 * 1024)
            if temporary.stat().st_size != spec["original_bytes"] or digest(temporary) != spec["original_sha256"]:
                raise ValueError("Reassembled array checksum mismatch")
            os.replace(temporary, primary)
        finally:
            if temporary.exists():
                temporary.unlink()
    for name in spec["outputs"][1:]:
        destination = directory / name
        if not destination.exists():
            shutil.copyfile(primary, destination)
        if digest(destination) != spec["original_sha256"]:
            raise ValueError("Output checksum mismatch: " + str(destination))
    return {"sha256": spec["original_sha256"], "outputs": spec["outputs"], "bytes": spec["original_bytes"]}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=HERE)
    print(json.dumps(restore(parser.parse_args().directory)))
