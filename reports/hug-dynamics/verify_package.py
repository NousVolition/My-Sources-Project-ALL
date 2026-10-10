"""Check SHA-256 hashes of the published source, reports and numerical data."""
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parent
manifest=json.loads((root/'checksums.json').read_text(encoding='utf-8'))
failed=[]
for name,expected in manifest.items():
    path=root/name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
        failed.append(name)
print(json.dumps(dict(files=len(manifest),passed=len(manifest)-len(failed),failed=failed),indent=2))
raise SystemExit(bool(failed))
