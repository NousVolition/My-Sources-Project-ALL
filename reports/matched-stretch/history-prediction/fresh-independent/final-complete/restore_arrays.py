"""Restore lossless prediction archive parts without overwriting different bytes."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent
def sha(b):return hashlib.sha256(b).hexdigest()
def restore():
    manifest=json.loads((HERE/'prediction-parts.json').read_text())
    for row in manifest['archives']:
        path=HERE.parent/row['name']
        if path.exists():
            assert sha(path.read_bytes())==row['sha256'],str(path)
            continue
        pieces=[]
        for part in row['parts']:
            data=(HERE/part['path']).read_bytes()
            assert len(data)==part['size'] and sha(data)==part['sha256'],part['path']
            pieces.append(data)
        data=b''.join(pieces)
        assert len(data)==row['size'] and sha(data)==row['sha256'],row['name']
        with path.open('xb') as stream:stream.write(data)
    print(str(len(manifest['archives']))+' prediction archives restored or verified byte for byte')
if __name__=='__main__':restore()
