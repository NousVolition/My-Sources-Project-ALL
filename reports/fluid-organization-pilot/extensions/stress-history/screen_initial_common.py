from pathlib import Path
import json,hashlib,csv
import numpy as np
ROOT=Path(__file__).resolve().parent
def save(name,obj):
    p=ROOT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(name,rows):
    p=ROOT/'results'/name;p.parent.mkdir(exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def read_protocol(screen=False):return json.loads((ROOT/('screen_protocol.json' if screen else 'protocol.json')).read_text(encoding='utf-8'))
