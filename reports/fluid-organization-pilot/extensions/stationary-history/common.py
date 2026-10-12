import hashlib,json,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,data):
 p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(path):return json.loads((ROOT/path).read_text(encoding='utf-8'))
def protocol():return read('protocol.json')
def sources():return {name:digest(ROOT/name) for name in ['solver.py','kernels.py','simulate.py','common.py']}

