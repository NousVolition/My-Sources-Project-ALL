"""Create explicit acceptance records from matching, completed recordings."""
from pathlib import Path
import json,hashlib
from accuracy_gate import assess
ROOT=Path(__file__).resolve().parent
def main():
    cases=[('gate-final-coarse-preview','mesh24-one','mesh48-one','mesh'),('gate-final-one','mesh48-one','mesh72-one','mesh'),('gate-final-two','mesh48-two','mesh72-two','mesh'),('gate-final-hard-one','hard48-one','hard72-one','mesh'),('gate-final-hard-two','hard48-two','hard72-two','mesh'),('gate-time-one','stress16-one','accurate16-one','time'),('gate-time-two','stress16-two','accurate16-two','time')]
    records=[]
    for name,a,b,kind in cases:
        ca=ROOT/'solver/data'/(a+'.json');fi=ROOT/'solver/data'/(b+'.json')
        if not ca.exists() or not fi.exists():raise ValueError('Every declared gate needs a recorded attempt: '+name)
        g=assess(ca,fi,kind);g['coarse']=ca.relative_to(ROOT).as_posix();g['fine']=fi.relative_to(ROOT).as_posix()
        g.update(coarse_sha256=hashlib.sha256(ca.read_bytes()).hexdigest(),fine_sha256=hashlib.sha256(fi.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        (ROOT/'data'/(name+'.json')).write_text(json.dumps(g,indent=2)+'\n');records.append(dict(name=name,**g))
    result=dict(checks=records,coarse_preview_rejected=not records[0]['passed'],standard_mesh_passed=all(g['passed'] for g in records[1:3]),hard_mesh_passed=all(g['passed'] for g in records[3:5]),time_checks_passed=all(g['passed'] for g in records[5:]))
    (ROOT/'data/final-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
