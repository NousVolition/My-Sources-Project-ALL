"""Safety checks against accepting under-resolution or incomplete evidence."""
from pathlib import Path
from unittest.mock import patch
import json,copy,hashlib
import adaptive_mesh as gate
ROOT=Path(__file__).resolve().parent

def main():
    paths=[ROOT/'parents'/f'mesh{n}-one.json' for n in [24,48,72]]
    checks=[]
    def check(name,condition):
        assert condition,name;checks.append(dict(name=name,passed=True))
    check('coarse-pair-rejected',not gate.select(paths[:2])['accepted'])
    selected=gate.select(paths)
    check('coarse-start-selects-actual-finer-mesh',selected['selected_mesh']==[72,7] and selected['selected']=='mesh72-one.json')
    original=Path.read_text
    def altered(change):
        def read(path,*args,**kwargs):
            text=original(path,*args,**kwargs)
            if path==paths[-1]:
                record=json.loads(text);change(record);return json.dumps(record)
            return text
        return read
    with patch.object(Path,'read_text',altered(lambda r:r.update(failure=dict(message='Test budget stop')))):
        check('incomplete-finer-run-rejected',not gate.select(paths)['accepted'])
    with patch.object(Path,'read_text',altered(lambda r:r.update(energy_balance_scaled=.1))):
        check('energy-failure-rejected',not gate.select(paths)['accepted'])
    with patch.object(Path,'read_text',altered(lambda r:r['audit'].update(minimum_swept_clearance=-1))):
        check('geometry-failure-rejected',not gate.select(paths)['accepted'])
    for label,change in [('different-interval',lambda r:r.update(requested_end=.1)),('different-initial-state',lambda r:r.update(initial=[0,0,0])),('different-physics',lambda r:r['parameters'].update(alpha=1))]:
        rejected=False
        with patch.object(Path,'read_text',altered(change)):
            try:gate.select(paths)
            except ValueError:rejected=True
        check(label+'-rejected',rejected)
    check('identical-mesh-is-not-refinement',not gate.compare(paths[-1],paths[-1])['passed'])
    result=dict(all_passed=True,checks=checks,test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),driver_sha256=gate.digest(ROOT/'adaptive_mesh.py'))
    (ROOT/'data/adaptive-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
