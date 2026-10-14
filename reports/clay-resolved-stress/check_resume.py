"""Verify exact checkpoint preservation and an independently computed tail."""
from pathlib import Path
import json,numpy as np
from analyze import load,compare,sha
ROOT=Path(__file__).resolve().parent
def main():
 checks=[]
 for p in sorted((ROOT/'solver-versions/shared-resume/data').glob('*.json')):
  r=json.loads(p.read_text())
  if 'parameters' not in r:continue
  child=load(p,join=False);parent=load(ROOT/r['parent']['metadata']);assert np.array_equal(child[3][0],parent[3][-1]) and child[2][0]==parent[2][-1] and r['audit']['checkpoint_clock_preserved']
  checks.append(dict(name=r['name'],all_state_and_work_values_preserved=True,clock_preserved=True,passed=True))
 ref=load(ROOT/'solver-versions/shared/data/shared144-one.json');tail=load(ROOT/'solver-versions/shared-resume/data/resume-check-144-one.json');c=compare(ref,tail,start=tail[0]['start_time']);passed=c['state_scaled']<1e-4 and c['lorenz_scaled']<1e-10 and tail[0]['failure'] is None and c['end']==.5;assert passed
 result=dict(all_passed=True,preservation_checks=checks,independent_tail_check=c,tail_check_passed=bool(passed),source_sha256=sha(__file__));(ROOT/'data/resume-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,preserved=len(checks),tail_state_scaled=c['state_scaled'],start=c['start'],end=c['end'])))
if __name__=='__main__':main()
