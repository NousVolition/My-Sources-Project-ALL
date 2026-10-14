"""Checks for false convergence caused by incomplete or incompatible evidence."""
import copy,json,hashlib
from pathlib import Path
from analyze import load,compare
ROOT=Path(__file__).resolve().parent
def main():
 a=load(ROOT/'parents/long72-two.json');b=load(ROOT/'solver/data/fine96-two.json');checks=[]
 def check(name,condition):assert condition,name;checks.append(dict(name=name,passed=True))
 check('actual fine reciprocal pair',compare(a,b)['full_mesh_check_passed'])
 check('same mesh rejected',not compare(b,b)['full_mesh_check_passed'])
 def mutate(key,value):
  parts=list(b);parts[0]=copy.deepcopy(b[0]);parts[0][key]=value;return tuple(parts)
 check('failure rejected',not compare(a,mutate('failure',dict(message='budget')) )['full_mesh_check_passed'])
 check('different requested end rejected',not compare(a,mutate('requested_end',.6))['full_mesh_check_passed'])
 check('different initial state rejected',not compare(a,mutate('initial',[1,1,1]))['full_mesh_check_passed'])
 p=copy.deepcopy(b[0]['parameters']);p['rho']=501
 check('different physics rejected',not compare(a,mutate('parameters',p))['full_mesh_check_passed'])
 check('bad energy rejected',not compare(a,mutate('energy_balance_scaled',.1))['full_mesh_check_passed'])
 audit=copy.deepcopy(b[0]['audit']);audit['minimum_swept_quality_ratio']=.3
 check('bad geometry rejected',not compare(a,mutate('audit',audit))['full_mesh_check_passed'])
 # Isolate input mismatch from the otherwise passing gap and geometry checks.
 aa=list(a);bb=list(b)
 for parts in [aa,bb]:parts[0]=copy.deepcopy(parts[0]);parts[0]['parameters']['feedback']=False
 bb[3]=bb[3].copy();bb[3][:,-7:-4]+=100
 mismatch=compare(tuple(aa),tuple(bb))
 check('unmatched one-way input rejected despite passing gaps',mismatch['gap_target_met'] and not mismatch['drive_agreement_passed'] and not mismatch['full_mesh_check_passed'])
 result=dict(all_passed=True,checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),analysis_sha256=hashlib.sha256(ROOT.joinpath('analyze.py').read_bytes()).hexdigest())
 (ROOT/'data/acceptance-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
