import json,math,hashlib
from pathlib import Path
from audit import ROOT
d=json.loads((ROOT/'results.json').read_text())
def finite(x):
 if isinstance(x,float):assert math.isfinite(x)
 elif isinstance(x,dict):
  for v in x.values():finite(v)
 elif isinstance(x,list):
  for v in x:finite(v)
finite(d)
assert {r['n'] for r in d['initial']}=={64,80,112,128,256}
assert not any(r['gate']['pass'] for r in d['initial'])
assert d['new_time_evolution'] is False and d['common_resolved_interval'] is None
errors=[];closure=[]
for run in d['runs'].values():
 for r in run['rows']:
  assert r['profile']['width_cells']>0
  errors += [abs(sum(r['spectra']['energy_sum'])-r['energy'])/r['energy'],abs(sum(r['spectra']['enstrophy_sum'])-r['enstrophy'])/r['enstrophy']]
  closure.append(r['peak_budget']['closure_relative_l2'])
assert max(errors)<1e-11
assert max(closure)<1e-9,max(closure)
for name in ('audit.py','report.py','check_diagnostics.py','source/numerics.py'):compile((ROOT/name).read_bytes(),name,'exec')
out=dict(status='passed',saved_fields=sum(len(r['rows']) for r in d['runs'].values()),initial_grids=5,
  spectral_parseval_max_relative=max(errors),vorticity_budget_max_relative_closure=max(closure),all_measurements_finite=True,
  code_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('audit.py','report.py','check_diagnostics.py','source/numerics.py')},
  interpretation='The numerical diagnostics are verified; all initial width gates fail. This is not a resolution pass.')
(ROOT/'verification.json').write_bytes(json.dumps(out,indent=2).encode());print(json.dumps(out))
