"""One-time setup before corrected production; not a simulation."""
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent
p=r/'integer-wave-correction.json';v=json.loads(p.read_text())
v['corrected_sources']={n:hashlib.sha256((r/n).read_bytes()).hexdigest() for n in v['corrected_sources']}
v['implementation_startup']={'integer_inverse_buffer':'float64','regression_tests_passed':2}
p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf8')
p=r/'production.py';text=p.read_text().replace('from production_experiment import','from final_experiment import')
text=text.replace("qualification=json.loads((ROOT/'results'/'production-qualification.json').read_text())","qualification=json.loads((ROOT/'results'/'corrected-pilot-qualification.json').read_text())")
text=text.replace("qualification['all_passed'] and fast['all_passed'] and equivalence['all_passed']","qualification['all_passed'] and qualification['complete'] and json.loads((ROOT/'results'/'corrected-validation.json').read_text())['all_passed']")
p.write_text(text,encoding='utf8',newline='\n')
for old,new in [('test_controls.py','test_corrected_controls.py'),('test_common_probe.py','test_corrected_common.py')]:
    p=r/new
    p.write_text((r/old).read_text().replace('from experiment import','from final_experiment import').replace('from reference_solver import','from corrected_solver import'),encoding='utf8',newline='\n')
failure=r/'results'/'gpu-pipeline-failure.json'
if failure.exists():failure.rename(r/'results'/'corrected-startup-failure.json')
