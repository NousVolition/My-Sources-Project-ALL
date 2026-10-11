"""Preserve failed implementations; create a version with exact integer waves."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
def copy(a,b,replacements):
    if (root/b).exists():raise RuntimeError('Refuse to replace corrected source')
    text=(root/a).read_text()
    for before,after in replacements:
        if before not in text:raise RuntimeError('Missing replacement '+before)
        text=text.replace(before,after)
    (root/b).write_text(text,encoding='utf8',newline='\n')
copy('reference_solver.py','corrected_solver.py',[
    ('m = fft.fftfreq(n)*n','m = np.rint(fft.fftfreq(n)*n).astype(np.int64)'),
    ('z = fft.rfftfreq(n)*n','z = np.rint(fft.rfftfreq(n)*n).astype(np.int64)'),
    ('np.max(abs(k), axis=0) < n/3','3*np.max(abs(k), axis=0) < n'),
    ('out=np.zeros_like(k[0])','out=np.zeros_like(k[0],dtype=float)')])
copy('fast_observation.py','corrected_observation.py',[('from reference_solver import Solver','from corrected_solver import Solver')])
copy('production_experiment.py','final_experiment.py',[
    ('from reference_solver import random_field, interventions','from corrected_solver import random_field, interventions'),
    ('from fast_observation import FastSolver as Solver','from corrected_observation import FastSolver as Solver'),
    ("['protocol.json','reference_solver.py','fast_observation.py','production_experiment.py']","['protocol.json','corrected_solver.py','corrected_observation.py','final_experiment.py']"),
    ("ROOT/'production-spinup'","ROOT/'final-spinup'"),("ROOT/'production-data'","ROOT/'final-data'"),
    ("if frozen.with_suffix('.json').exists():","if False and frozen.with_suffix('.json').exists():"),
    ("if original.with_suffix('.json').exists():","if False and original.with_suffix('.json').exists():"),
    ('batch=16 if n>=64 else len(labels)','batch=4 if n>=128 else (16 if n>=64 else len(labels))')])
for name in ['common_probe.py','unforced_readout.py','common_unforced.py','independent_check_full.py']:
    p=root/name;text=p.read_text().replace('from production_experiment import','from final_experiment import')
    text=text.replace("'production_experiment.py'","'final_experiment.py'").replace("'reference_solver.py'","'corrected_solver.py'").replace("'fast_observation.py'","'corrected_observation.py'")
    p.write_text(text,encoding='utf8',newline='\n')
p=root/'analysis.py';p.write_text(p.read_text().replace("folder='production-data'","folder='final-data'"),encoding='utf8',newline='\n')
note=dict(date='2026-10-11',production_effects_inspected=False,failed_attempts_preserved=True,
    correction='Round FFT frequencies to exact integers before all spectral operations; use the exact strict cutoff 3*max(abs(k))<n. The reset reference has exactly parallel integer vectors on some modes; tiny floating-wave errors incorrectly bypassed its alternate-vector fallback and were amplified by normalization.',
    original_sources={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ['reference_solver.py','production_experiment.py']},
    corrected_sources={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ['corrected_solver.py','corrected_observation.py','final_experiment.py']},
    spinup='Repeat with corrected solver; do not reuse original spinup as an interchangeable result.',
    validation='Grid-independent reset regression at n=64/80/96 and repeated numerical pilot matrix. No thresholds relaxed.')
(root/'integer-wave-correction.json').write_text(json.dumps(note,indent=2)+'\n',encoding='utf8')
