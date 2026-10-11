"""Run dependent stages after numerical prerequisites; no external publication."""
import json,time
from experiment import ROOT,dump,digest
deadline=time.time()+7200
while True:
    ready=ROOT/'results'/'accelerated-pilot-equivalence.json'
    if ready.exists():break
    if time.time()>deadline:raise RuntimeError('Pilot completion timeout; inspect running job')
    time.sleep(10)
from qualify_production import qualify
qualify()
from production import execute,tasks
for stage in ['main','resolution_audit','translation','common_probe','common_resolution','unforced_readout','common_unforced_readout']:
    print('START '+stage,flush=True);execute(stage);print('FINISH '+stage,flush=True)
dump(ROOT/'results'/'gpu-pipeline-complete.json',dict(status='complete',controller_sha256=digest(ROOT/'gpu_pipeline.py'),counts={k:len(v) for k,v in tasks().items()}))
