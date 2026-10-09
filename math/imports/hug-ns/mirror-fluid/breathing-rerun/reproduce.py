from pathlib import Path
import json, shutil, subprocess, sys
root=Path(__file__).resolve().parent
run=root/'rerun';run.mkdir(exist_ok=True)
for name in ('breathing_rerun.py','clay_hug.py'):shutil.copyfile(root/name,run/name)
subprocess.run([sys.executable,str(run/'breathing_rerun.py')],check=True)
reference=json.loads((root/'reference-results.json').read_text())
actual=json.loads((run/'results.json').read_text())
assert set(reference)==set(actual)
largest=0.
for name in reference:
    assert len(reference[name])==len(actual[name])==2
    for a,b in zip(reference[name],actual[name]):
        for key in ('t','E','D'):
            error=abs(a[key]-b[key]);largest=max(largest,error)
            assert error<1e-12,(name,key,error)
print('All six cases reproduced. Largest absolute difference:',largest)
