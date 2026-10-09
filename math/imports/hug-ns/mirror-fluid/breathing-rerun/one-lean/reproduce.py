from pathlib import Path
import json,shutil,subprocess,sys
p=Path(__file__).resolve().parent;r=p/'rerun';r.mkdir(exist_ok=True)
for n in ('one_lean.py','clay_hug.py'):shutil.copyfile(p/n,r/n)
subprocess.run([sys.executable,str(r/'one_lean.py')],check=True)
a=json.loads((p/'reference-one-table.json').read_text());b=json.loads((r/'one-table.json').read_text())
assert set(a)==set(b)
for k in a:
    assert len(a[k])==len(b[k])
    for x,y in zip(a[k],b[k]):
        assert len(x)==len(y)
        assert all(abs(v-w)<1e-12 for v,w in zip(x,y)),k
print('All four series reproduced.')
