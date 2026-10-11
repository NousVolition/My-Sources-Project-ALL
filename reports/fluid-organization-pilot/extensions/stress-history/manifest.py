import json
from common import ROOT,digest,save
files=[]
for q in sorted(ROOT.rglob('*')):
    rel=q.relative_to(ROOT)
    if not q.is_file() or '__pycache__' in rel.parts or '.pytest_cache' in rel.parts or rel.as_posix()=='manifest.json':continue
    if q.suffix not in ('.py','.txt','.md','.json','.npz','.csv','.html','.png','.xml'):continue
    files.append({'path':rel.as_posix(),'size':q.stat().st_size,'sha256':digest(q)})
save('manifest.json',{'note':'Exact executed source and completed numerical outputs; excludes this inventory, cache files and later publication archive wrappers. Full fields regenerated from seeds, observables saved losslessly.','files':files})
print('Manifest files:',len(files))
