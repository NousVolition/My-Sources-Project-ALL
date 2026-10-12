from common import ROOT,digest,save
files=[]
for p in sorted(ROOT.rglob('*')):
 rel=p.relative_to(ROOT)
 if not p.is_file() or any(x in rel.parts for x in ['__pycache__','.pytest_cache','scratch']) or rel.as_posix()=='manifest.json':continue
 if p.suffix not in ['.py','.txt','.md','.json','.npz','.csv','.html','.png','.xml','.log']:continue
 files.append(dict(path=rel.as_posix(),size=p.stat().st_size,sha256=digest(p)))
save('manifest.json',dict(note='Executed source, protocols, tests, completed measurements and retained screening outcomes. Excludes this manifest and publication wrappers. Full baseline fields are reproducible from seeds; selected restart fields are retained.',files=files))
print('Files inventoried:',len(files))
