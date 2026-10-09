"""Verify archived constructions and tracks without redoing any integration."""
from pathlib import Path
import json,hashlib
import numpy as np
from menger_tracks import construction,initial_positions
from quarter_connected import build
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    original=json.loads((ROOT/'other-results.json').read_text()) if (ROOT/'other-results.json').exists() else None
    m=json.loads((ROOT/'menger-results.json').read_text());mc=np.load(ROOT/'menger-construction.npz')
    for d,p in construction().items():assert np.array_equal(p,mc['menger'+str(d)])
    q,flips=build();qc=np.load(ROOT/'quarter-construction.npz')
    for d,p in q.items():assert np.array_equal(p,qc['level'+str(d)])
    fields={};files={};extraction=[]
    for result,geom,prefix,script in [('results.json','construction.npz','','run_fractal.py'),('quarter-results.json','quarter-construction.npz','quarter-','quarter_connected.py'),('menger-results.json','menger-construction.npz','menger-',None)]:
        data=json.loads((ROOT/result).read_text());assert sha(ROOT/geom)==data.get('geometry_sha256',data.get('construction_sha256'))
        if script:assert sha(ROOT/script)==data['run_script_sha256']
        for name,digest in data['fields'].items():
            if name in fields:assert fields[name]==digest
            fields[name]=digest
        for name,r in data['runs'].items():
            file=ROOT/(prefix+name+'.npz');assert sha(file)==r['trajectory_sha256'];z=np.load(file,allow_pickle=False)
            assert np.isfinite(z['positions']).all() and z['time'][0]==0 and z['time'][-1]==.1
            assert np.all(np.diff(z['time'])>0);files[file.name]=sha(file)
            if prefix=='menger-':
                assert np.array_equal(z['positions'][0],initial_positions())
                if original:
                    fullpath=ROOT/('other-'+name+'.npz');assert sha(fullpath)==r['original_full_trajectory_sha256']
                    full=np.load(fullpath);start=next(g['start'] for g in original['groups'] if g['name']=='menger2')
                    assert np.array_equal(z['positions'],full['positions'][:,start:]) and np.array_equal(z['time'],full['time'])
                    extraction.append(name)
    a=json.loads((ROOT/'analysis.json').read_text());assert a['verification']['raw_files']==files
    # Local source availability permits a fresh check; portable packages retain hashes.
    source=ROOT.parents[1]/'runs';available=source.exists()
    if available:
        for name,digest in fields.items():assert sha(source/name)==digest
        for name,digest in m['source_hashes'].items():assert sha(ROOT.parents[1]/name)==digest
    out={'status':'passed','trajectory_files':files,'source_fields':fields,'source_hashes':m['source_hashes'],'constructions_regenerated_exactly':True,'correct_quarter_flips':flips,'menger_extraction_verified_locally':extraction,'source_fields_rechecked_locally':available,'menger_reproducer_sha256':sha(ROOT/'menger_tracks.py'),'menger_reproducer_note':'Geometry matches saved labels exactly. Historical independent-marker RK4 tracks are preserved; no redundant integration was run.','scope':'Finite archive and construction consistency; no physical fractal convergence claim.'}
    (ROOT/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','tracks':len(files),'fields':len(fields),'menger_lossless_subsets':len(extraction)}))
if __name__=='__main__':main()
