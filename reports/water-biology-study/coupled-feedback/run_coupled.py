"""Execute the frozen matrix; preserve completed cases and reject source drift."""
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
from coupled_model import run_pair

ROOT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def case_id(n,dt,seed,organization,feedback):return f'n{n}-dt{dt:g}-s{seed}-{organization}-'+('on' if feedback else 'off')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'output');args=ap.parse_args();out=args.out
    out.mkdir(parents=True,exist_ok=True);(out/'trajectories').mkdir(exist_ok=True)
    protocol=json.loads((ROOT/'protocol.json').read_text(encoding='utf-8'))
    sources={p:digest(ROOT/p) for p in ['coupled_model.py','protocol.json','run_coupled.py']}
    provenance={'source_sha256':sources,'python':platform.python_version(),'numpy':np.__version__}
    if (out/'provenance.json').exists():
        old=json.loads((out/'provenance.json').read_text(encoding='utf-8'))
        if old!=provenance:raise RuntimeError('Existing output has different source or environment; use a separate output folder')
    else:write(out/'provenance.json',provenance)
    cases=[(48,.004,s,o,f) for s in protocol['primary_seeds'] for o in ['x','y'] for f in [True,False]]
    cases += [(n,dt,100,o,f) for n,dt in protocol['additional_refinements'] for o in ['x','y'] for f in [True,False]]
    cases += [(48,.004,100,'uniform',f) for f in [True,False]]
    assert len(cases)==protocol['paired_jobs']
    results=[];final={}
    for case in cases:
        key=case_id(*case);jp=out/(key+'.json');fp=out/'trajectories'/(key+'.npz')
        if jp.exists() and fp.exists():
            r=json.loads(jp.read_text(encoding='utf-8'));h=np.load(fp)['h']
            assert r['trajectory_sha256']==digest(fp)
        else:
            started=time.perf_counter()
            try:r,h=run_pair(*case[:2],duration=protocol['duration'],seed=case[2],organization=case[3],feedback=case[4])
            except Exception as exc:
                write(out/(key+'-failure.json'),{'case':key,'error':repr(exc)});raise
            np.savez_compressed(fp,h=h);r['case']=key;r['trajectory_sha256']=digest(fp);write(jp,r)
            print(key,'complete',round(time.perf_counter()-started,2),'s',flush=True)
        results.append(r);final[key]=h[-1]
    write(out/'results.json',results);np.savez_compressed(out/'final-fields.npz',**final)
    print('COMPLETE',len(results),'paired jobs',flush=True)

if __name__=='__main__':main()
