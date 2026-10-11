"""Additional barrier-model sensitivity, declared separately from main protocol."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict,replace
import hashlib,json
from contact_hug import Parameters
from run_experiments import run,write,ROOT,DATA
def main():
    base=Parameters(preload=.35,alpha=0);configs=[]
    for k in [1000,10000,100000]:
        configs.append(dict(name='barrier-strength'+str(k),parameters=asdict(replace(base,contact=k,wall=k)),initial=[-8,8,27],end=4,dt=.004,group='barrier-model'))
    for width in [.02,.04,.08]:
        configs.append(dict(name='barrier-width'+str(width),parameters=asdict(replace(base,barrier_width=width)),initial=[-8,8,27],end=4,dt=.004,group='barrier-model'))
    write('sensitivity-protocol.json',dict(configurations=configs,sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['contact_hug.py','run_experiments.py','run_sensitivity.py']}))
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run,c) for c in configs if not (DATA/(c['name']+'.json')).exists()]
        for f in as_completed(futures):
            r=f.result();print(r['name'],r['failure'],'clearance',r['audit']['minimum_swept_clearance'],flush=True)
    write('sensitivity.json',dict(runs=[json.loads((DATA/(c['name']+'.json')).read_text()) for c in configs]))
if __name__=='__main__':main()
