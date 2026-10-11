"""Dense/sparse agreement control and tighter historical thin-cell control."""
import json,hashlib
from pathlib import Path
from run_experiments import run,Parameters,DATA,hashes
def main():
    configs=[dict(name='sparse-agreement',p=Parameters(preload=.35,alpha=16,rho=160,viscosity=.15,feedback=False),end=.2,group='storage-control'),dict(name='old-method-fine-one',p=Parameters(preload=.35,alpha=16,rho=160,viscosity=.15,feedback=False),old=True,rtol=2e-9,atol=2e-11,group='method-control')]
    declared=[dict(**{k:v for k,v in c.items() if k!='p'},parameters=c['p'].__dict__) for c in configs]
    (DATA/'solver-upgrade-protocol.json').write_text(json.dumps(dict(configurations=declared,sources=hashes(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    for c in configs:run(**c)
if __name__=='__main__':main()
