"""Run the frozen design in a fresh directory, CPU by default (GPU optional)."""
import argparse,os,shutil,subprocess,sys
from pathlib import Path


def main():
    source=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser();ap.add_argument('--destination',type=Path,required=True);ap.add_argument('--backend',choices=['cpu','gpu'],default='cpu');args=ap.parse_args()
    dest=args.destination.resolve()
    if dest.exists():raise RuntimeError('Choose a new nonexistent destination; existing work is never overwritten')
    dest.mkdir(parents=True)
    for p in source.iterdir():
        if p.is_file() and p.suffix in ['.py','.json','.txt','.md']:shutil.copy2(p,dest/p.name)
    env=os.environ.copy();env['FLUID_BACKEND']=args.backend;env['OPENBLAS_NUM_THREADS']='1'
    def run(*cmd):subprocess.run([sys.executable,'-B',*cmd],cwd=dest,env=env,check=True)
    run('-m','unittest','test_solver','test_analysis','test_matched_sham')
    import json
    p=json.loads((dest/'protocol.json').read_text())
    for s in p['pilot_seeds']:
        for n in p['grids']:run('run.py','--seed',str(s),'--n',str(n),'--backend',args.backend,'--out','pilot-data')
    run('analyze.py','--pilot','--data','pilot-data','--out','pilot-results')
    run('production.py','--stage','primary','--backend',args.backend)
    run('freeze_models.py')
    run('production.py','--stage','refinement','--backend',args.backend)
    run('production.py','--stage','sensitivity','--backend',args.backend)
    run('matched_sham.py','--backend',args.backend)
    run('delayed_probe.py','--backend',args.backend)
    run('analyze.py');run('followup_analysis.py');run('template_analysis.py');run('impulse_energy.py');run('verify.py');run('report.py')
    print('Completed fresh reproduction in '+str(dest))


if __name__=='__main__':main()
