"""Fresh reproduction that cannot mistake compact metadata for cached raw data.
Usage: python reproduce.py --out path-to-new-directory --backend cpu|gpu
Install requirements.txt (CPU) or requirements-gpu.txt (GPU) first.
"""
from pathlib import Path
import argparse, shutil, subprocess, sys

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--backend',choices=['cpu','gpu'],default='cpu');a=p.parse_args()
    src=Path(__file__).resolve().parent;dest=a.out.resolve()
    if dest.exists():raise SystemExit('Destination must be a new directory; existing studies are never overwritten.')
    dest.mkdir(parents=True)
    for f in src.glob('*.py'):shutil.copy2(f,dest/f.name)
    for name in ['protocol.json','gpu-amendment.json','analysis-amendment.json','refinement-plan.json','requirements.txt','requirements-gpu.txt']:shutil.copy2(src/name,dest/name)
    def run(*args):subprocess.run([sys.executable,*args],cwd=dest,check=True)
    run('test_solver.py')
    run('test_analysis.py')
    run('-c','from run import dump,ROOT,tasks;dump(ROOT/"planned_tasks.json",tasks())')
    if a.backend=='gpu':
        run('gpu_solver.py');run('run_gpu.py');run('refine.py')
    else:
        run('run.py','--workers','2')
        run('-c','import json;from run import ROOT,run;[print(run(t),flush=True) for t in json.loads((ROOT/"refinement-plan.json").read_text())["cases"]]')
    run('analyze.py');run('report.py')
    print('Fresh reproduction complete:',dest/'report.html')

if __name__=='__main__':main()
