"""Wait for production, overlap CPU analysis with sequential GPU controls."""
import json,os,subprocess,sys,time
from pathlib import Path

if __name__=='__main__':
    root=Path(__file__).resolve().parent;deadline=time.monotonic()+1200
    while True:
        paths=[root/'data'/f'production-execution-{s}.json' for s in ['primary','refinement','sensitivity']]
        if all(p.exists() and (lambda d:d['completed']==d['planned'])(json.loads(p.read_text())) for p in paths):break
        if time.monotonic()>deadline:raise RuntimeError('Production did not complete; stopping follow-ups')
        time.sleep(2)
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1'
    analysis=subprocess.Popen([sys.executable,'-B','analyze.py'],cwd=root,env=env)
    def run(name):
        print('Running '+name,flush=True)
        subprocess.run([sys.executable,'-B',name],cwd=root,env=env,check=True)
    run('matched_sham.py');run('delayed_probe.py')
    if analysis.wait()!=0:raise RuntimeError('Core analysis failed')
    run('followup_analysis.py');run('template_analysis.py');run('impulse_energy.py');run('verify.py');run('report.py')
