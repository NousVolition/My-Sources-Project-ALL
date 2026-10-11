"""Own the bounded matrix; conditional extensions follow frozen numerical screens."""
from pathlib import Path
import os
import sys
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
from run_history import sources, validate, name_for, sha, save_json

HERE=Path(__file__).resolve().parent


def subprocess_job(seed,config,frozen):
    if sources()!=frozen: raise ValueError('Frozen sources changed')
    name=name_for(seed,config)
    dest=HERE/'recorded-data'/(name+'.npz')
    if config['n']==56 and config['dt']==.0025 and seed in [1024,1025]:
        dest=HERE.parent/'resolution-followup'/'recorded-data'/dest.name
        validate(dest,config)
    elif dest.exists():
        validate(dest,config)
    else:
        with (HERE/(name+'.log')).open('w') as log:
            child=subprocess.run([sys.executable,'-B','-u',str(HERE/'run_history.py'),
                '--seed',str(seed),'--grid',str(config['n']),'--dt',str(config['dt'])],
                stdout=log,stderr=subprocess.STDOUT,check=False,creationflags=subprocess.CREATE_NO_WINDOW)
        if child.returncode:raise RuntimeError(name+f' failed ({child.returncode}); preserved log/checkpoints')
        validate(dest,config)
    m=json.loads(dest.with_suffix('.json').read_text())
    return {'path':os.path.relpath(dest,HERE).replace('\\','/'),'name':name,'seed':seed,'config':config,
            'sha256':m['sha256'],'seconds':m['seconds'],'verified':True,'reused':dest.parent.parent!=HERE}


def step(script,log):
    with (HERE/log).open('w') as stream:
        c=subprocess.run([sys.executable,'-B',str(HERE/script)],stdout=stream,stderr=subprocess.STDOUT,check=False)
    if c.returncode:raise RuntimeError(script+f' failed ({c.returncode}); preserved '+log)


def main():
    p=json.loads((HERE/'protocol.json').read_text()); frozen=sources()
    lock=HERE/'controller.json'
    with lock.open('x') as stream:json.dump({'pid':os.getpid(),'sources':frozen,'status':'running'},stream)
    try:
        # Pilot is launched and verified separately before this controller starts.
        pilot=HERE/'recorded-data'/(name_for(p['pilot']['seed'],{**p['base'],'n':p['pilot']['grid']})+'.npz')
        pm=validate(pilot,{**p['base'],'n':p['pilot']['grid']})
        save_json(HERE/'pilot.json',{'verified':True,'sources':frozen,'seed':pm['seed'],'config':pm['config'],
            'sha256':pm['sha256'],'seconds':pm['seconds'],'max_cfl':pm['max_cfl']})
        records=[]; grids=[]; stages=[p['initial_grids']]+[[n] for n in p['conditional_extension_grids']]
        for stage in stages:
            if sources()!=frozen:raise ValueError('Source change before next stage')
            jobs=[(seed,{**p['base'],'n':n},'main') for n in stage for seed in p['test_seeds']]
            jobs += [(seed,{**p['base'],'n':n,'dt':p['half_step_dt']},'half')
                for n in stage if n in p['half_step_grids'] for seed in p['half_step_seeds']]
            failures=[]; new=[]
            with ThreadPoolExecutor(max_workers=p['workers']) as pool:
                pending={pool.submit(subprocess_job,seed,config,frozen):(seed,config,kind) for seed,config,kind in jobs}
                for future in as_completed(pending):
                    seed,config,kind=pending[future]
                    try:
                        record=future.result();record['kind']=kind;new.append(record)
                        print(record['name']+' verified',flush=True)
                    except Exception as exc:
                        failures.append({'seed':seed,'config':config,'kind':kind,'error':repr(exc)})
            if failures:
                save_json(HERE/'stage-failures.json',{'failures':failures,'completed':new,'sources':frozen})
                raise RuntimeError('Stage failures retained after unrelated jobs finished')
            records+=new;records.sort(key=lambda r:(r['kind'],r['config']['n'],r['seed']))
            grids+=stage;finest=grids[-1]
            execution={'complete':True,'records':records,'grids':grids,'sources':frozen}
            save_json(HERE/'execution.json',execution);save_json(HERE/f'execution-{finest}.json',execution)
            step('analyze_history.py',f'analysis-{finest}.log')
            step('build_report.py',f'report-{finest}.log')
            result=json.loads((HERE/'results.json').read_text())
            for filename in ['results.json','predictions-and-features.npz','convergence.png','convergence.svg','per-run-comparisons.csv','README.md']:
                x=HERE/filename;shutil.copyfile(x,x.with_name(x.stem+f'-{finest}'+x.suffix))
            decision={'grid':finest,'screens_pass':result['all_screens_pass'],'failures':result['failures'],
                      'protocol_sha256':sha(HERE/'protocol.json'),'results_sha256':sha(HERE/'results.json')}
            save_json(HERE/f'decision-{finest}.json',decision)
            print(f'Stage{finest}: {sum(s["pass"] for s in result["screens"])}/{len(result["screens"])} screens pass',flush=True)
            if result['all_screens_pass']:break
        save_json(HERE/'completion.json',{'complete':True,'sources':frozen,'grids':grids,'records':len(records),
            'new_simulations':sum(not r['reused'] for r in records),'all_screens_pass':result['all_screens_pass'],
            'results_sha256':sha(HERE/'results.json'),'report_sha256':sha(HERE/'README.md'),
            'analysis_sha256':sha(HERE/'analyze_history.py'),'frozen_models_sha256':p['models']['sha256']})
        print('Bounded full-history study complete; publication pending',flush=True)
    except Exception as exc:
        save_json(HERE/'controller-failure.json',{'error':repr(exc),'sources':frozen})
        raise
    finally:lock.unlink()


if __name__=='__main__':main()
