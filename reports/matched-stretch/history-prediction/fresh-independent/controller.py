"""Bounded pilot → train/validation → published model lock → untouched test phase."""
from pathlib import Path
import sys,os,json,subprocess,datetime,argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from record_history import sources,validate,name_for,sha,save_json
from fit_models import assert_test_absent,load_models
HERE=Path(__file__).resolve().parent


def review_gate(p,test_phase=False):
    review=json.loads((HERE/'launch-review.json').read_text())
    if review['protocol_sha256']!=sha(HERE/'protocol.json'):raise ValueError('Published protocol changed')
    if not review['predecessor_verified_published_checks_success'] or not review['preregistration_published_checks_success']:
        raise ValueError('Predecessor/preregistration publication not verified')
    predecessor=HERE.parent/'full-history-refinement'
    if (predecessor/'controller.json').exists():raise ValueError('Predecessor still owns numerical work')
    completion=json.loads((predecessor/'completion.json').read_text())
    if not completion['complete'] or sha(predecessor/'completion.json')!=review['predecessor_completion_sha256']:
        raise ValueError('Predecessor completion changed')
    if test_phase:
        lock=json.loads((HERE/'test-lock.json').read_text());publication=json.loads((HERE/'model-lock-publication.json').read_text())
        if not publication['checks_success'] or publication['test_lock_sha256']!=sha(HERE/'test-lock.json'):
            raise ValueError('Frozen model lock not independently published/checked')
        manifest=HERE/lock['manifest']
        if sha(manifest)!=lock['manifest_sha256']:raise ValueError('Fitted manifest changed')
        load_models(manifest)
    else:assert_test_absent(p)


def script(file,*args):
    name=file.replace('.py','')+'-'+('-'.join(str(a) for a in args).replace('--',''))+'.log'
    with (HERE/name).open('w') as log:
        c=subprocess.run([sys.executable,'-B','-u',str(HERE/file),*map(str,args)],stdout=log,stderr=subprocess.STDOUT,
                         creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    if c.returncode:raise RuntimeError(file+' failed; preserved '+name)


def job(seed,c,frozen):
    if sources()!=frozen:raise ValueError('Active source changed')
    name=name_for(seed,c);path=HERE/'recorded-data'/(name+'.npz')
    if not path.exists():
        with (HERE/(name+'.log')).open('w') as stream:
            result=subprocess.run([sys.executable,'-B','-u',str(HERE/'record_history.py'),'--seed',str(seed),
                '--grid',str(c['n']),'--dt',str(c['dt']),'--nu',str(c['nu'])],stdout=stream,stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        if result.returncode:raise RuntimeError(name+' failed; checkpoint preserved')
    meta=validate(path,c)
    return {'name':name,'seed':seed,'config':c,'sha256':meta['sha256'],'seconds':meta['seconds']}


def batch(jobs,p,frozen,kind):
    saved=HERE/'execution.json';execution=json.loads(saved.read_text()) if saved.exists() else {'records':[],'sources':frozen}
    failures=[]
    with ThreadPoolExecutor(max_workers=p['workers']) as pool:
        pending={pool.submit(job,seed,c,frozen):(seed,c) for seed,c in jobs}
        for future in as_completed(pending):
            seed,c=pending[future]
            try:
                record=future.result();record['kind']=kind
                if not any(r['name']==record['name'] for r in execution['records']):execution['records'].append(record)
                if len(execution['records'])>p['max_new_simulations']:raise ValueError('Simulation budget exceeded')
                save_json(saved,execution);print(record['name']+' verified',flush=True)
            except Exception as exc:failures.append({'seed':seed,'config':c,'error':repr(exc)})
    if failures:
        save_json(HERE/(kind+'-failures.json'),{'failures':failures,'sources':frozen})
        raise RuntimeError('Batch failures preserved after unrelated runs finished')


def pilot_jobs(grids,p):
    ret=[(s,{**p['base'],'n':n}) for n in grids for s in p['pilot']['seeds']]
    ret += [(s,{**p['base'],'n':n,'dt':p['pilot']['half_step_dt']}) for n in grids if n in p['pilot']['half_step_grids'] for s in p['pilot']['half_step_seeds']]
    return ret


def audit(n,mode):
    script('numeric_audit.py','--grid',n,'--mode',mode)
    return json.loads((HERE/f'{mode}-audit-{n}.json').read_text())['all_screens_pass']


def stop(reason):
    save_json(HERE/'numerical-stop.json',{'bounded_stop':True,'reason':reason,'protocol_sha256':sha(HERE/'protocol.json')})
    print('Numerical gate failed; no untouched test simulations launched',flush=True)


def preparation(p,frozen):
    batch(pilot_jobs(p['pilot']['initial_grids'],p),p,frozen,'pilot-initial')
    grid=80
    if not audit(80,'pilot'):
        batch(pilot_jobs([160],p),p,frozen,'pilot-conditional')
        if not audit(112,'pilot'):stop('Neither predeclared candidate grid passes pilot accuracy');return
        grid=112
    for candidate in ([80,112] if grid==80 else [112]):
        if candidate==112 and grid==80:
            batch(pilot_jobs([160],p),p,frozen,'pilot-conditional')
            if not audit(112,'pilot'):stop('Higher candidate also fails after fitted-model accuracy failure');return
        jobs=[(s,{**p['base'],'n':candidate}) for s in p['train_seeds']+p['validation_seeds']]
        batch(jobs,p,frozen,f'train-validation-{candidate}')
        script('fit_models.py','--grid',candidate)
        if not audit(candidate,'model'):continue
        assert_test_absent(p)
        manifest=HERE/f'models-{candidate}/manifest.json'
        lock={'grid':candidate,'manifest':str(manifest.relative_to(HERE)).replace('\\','/'),
            'manifest_sha256':sha(manifest),'protocol_sha256':sha(HERE/'protocol.json'),'sources':frozen,
            'locked_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'test_absent_at_lock':True,
            'model_numeric_audit_sha256':sha(HERE/f'model-audit-{candidate}.json')}
        with (HERE/'test-lock.json').open('x') as stream:json.dump(lock,stream,indent=2)
        save_json(HERE/'waiting-model-publication.json',{'status':'Publish and independently verify models/lock before --test-phase','test_lock_sha256':sha(HERE/'test-lock.json')})
        print('All models locked; untouched tests remain absent until publication verification',flush=True)
        return
    stop('Fresh fitted predictors fail the predefined numerical screens at both candidates')


def tests(p,frozen):
    lock=json.loads((HERE/'test-lock.json').read_text());n=lock['grid'];finer={80:112,112:160}[n]
    jobs=[(s,{**p['base'],'n':n}) for s in p['test_seeds']]
    jobs += [(s,{**p['base'],'n':finer}) for s in p['test_numerical_controls']['seeds']]
    jobs += [(s,{**p['base'],'n':n,'dt':p['test_numerical_controls']['half_dt']}) for s in p['test_numerical_controls']['seeds']]
    jobs += [(s,{**p['base'],'n':n,'nu':nu}) for nu in p['viscosity_controls']['nu'] for s in p['viscosity_controls']['seeds']]
    batch(jobs,p,frozen,'untouched-tests-and-controls')
    script('numeric_audit.py','--grid',n,'--mode','test')
    script('analyze_test.py');script('build_report.py');script('verify_completion.py')
    save_json(HERE/'completion.json',{'complete':True,'sources':frozen,'protocol_sha256':sha(HERE/'protocol.json'),
        'results_sha256':sha(HERE/'results.json'),'report_sha256':sha(HERE/'README.md'),
        'verification_sha256':sha(HERE/'verification.json'),'publication_pending':True})
    print('Fresh bounded study complete; verified final publication pending',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--test-phase',action='store_true');args=ap.parse_args()
    p=json.loads((HERE/'protocol.json').read_text());review_gate(p,args.test_phase);frozen=sources()
    lock=HERE/'controller.json'
    with lock.open('x') as stream:json.dump({'pid':os.getpid(),'sources':frozen,'status':'running','phase':'test' if args.test_phase else 'preparation'},stream)
    try:tests(p,frozen) if args.test_phase else preparation(p,frozen)
    except Exception as exc:
        save_json(HERE/'controller-failure.json',{'error':repr(exc),'sources':frozen});raise
    finally:lock.unlink()


if __name__=='__main__':main()
