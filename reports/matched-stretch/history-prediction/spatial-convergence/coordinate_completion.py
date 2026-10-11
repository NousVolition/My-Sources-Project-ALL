"""Await owned jobs, verify caches, apply frozen extension rule and render results."""
from pathlib import Path
import os
import sys
import json
import time
import subprocess
import shutil
from run_refinement import sources, sha, save_json
HERE = Path(__file__).resolve().parent


def step(script, log):
    with (HERE/log).open('w') as handle:
        result = subprocess.run([sys.executable, '-B', str(HERE/script)], stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode: raise RuntimeError(f'{script} exited {result.returncode}; preserved {log}')


def verify_complete(name, source):
    path = HERE/'recorded-data'/(name+'.npz')
    m = json.loads(path.with_suffix('.json').read_text())
    if not m['complete'] or m['sources'] != source or sha(path) != m['sha256']:
        raise ValueError('Changed/incomplete cache '+name)
    if sha(path.with_name(name+'_endpoint-spectrum.npz')) != m['spectrum_sha256']:
        raise ValueError('Changed endpoint spectrum '+name)
    return m


def main():
    p = json.loads((HERE/'protocol.json').read_text()); frozen = sources()
    lock = HERE/'coordinator.json'
    with lock.open('x') as handle: json.dump({'pid': os.getpid(), 'sources': frozen, 'status': 'waiting'}, handle)
    try:
        while True:
            claims = [json.loads(x.read_text()) for x in (HERE/'parallel-claims').glob('*.json')]
            if any(c['status'] == 'failed' for c in claims): raise ValueError('A helper failed; no automatic restart')
            if sources() != frozen: raise ValueError('Active numerical source changed')
            if len(claims) == 9 and all(c['status'] == 'complete' for c in claims) and not (HERE/'running.json').exists(): break
            time.sleep(10)
        initial = [(s, n, p['dt']) for n in p['new_grids'] for s in p['seeds']]
        initial += [(p['time_control_seed'], n, p['half_step_dt']) for n in p['time_control_grids']]
        for s, n, dt in initial: verify_complete(f's{s}_n{n}_dt{dt}', frozen)
        print('All18 initial caches verified; collecting without reintegration', flush=True)
        step('run_refinement.py', 'run-verified-resume.log')
        step('analyze.py', 'analysis-160.log')
        for filename in ['results.json', 'targets.npz', 'convergence.png', 'convergence.svg', 'grid-comparisons.csv']:
            path = HERE/filename; shutil.copyfile(path, path.with_name(path.stem+'-160'+path.suffix))
        decision = json.loads((HERE/'extension-decision.json').read_text())
        if decision['extend']:
            print(f'Frozen screens failed ({len(decision["screen_failures"])}); starting five disjoint224 pairs', flush=True)
            jobs = [(s, p['dt']) for s in p['seeds']]+[(p['time_control_seed'], p['half_step_dt'])]
            children = []
            for seed, dt in jobs:
                name = f's{seed}_n{p["extension_grid"]}_dt{dt}'
                log = (HERE/(name+'.log')).open('w')
                child = subprocess.Popen([sys.executable, '-B', str(HERE/'extension_job.py'), '--seed', str(seed), '--dt', str(dt)],
                          stdout=log, stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW)
                children.append((child, log, name))
            failures = []
            for child, log, name in children:
                code = child.wait(); log.close()
                if code: failures.append({'name': name, 'exit_code': code})
            if failures:
                save_json(HERE/'extension-failures.json', {'failures': failures})
                raise RuntimeError('Extension failures preserved; unrelated pairs allowed to finish')
            records = []
            for seed, dt in jobs:
                m = verify_complete(f's{seed}_n{p["extension_grid"]}_dt{dt}', frozen)
                records.append({k: m[k] for k in ['name', 'sha256', 'seconds', 'complete']})
            save_json(HERE/'extension-execution.json', {'jobs': records, 'complete': True, 'sources': frozen})
            step('analyze.py', 'analysis-final.log')
        step('build_report.py', 'report.log')
        save_json(HERE/'completion.json', {'complete': True, 'sources': frozen,
            'extended': decision['extend'], 'analysis_sha256': sha(HERE/'analyze.py'),
            'results_sha256': sha(HERE/'results.json'), 'report_sha256': sha(HERE/'README.md')})
        print('Numerical matrices, frozen screens and report complete; publication still pending', flush=True)
    except Exception as exc:
        save_json(HERE/'coordinator-failure.json', {'error': repr(exc), 'sources': frozen})
        raise
    finally: lock.unlink()


if __name__ == '__main__': main()
