"""Run one disjoint future job using the frozen numerical runner, unchanged."""
from pathlib import Path
import argparse
import json
import os
import time
from run_refinement import run, sources, save_json
HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--dt', type=float, required=True)
    a = ap.parse_args(); p = json.loads((HERE/'protocol.json').read_text())
    authorized = [(s, 160, p['dt']) for s in p['seeds']]
    authorized += [(p['time_control_seed'], n, p['half_step_dt']) for n in p['time_control_grids']]
    if (a.seed, a.n, a.dt) not in authorized: raise ValueError('Outside initial parallel assignment')
    name = f's{a.seed}_n{a.n}_dt{a.dt}'
    control = json.loads((HERE/'running.json').read_text())
    execution = json.loads((HERE/'execution.json').read_text())
    # Launch well before the sequential controller can reach the assigned grids.
    if any('_n112_' in j['name'] or '_n160_' in j['name'] for j in execution['jobs']):
        raise ValueError('Controller has advanced; do not launch overlapping work')
    if control['sources'] != sources(): raise ValueError('Active source hash changed')
    claims = HERE/'parallel-claims'; claims.mkdir(exist_ok=True)
    claim = claims/(name+'.json')
    record = {'name': name, 'pid': os.getpid(), 'started_unix': time.time(),
              'controller_pid': control['pid'], 'sources': sources(), 'status': 'running'}
    with claim.open('x') as handle: json.dump(record, handle)
    try:
        meta = run(a.seed, a.n, a.dt, p)
        record.update(status='complete', sha256=meta['sha256'], finished_unix=time.time())
    except Exception as exc:
        record.update(status='failed', error=repr(exc)); raise
    finally:
        save_json(claim, record)


if __name__ == '__main__': main()
