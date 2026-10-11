"""Disjoint grid112 main jobs; seed1024 remains reserved for the controller."""
from pathlib import Path
import argparse
import json
import os
import time
from run_refinement import run, sources, save_json
HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--seed', type=int, required=True)
    args = ap.parse_args(); p = json.loads((HERE/'protocol.json').read_text())
    if args.seed not in p['seeds'][1:]: raise ValueError('Reserved or unauthorized seed')
    name = f's{args.seed}_n112_dt{p["dt"]}'
    control = json.loads((HERE/'running.json').read_text())
    e = json.loads((HERE/'execution.json').read_text())
    if any('_n112_' in j['name'] or '_n160_' in j['name'] for j in e['jobs']):
        raise ValueError('Controller advanced; overlapping launch refused')
    if control['sources'] != sources(): raise ValueError('Active source changed')
    claim = HERE/'parallel-claims'/(name+'.json')
    record = {'name': name, 'pid': os.getpid(), 'started_unix': time.time(),
              'controller_pid': control['pid'], 'sources': sources(), 'status': 'running'}
    with claim.open('x') as handle: json.dump(record, handle)
    try:
        meta = run(args.seed, 112, p['dt'], p)
        record.update(status='complete', sha256=meta['sha256'], finished_unix=time.time())
    except Exception as exc:
        record.update(status='failed', error=repr(exc)); raise
    finally: save_json(claim, record)


if __name__ == '__main__': main()
