"""One disjoint grid224 job, enabled only by the frozen screen decision."""
from pathlib import Path
import argparse
import json
import os
import time
from run_refinement import run, sources, sha, save_json
HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--dt', type=float, required=True); args = ap.parse_args()
    p = json.loads((HERE/'protocol.json').read_text())
    decision = json.loads((HERE/'extension-decision.json').read_text())
    allowed = [(s, p['dt']) for s in p['seeds']]+[(p['time_control_seed'], p['half_step_dt'])]
    if not decision['extend'] or decision['frozen_protocol_sha256'] != sha(HERE/'protocol.json'):
        raise ValueError('Extension decision absent or changed')
    if (args.seed, args.dt) not in allowed or (HERE/'running.json').exists():
        raise ValueError('Unauthorized or overlapping extension')
    name = f's{args.seed}_n{p["extension_grid"]}_dt{args.dt}'
    claims = HERE/'extension-claims'; claims.mkdir(exist_ok=True)
    claim = claims/(name+'.json')
    record = {'name': name, 'pid': os.getpid(), 'started_unix': time.time(),
              'sources': sources(), 'status': 'running'}
    with claim.open('x') as handle: json.dump(record, handle)
    try:
        m = run(args.seed, p['extension_grid'], args.dt, p)
        record.update(status='complete', sha256=m['sha256'], finished_unix=time.time())
    except Exception as exc:
        record.update(status='failed', error=repr(exc)); raise
    finally: save_json(claim, record)


if __name__ == '__main__': main()
