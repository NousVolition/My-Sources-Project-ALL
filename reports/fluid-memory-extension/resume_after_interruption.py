"""Resume the frozen matrix after a stopped host session, then verify it.

Completed recordings are reused only after the original runners verify their
source and raw-data fingerprints. This wrapper changes no scientific settings.
"""
import hashlib
import json
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, value):
    (ROOT / 'results' / name).write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def main():
    freeze = json.loads((ROOT / 'refinement-analysis-freeze.json').read_text(encoding='utf-8'))
    for name, expected in freeze['sources'].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError('Frozen source mismatch: ' + name)
    prior = json.loads((ROOT / 'results' / 'refinement-progress.json').read_text(encoding='utf-8'))
    audit = dict(time_utc=datetime.now(timezone.utc).isoformat(),
                 reason='Continue after the previous host processes stopped during an account-limit interruption; exact stop cause not established.',
                 prior_refinements_completed=len(prior['completed']),
                 planned_refinements=prior['planned'],
                 fresh_memory_effects_inspected=False,
                 scientific_settings_changed=False,
                 frozen_sources_verified=True,
                 source_sha256=digest(Path(__file__)),
                 request_id='2300ecd1-d208-4695-9b39-62b93d705892')
    save('resumption-2026-10-13.json', audit)
    for script in ['refinement_gpu_pipeline.py', 'validate_refined_observation.py', 'final_analysis.py']:
        print('START ' + script, flush=True)
        subprocess.run([sys.executable, '-u', script], cwd=ROOT, check=True)
        print('FINISH ' + script, flush=True)
    save('resumption-complete.json', dict(time_utc=datetime.now(timezone.utc).isoformat(),
                                        source_sha256=digest(Path(__file__)),
                                        status='Verified locally; visual review and GitHub publication remain'))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        save('resumption-failure.json', dict(time_utc=datetime.now(timezone.utc).isoformat(),
                                             traceback=traceback.format_exc()))
        raise
