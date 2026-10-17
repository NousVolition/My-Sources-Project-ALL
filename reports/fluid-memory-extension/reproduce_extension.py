"""Explicit full-study rerun, or analysis from verified released recordings."""
import argparse,subprocess,sys,json,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def command(script,*args):subprocess.run([sys.executable,script,*args],cwd=ROOT,check=True)

def preserve_metadata_without_recordings():
    # The GitHub source tree includes small metadata files, not the large NPZ
    # arrays. Such files must not be mistaken for a usable simulation cache.
    from verify_extension import FOLDERS
    missing=[p for folder in FOLDERS for p in (ROOT/folder).glob('*.json') if not p.with_suffix('.npz').exists()]
    if not missing:return
    saved=ROOT/('metadata-only-cache-'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    if saved.exists():raise RuntimeError('Metadata preservation destination already exists')
    for p in missing:
        if not p.resolve().is_relative_to(ROOT.resolve()):raise RuntimeError('Metadata path escapes study directory')
        destination=saved/p.relative_to(ROOT);destination.parent.mkdir(parents=True,exist_ok=True)
        p.rename(destination)
    print(f'Preserved {len(missing)} metadata-only files separately; missing recordings will be computed.',flush=True)
def main():
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True);mode.add_argument('--run',action='store_true');mode.add_argument('--analysis-only',action='store_true');a=ap.parse_args()
    if a.run:
        # GPU production first; at most two independent CPU checks run together.
        preserve_metadata_without_recordings()
        command('final_gpu_pipeline.py');command('parallel_cpu_pipeline.py')
        if (ROOT/'fresh-grid-refinement.json').exists():
            command('refinement_gpu_pipeline.py');command('validate_refined_observation.py')
    command('-m','unittest','test_integer_waves','test_corrected_controls','test_corrected_common','test_inference','test_spatial','-v')
    if (ROOT/'fresh-grid-refinement.json').exists():command('final_analysis.py')
    else:
        command('followup_analysis.py');command('report_extension.py');command('verify_extension.py')
if __name__=='__main__':main()
