"""Explicit full-study rerun, or analysis from verified released recordings."""
import argparse,subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def command(script,*args):subprocess.run([sys.executable,script,*args],cwd=ROOT,check=True)
def main():
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True);mode.add_argument('--run',action='store_true');mode.add_argument('--analysis-only',action='store_true');a=ap.parse_args()
    if a.run:
        # GPU production first; CPU checks are then reproducible sequentially.
        command('final_gpu_pipeline.py');command('cpu_pipeline.py')
    command('-m','unittest','test_integer_waves','test_corrected_controls','test_corrected_common','test_inference','-v')
    command('followup_analysis.py');command('report_extension.py');command('verify_extension.py')
if __name__=='__main__':main()
