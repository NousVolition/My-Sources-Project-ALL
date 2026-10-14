"""Replay exact protocols in separate processes with recorded BLAS settings."""
import argparse,json,sys,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--only',nargs='+');a=ap.parse_args();protocol=json.loads((ROOT/'data/executed-protocol.json').read_text());records=protocol['attempts']
 if a.only:assert set(a.only)<=set(v['name'] for v in records)
 for r in records:
  if a.only and r['name'] not in a.only:continue
  profile=protocol['runtime_profiles'][r['name']];env=os.environ.copy();env['OPENBLAS_NUM_THREADS']=str(profile['numpy_blas_threads']);env['OMP_NUM_THREADS']='1';env['MKL_NUM_THREADS']='1'
  subprocess.run([sys.executable,str(ROOT/'replay_one.py'),r['name']],env=env,check=True)
if __name__=='__main__':main()
