"""Reproduce the prospective design in a fresh copy, preserving the freeze gate."""
import argparse,subprocess,sys,concurrent.futures
from common import ROOT,protocol,save

def call(args):
    subprocess.run([sys.executable,*args],cwd=ROOT,check=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['pilots','calibration','confirmation']);ap.add_argument('--workers',type=int,default=1);a=ap.parse_args()
    if a.stage=='pilots':
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
            list(pool.map(lambda i:call(['pilot.py',str(i)]),range(len(protocol()['pilot_cases']))))
        from checks import screening
        s=screening();save('results/pilot_response_checks.json',s);assert s['passed'],s
    elif a.stage=='calibration':
        from checks import screening
        assert screening()['passed']
        call(['batch.py','calibration','--workers',str(a.workers)])
        call(['forecast.py','freeze'])
    else:
        from forecast import verify_freeze
        verify_freeze()
        call(['batch.py','test','--workers',str(a.workers)])
        call(['batch.py','controls','--workers',str(a.workers)])
        call(['checks.py']);call(['forecast.py','evaluate']);call(['analyze.py'])

if __name__=='__main__':main()
