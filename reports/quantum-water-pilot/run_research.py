"""Full serial reproduction into a fresh directory, with no external publishing."""
import argparse,shutil,subprocess,sys
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--xtb',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--platform',default='OpenCL');a=ap.parse_args()
    root=a.out.resolve();source=Path(__file__).resolve().parent
    if not a.xtb.is_file():raise FileNotFoundError(a.xtb)
    root.mkdir(parents=True,exist_ok=False)
    for name in ['protocol.json','amendment.json','followup-amendments.json','molecular-protocol.json','flow-protocol.json','provenance.json','requirements-lock.txt','later-viscosity-reference.json']:
        shutil.copyfile(source/name,root/name)
    def call(script,*args):
        command=[sys.executable,str(source/script),*map(str,args)];print('RUN',script,flush=True)
        subprocess.run(command,check=True)
    call('run_vibrations.py','--xtb',a.xtb.resolve(),'--out',root/'vibration-final')
    call('check_physics.py','--data',root,'--out',root/'physics-checks.json')
    call('run_molecular.py','--validate','--out',root/'molecular-force-validation.json')
    for stage in ['main','controls','temperature']:
        call('run_molecular_ensemble.py','--out',root/'molecular-data','--stage',stage,'--platform',a.platform)
    call('analyze_molecular.py','--data',root/'molecular-data')
    call('run_fluid.py','--out',root/'fluid-data')
    call('extend_fluid_checks.py','--data',root/'fluid-data')
    call('channel_error_budget.py','--data',root)
    call('stress_more.py','--data',root,'--fluid')
    call('refine_fast_flow.py','--data',root)
    call('compare_later_viscosity.py','--data',root)
    call('build_report.py','--data',root)
    call('verify_results.py','--data',root,'--write-manifest')
    print('Completed. Open',root/'report.html',flush=True)

if __name__=='__main__':main()
