"""Execute prerequisites and frozen matrix; preserve failures explicitly."""
import json,traceback,subprocess,sys
from final_experiment import ROOT,dump,digest
def main():
    tests=subprocess.run([sys.executable,'-m','unittest','test_integer_waves','test_corrected_controls','test_corrected_common','-v'],cwd=ROOT,capture_output=True,text=True)
    dump(ROOT/'results'/'corrected-unit-tests.json',dict(exit_code=tests.returncode,output=tests.stdout+tests.stderr,request_id='2300ecd1-d208-4695-9b39-62b93d705892'))
    print(tests.stdout+tests.stderr,flush=True)
    if tests.returncode:raise RuntimeError('Regression tests failed')
    from validate_corrected import validate
    validate()
    from corrected_pilot import qualify
    qualify()
    from production import execute,tasks
    for stage in ['main','resolution_audit','translation','common_probe','common_resolution','unforced_readout','common_unforced_readout']:
        print('START '+stage,flush=True);execute(stage);print('FINISH '+stage,flush=True)
    dump(ROOT/'results'/'gpu-pipeline-complete.json',dict(status='complete',controller_sha256=digest(ROOT/'final_gpu_pipeline.py'),counts={k:len(v) for k,v in tasks().items()}))
if __name__=='__main__':
    try:main()
    except Exception:
        dump(ROOT/'results'/'gpu-pipeline-failure.json',dict(status='failed',traceback=traceback.format_exc()))
        raise
