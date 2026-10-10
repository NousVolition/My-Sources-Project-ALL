"""Run the declared pilot and regenerate its tables, plots and report."""
import subprocess,sys
from pathlib import Path
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
    import runpy
    import pytest
    if pytest.main([str(ROOT/'test_extensions.py'),'-q','-p','no:cacheprovider','--junitxml='+str(ROOT/'results'/'tests.xml')]):
        raise SystemExit('Automated model validation failed')
    # In-process scripts also work through the supplied offline runtime wrapper.
    with threadpool_limits(limits=1):
        for name in ('recovery.py','network.py','recurrence.py','make_report.py','verify.py'):
            print('Running '+name,flush=True);sys.argv=[str(ROOT/name)];runpy.run_path(str(ROOT/name),run_name='__main__')
