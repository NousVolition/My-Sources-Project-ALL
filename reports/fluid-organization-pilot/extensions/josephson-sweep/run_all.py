"""Run accepted pilot, diagnostics, analysis, figures, tests, then manifest."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent


def main():
    for script in ['run_sweep.py','validate.py','diagnose_precision.py','analyze.py','refine.py','analyze_reference.py','make_report.py']:
        subprocess.run([sys.executable,str(ROOT/script)],check=True,cwd=ROOT)
    subprocess.run([sys.executable,'-m','pytest','test_model.py','-q','-p','no:cacheprovider',
                    '--junitxml=tests.xml'],check=True,cwd=ROOT)
    subprocess.run([sys.executable,str(ROOT/'verify.py')],check=True,cwd=ROOT)


if __name__=='__main__':main()
