"""Generate and audit the completed report from the original hosted recordings."""
import argparse,subprocess,sys,zipfile,json,os
from pathlib import Path
from common import ROOT,read,save,digest,now
from cloud_stage import unpack

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='downloaded');a=ap.parse_args()
    unpack(a.input)
    from forecast import verify_freeze
    verify_freeze()
    from checks import verify_all
    checks=verify_all(write=False);assert checks==read('results/numerical_checks.json')
    save('results/analysis_receipt.json',dict(started_utc=now(),analysis_commit=os.environ.get('GITHUB_SHA'),source_acquisition_commit='3a6825e9e5f103f6632ac42bfe359a46751dfc1f',workflow_run_id=os.environ.get('GITHUB_RUN_ID'),selection='All original independent runs and dependent controls; no model refitting or outcome selection.'))
    for script in ['environment.py','execution_summary.py','analyze.py','report.py','manifest.py','audit.py']:
        subprocess.run([sys.executable,script],cwd=ROOT,check=True)
    done=read('results/completion.json')
    print(json.dumps({k:done[k] for k in ['numerical_pass','test_stationarity_pass','history_supported']}),flush=True)
    out=ROOT/'higher_re_finished_snapshot.zip'
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path in sorted(ROOT.rglob('*')):
            rel=path.relative_to(ROOT)
            if not path.is_file() or rel.parts[0] in ['downloaded','__pycache__','.pytest_cache'] or path.suffix=='.zip':continue
            z.write(path,rel.as_posix())
    print(json.dumps({'completed_snapshot':out.name,'size':out.stat().st_size,'sha256':digest(out)}),flush=True)

if __name__=='__main__':main()
