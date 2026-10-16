"""Exercise a fresh restore and verify preservation of matching result files."""
from pathlib import Path
import sys,argparse,subprocess,json,hashlib
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,required=True);ap.add_argument('--destination',type=Path,required=True);ap.add_argument('--write',action='store_true');a=ap.parse_args();dest=a.destination.resolve();assert not dest.exists(),'Use a new destination'
def run(args):
 p=subprocess.run([sys.executable,str(R/'restore_data.py'),*args],capture_output=True,text=True,check=True);return json.loads(p.stdout)
first=run(['--archive',str(a.archive.resolve()),'--destination',str(dest)]);second=run(['--verify-only','--destination',str(dest)]);assert first['passed'] and second['passed'] and first['restored']==first['files'] and second['restored']==0
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();result=dict(passed=True,files=first['files'],fresh_restore=True,matching_files_preserved=True,archive_sha256=first['archive_sha256'],catalog_sha256=sha(R/'data-release.json'),restorer_sha256=sha(R/'restore_data.py'),check_sha256=sha(Path(__file__)))
out=R/'data/archive-restore-check.json'
if a.write or not out.exists():out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
