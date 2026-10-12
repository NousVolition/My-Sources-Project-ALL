"""Compare fresh replay recordings with this package without reusing states."""
from pathlib import Path
import argparse,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('copy',type=Path);ap.add_argument('--only',nargs='+',required=True);a=ap.parse_args();checks=[]
    for name in a.only:
        base=ROOT/'solver/data'/name;replay=a.copy/'solver/data'/name
        expected=json.loads(base.with_suffix('.json').read_text());actual=json.loads(replay.with_suffix('.json').read_text())
        x=np.load(base.with_suffix('.npz'));y=np.load(replay.with_suffix('.npz'))
        assert expected['parameters']==actual['parameters'] and expected['source_hashes']==actual['source_hashes']
        assert not actual['failure'] and np.array_equal(x['time'],y['time'])
        error=float(np.max(abs(x['state']-y['state'])/np.maximum(1,abs(x['state']))))
        checks.append(dict(name=name,source_hashes_equal=True,state_scaled_difference=error,bitwise_equal=bool(np.array_equal(x['state'],y['state'])),passed=error<1e-8))
    result=dict(all_passed=all(c['passed'] for c in checks),checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'data/reproduction-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['all_passed']
if __name__=='__main__':main()
