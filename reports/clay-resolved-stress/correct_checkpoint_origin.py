"""Transparent editorial correction: local checkpoint attribution, raw untouched.

The first local continuation inherited a model-baseline commit label from the
published-checkpoint driver. Preserve that metadata, correct only the origin,
and record before/after metadata hashes. No numerical state or source changes.
"""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=ROOT/'solver/data/extend120-two.json';raw=p.with_suffix('.npz');assert p.exists();r=json.loads(p.read_text());assert r['parent']['upstream_commit']=='fbf4ab15bc1e0875929955e5f1720748464768b6'
 target=ROOT/r['parent']['metadata'];source=ROOT/'solver/data/fine120-two.json';assert target.read_bytes()==source.read_bytes() and sha(target)==r['parent']['metadata_sha256'] and sha(target.with_suffix('.npz'))==r['parent']['raw_sha256']
 history=ROOT/'history';history.mkdir(exist_ok=True);archive=history/'extend120-two-original-metadata.json';assert not archive.exists();archive.write_bytes(p.read_bytes());before=sha(p);state=sha(raw)
 r['parent'].pop('upstream_commit');r['parent'].update(checkpoint_origin='local-recording',source_record='solver/data/fine120-two.json')
 r['metadata_correction']=dict(kind='local-checkpoint-origin',original_metadata='history/'+archive.name,original_metadata_sha256=before,raw_sha256_unchanged=state)
 p.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');assert sha(raw)==state
 correction=dict(record='solver/data/'+p.name,before_metadata_sha256=before,after_metadata_sha256=sha(p),raw_sha256=state,reason='The checkpoint was computed in this study; the former commit label identified the model baseline, not this input recording.')
 (ROOT/'data/provenance-correction.json').write_text(json.dumps(correction,indent=2)+'\n');print(json.dumps(correction,indent=2))
if __name__=='__main__':main()
