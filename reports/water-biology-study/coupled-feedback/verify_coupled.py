"""Verify the delivered package and optionally compare an independent rerun."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from coupled_model import CoupledFlow
from analyze_coupled import embedded, norm, delta

ROOT=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--actual',type=Path);args=ap.parse_args()
    manifest=read(ROOT/'manifest.json')
    for name,expected in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    result=read(ROOT/'data/results.json');analysis=read(ROOT/'data/analysis.json');fields=np.load(ROOT/'data/final-fields.npz')
    expected={r['case'] for r in result}
    assert len(result)==30 and len(expected)==30 and set(fields.files)==expected
    assert analysis['passed'] and all(analysis['gates'].values())
    for r in result:
        h=fields[r['case']];m=CoupledFlow(r['n'],feedback=r['feedback'])
        assert h.shape==(2,3,r['n'],r['n']) and np.isfinite(h).all()
        assert len(r['rows'])==51 and r['rows'][-1]['t']==1.
        for j,name in enumerate(['baseline','disturbed']):
            d=m.diagnostics(h[j]);s=r['rows'][-1][name]
            np.testing.assert_allclose([d[k] for k in d],[s[k] for k in d],rtol=1e-11,atol=1e-12)
        response=float(norm(embedded(delta(h[None])))[0]/r['initial_impulse_L2'])
        assert abs(response-r['rows'][-1]['disturbance_velocity_L2_relative_to_initial_impulse'])<1e-11
    if args.actual:
        actual=read(args.actual/'results.json');af=np.load(args.actual/'final-fields.npz');aa=read(args.actual/'analysis.json')
        assert len(actual)==30 and {r['case'] for r in actual}==expected and set(af.files)==expected
        assert aa['passed'] and all(aa['gates'].values())
        byid={r['case']:r for r in actual}
        for r in result:
            np.testing.assert_allclose(af[r['case']],fields[r['case']],rtol=1e-10,atol=1e-11)
            np.testing.assert_allclose([x['disturbance_velocity_L2_relative_to_initial_impulse'] for x in r['rows']],
                [x['disturbance_velocity_L2_relative_to_initial_impulse'] for x in byid[r['case']]['rows']],rtol=1e-10,atol=1e-11)
        assert len(aa['contrasts'])==len(analysis['contrasts'])
        for a,b in zip(aa['contrasts'],analysis['contrasts']):
            assert (a['n'],a['dt'],a['seed'],a['feedback'])==(b['n'],b['dt'],b['seed'],b['feedback'])
            np.testing.assert_allclose(a['curve'],b['curve'],rtol=1e-9,atol=1e-11)
    print('PASS: hashes, all 30 case identities, endpoint fields, response diagnostics'+(' and independent reproduction' if args.actual else ''))

if __name__=='__main__':main()
