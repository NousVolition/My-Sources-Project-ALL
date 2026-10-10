"""Independent same-grid comparisons from saved fields; no time evolution."""
import gc
import sys
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE.parent
STUDY = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else WORK / 'adversarial-vortex-study'
BASE = 'exodus-fd4-n160-base'
HALF = 'exodus-fd4-n160-half'
read = lambda p: json.loads(p.read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

audit = read(HERE/'verification.json')
previous = read(HERE/'verification.json')
assert audit['status'] == previous['status'] == 'passed'
assert audit['source_hashes'] == previous['source_hashes']
records = {rid:read(STUDY/'runs'/rid/'result.json') for rid in [BASE,HALF]}
for rid, proof in [(BASE, previous), (HALF, audit)]:
    check = next(x for x in proof['checks'] if x['id']==rid)
    assert sha(STUDY/'runs'/rid/'result.json') == check['result_sha256']
    assert records[rid]['source_hashes'] == audit['source_hashes']
    assert records[rid]['status'] == 'complete'
assert records[BASE]['job']['dt'] == 2*records[HALF]['job']['dt']
comparison = read(STUDY/'comparisons.json')[BASE+'__'+HALF]
curve_checks=[]
for key in ['Wmax','W_central_roi','I','energy','spin_ratio','origin_axial_strain','central_marker_parallel_stretch_spin_weighted']:
    a=np.array([r[key] for r in records[BASE]['rows']])
    b=np.array([r[key] for r in records[HALF]['rows']])
    difference=float(np.max(abs(a-b))/np.max(abs(b)))
    assert np.isclose(difference,comparison['diagnostics'][key]['max_curve_difference_relative_to_reference_peak'],rtol=1e-12,atol=1e-15)
    curve_checks.append({'diagnostic':key,'max_curve_difference_relative_to_half_peak':difference})

n=160
indices=np.rint(np.fft.fftfreq(n)*n).astype(int)
keep=(3*abs(indices[:,None,None])<n)&(3*abs(indices[None,:,None])<n)&(3*np.arange(n//2+1)[None,None,:]<n)
weights=np.ones(n//2+1)*2;weights[0]=weights[-1]=1
field_checks=[]
for t in [0.,.1,.2,.3,.4]:
    fields=[]; hashes={}
    for rid,proof in [(BASE,previous),(HALF,audit)]:
        path=STUDY/'runs'/rid/f'field-t{t:.2f}.npz'
        expected=next(f['sha256'] for c in proof['checks'] if c['id']==rid for f in c['fields'] if f['file']==path.name)
        hashes[rid]=sha(path);assert hashes[rid]==expected
        with np.load(path,allow_pickle=False) as a:
            assert float(a['t'])==t
            fields.append(a['h'])
    a,b=fields
    for h in fields:
        assert np.max(abs(h[:,~keep]))<1e-12*max(1,np.max(abs(h)))
    numerator=np.sum(abs(a-b)**2*weights[None,None,None,:])
    denominator=np.sum(abs(b)**2*weights[None,None,None,:])
    l2=float(np.sqrt(numerator/denominator))
    ua=np.fft.irfftn(a,s=(n,)*3,axes=(-3,-2,-1))
    ub=np.fft.irfftn(b,s=(n,)*3,axes=(-3,-2,-1))
    linf=float(np.linalg.norm(ua-ub,axis=0).max()/np.linalg.norm(ub,axis=0).max())
    physical_l2=float(np.linalg.norm(ua-ub)/np.linalg.norm(ub))
    assert np.isclose(l2,physical_l2,atol=1e-14,rtol=1e-7)
    if t==0:
        assert np.array_equal(a,b)
    else:
        saved=next(r for r in comparison['fields'] if r['t']==t)
        assert np.isclose(l2,saved['relative_velocity_L2_common_modes'],atol=1e-13,rtol=1e-7)
        assert np.isclose(linf,saved['relative_velocity_Linf_common_modes'],atol=1e-13,rtol=1e-7)
    field_checks.append({'t':t,'relative_velocity_L2':l2,'relative_velocity_Linf':linf,'field_sha256':hashes})
    del a,b,h,ua,ub,fields
    gc.collect()

# Final checkpoint is the identical saved final field, including metadata.
assert sha(STUDY/'runs'/HALF/'checkpoint.npz')==sha(STUDY/'runs'/HALF/'field-t0.40.npz')
out={'status':'passed','new_simulations':0,'initial_fields_bitwise_identical':True,'curve_checks':curve_checks,'field_checks':field_checks,'comparison':comparison,'scope':'Same grid and retained band: weighted RFFT Parseval L2 and physical-space L2/Linf were independently checked at five saved times. Previous ordinary-run audit reused after all five field and result hashes matched. Final checkpoint equals the saved final field.'}
(HERE/'comparison-verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'curves':len(curve_checks),'field_pairs':len(field_checks),'endpoint':field_checks[-1],'status':'passed'}))
