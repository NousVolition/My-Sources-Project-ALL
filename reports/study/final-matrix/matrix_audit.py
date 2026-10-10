"""Verify completed matrix records and checkpoint identities; perform no evolution."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, zipfile, sys
import numpy as np

HERE=Path(__file__).resolve().parent
WORK=HERE.parent
STUDY=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else WORK/'adversarial-vortex-study'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def blob(p):
    b=p.read_bytes();return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def finite(v):
    if isinstance(v,dict):return all(finite(x) for x in v.values())
    if isinstance(v,list):return all(finite(x) for x in v)
    return bool(np.isfinite(v)) if isinstance(v,(int,float)) else True
prot=read(STUDY/'protocol.json');state=read(STUDY/'task-state.json')
base=read(HERE/'base-blob-map.json');hashes=prot['source_hashes']
for name,h in hashes.items():assert sha(STUDY/name)==h
assert state['status']=='complete' and len(state['completed'])==48
jobs={j['id']:j for j in prot['jobs']};assert len(jobs)==48
records={p.parent.name:read(p) for p in (STUDY/'runs').glob('*/result.json')}
assert records.keys()==jobs.keys()
checks=[];summary={}
for rid,d in sorted(records.items()):
    folder=STUDY/'runs'/rid
    assert d['status']=='complete' and d['job']==jobs[rid] and d['source_hashes']==hashes
    rows=d['rows'];assert len(rows)==21 and finite(rows)
    assert np.allclose([r['t'] for r in rows],np.arange(21)*.02,rtol=0,atol=1e-12)
    assert d['maximum_stage_cfl']<.75
    assert max(abs(r['energy_balance_relative']) for r in rows)<.005
    assert max(r['scaled_native_divergence'] for r in rows)<1e-10
    assert max(abs(x) for r in rows for x in r['mean_velocity'])<1e-10
    assert np.all(np.diff([r['I'] for r in rows])>=0)
    published='reports/study/runs/'+rid+'/result.json'
    old=published in base['blobs']
    if old:assert blob(folder/'result.json')==base['blobs'][published]
    else:assert rid=='exodus-fd4-n160-half'
    cp=folder/'checkpoint.npz';last=folder/'field-t0.40.npz'
    cp_sha=sha(cp);assert cp_sha==sha(last)
    with np.load(cp,allow_pickle=False) as a:
        assert float(a['t'])==.4
        assert json.loads(str(a['job_json'].item()))==jobs[rid]
        assert json.loads(str(a['source_hashes'].item()))==hashes
    with zipfile.ZipFile(cp) as z:
        with z.open('h.npy') as f:
            ver=np.lib.format.read_magic(f)
            shape,order,dtype=(np.lib.format.read_array_header_1_0(f) if ver==(1,0) else np.lib.format.read_array_header_2_0(f))
    n=jobs[rid]['n'];assert shape==(3,n,n,n//2+1) and dtype==np.complex128
    checks.append({'id':rid,'result_sha256':sha(folder/'result.json'),'previous_published_blob_verified':old,'final_checkpoint_sha256':cp_sha,'checkpoint_equals_final_field':True,'shape':shape,'dtype':str(dtype)})
    w=np.array([r['Wmax'] for r in rows]);changes=np.diff(w)
    peaks=[i for i in range(1,len(w)-1) if w[i]>w[i-1] and w[i]>w[i+1]]
    summary[rid]={'job':jobs[rid],'endpoint':rows[-1],'peak_W':float(max(w)),'peak_saved_time':rows[int(np.argmax(w))]['t'],'minimum_global_width_cells':min(r['global_width']['minimum_chord_cells'] for r in rows),'maximum_high_band_fraction':max(r['high_band_enstrophy_fraction'] for r in rows),'first_spectral_warning':next((r['t'] for r in rows if r['high_band_enstrophy_fraction']>.01),None),'first_width_warning':next((r['t'] for r in rows if r['global_width']['minimum_chord_cells']<6),None),'max_abs_energy_budget_relative':max(abs(r['energy_balance_relative']) for r in rows),'max_abs_native_enstrophy_budget_relative':max(abs(r['enstrophy_balance_relative']) for r in rows),'maximum_stage_cfl':d['maximum_stage_cfl'],'minimum_marker_distance_start':min(x['minimum'] for x in rows[0]['outer_marker_distances']),'minimum_marker_distance_end':min(x['minimum'] for x in rows[-1]['outer_marker_distances']),'recurrence_diagnostic':{'saved_samples':21,'strict_interior_local_maxima_times':[rows[i]['t'] for i in peaks],'successive_W_pairs':[[float(a),float(b)] for a,b in zip(w[:-1],w[1:])],'scope':'Consecutive scalar samples and sampled turning points only; not a Poincare section, field recurrence, period estimate or asymptotic Lyapunov exponent.'}}
    print('Checked completed record and final checkpoint: '+rid,flush=True)
comparisons=read(STUDY/'comparisons.json');assert len(comparisons)==77
curve_checks=0
for key,c in comparisons.items():
    a,b=records[c['a']],records[c['b']]
    assert c['fingerprint']==hashlib.sha256(json.dumps([a,b],sort_keys=True).encode()).hexdigest()
    for metric in ['Wmax','W_central_roi','I','energy','spin_ratio','origin_axial_strain','central_marker_parallel_stretch_spin_weighted']:
        av=np.array([r[metric] for r in a['rows']]);bv=np.array([r[metric] for r in b['rows']])
        expected=float(max(abs(av-bv))/max(max(abs(bv)),1e-300))
        assert np.isclose(expected,c['diagnostics'][metric]['max_curve_difference_relative_to_reference_peak'],rtol=1e-12,atol=1e-15)
        curve_checks+=1
out={'status':'passed','checked_utc':datetime.now(timezone.utc).isoformat(),'new_simulations':0,'completed_runs':48,'saved_observations':48*21,'previous_published_records_hash_verified':47,'final_checkpoints_verified':48,'comparison_records_verified':77,'scalar_comparisons_recomputed':curve_checks,'source_hashes':hashes,'base_commit':base['head'],'checks':checks,'scope':'All completed result records match protocol, times and numerical gates; prior 47 result blobs match the verified publication base. All final checkpoint metadata and headers are checked and hashes equal final archived fields. All 77 comparison fingerprints and seven scalar comparisons per pair are recomputed. Earlier full-field audits are reused; only the newly completed FD4 run receives a new five-field physical audit and five-pair timestep audit in companion records. This audit does not rerun trajectories or prove continuum convergence.'}
(HERE/'matrix-verification.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
(HERE/'matrix-summary.json').write_text(json.dumps({'through':.4,'completed':48,'planned':48,'remaining':0,'hard_failures':0,'runs':summary},indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k not in ['checks','source_hashes','scope']}))
