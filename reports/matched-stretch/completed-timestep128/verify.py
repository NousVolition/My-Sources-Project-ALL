"""Remeasure completed 128-grid timestep controls without evolving fields.

Usage: python verify.py PATH_TO_LOCAL_MATCHED_STRETCH_STUDY
Saved restart arrays are local and are required for this audit.
"""
from pathlib import Path
from datetime import datetime, timezone
import sys, json, hashlib, gc
import numpy as np

HERE=Path(__file__).resolve().parent
STUDY=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'matched-stretch-study'
sys.path.insert(0,str(STUDY))
from numerics import Flow
IDS=['baseline-n128-base','baseline-n128-half']
def read(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(STUDY/'protocol.json'); hashes=protocol['source_hashes']
assert all(sha(STUDY/k)==v for k,v in hashes.items())
records={rid:read(STUDY/'runs'/rid/'result.json') for rid in IDS}
for rid,r in records.items():
    assert r['status']=='complete' and len(r['series'])==41 and r['series'][-1]['t']==.4
    assert r['source_hashes']==hashes and r['job']['n']==128 and r['job']['nu']==.001
    assert r['max_speed_cfl']<.5 and max(abs(z['energy_budget_relative_residual']) for z in r['series'])<.005
assert records[IDS[0]]['dt']==2*records[IDS[1]]['dt']
initial={rid:sha(STUDY/'runs'/rid/'field-000.npy') for rid in IDS}
assert len(set(initial.values()))==1
n=128; dx=6/n; axis=np.arange(n)*dx-3; flow=Flow(n,.001,workers=1)
h0=np.load(STUDY/'runs'/IDS[0]/'field-000.npy',allow_pickle=False)
assert np.array_equal(h0,flow.initial()); del h0
ks=[2*np.pi*np.fft.fftfreq(n,d=dx),2*np.pi*np.fft.fftfreq(n,d=dx),2*np.pi*np.fft.rfftfreq(n,d=dx)]
k=[ks[0][:,None,None],ks[1][None,:,None],ks[2][None,None,:]]
checks=[]; differences=[]; endpoints={}; max_error=0.; mean_drift=0.
for i in range(41):
    velocities=[]
    for rid in IDS:
        r=records[rid]; row=r['series'][i]; assert abs(row['t']-.01*i)<1e-12
        path=STUDY/'runs'/rid/f'field-{i:03d}.npy'; h=np.load(path,allow_pickle=False)
        assert np.isfinite(h).all()
        wh=1j*np.stack([k[1]*h[2]-k[2]*h[1],k[2]*h[0]-k[0]*h[2],k[0]*h[1]-k[1]*h[0]])
        u=np.fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1)); w=np.fft.irfftn(wh,s=(n,)*3,axes=(-3,-2,-1))
        mag=np.linalg.norm(w,axis=0); peak=np.unravel_index(np.argmax(mag),mag.shape)
        measured={'Wmax':float(mag.max()),'energy':float(.5*np.sum(u*u)*dx**3),'enstrophy':float(.5*np.sum(w*w)*dx**3)}
        errors={key:abs(value-row[key])/max(abs(row[key]),1) for key,value in measured.items()}
        assert max(errors.values())<1e-11; max_error=max(max_error,max(errors.values()))
        assert [float(axis[j]) for j in peak]==row['peak_location']
        mean_drift=max(mean_drift,float(np.max(abs(u.mean(axis=(1,2,3))-r['series'][0]['mean_velocity']))))
        checks.append({'run':rid,'t':row['t'],'field_sha256':sha(path),'independent':measured,'relative_errors':errors})
        if i==40:
            observed=flow.observe(h,.4,r['accum'],r['series'][0])
            keys=['high_band_energy_fraction','high_band_enstrophy_fraction','enstrophy_production','energy_budget_relative_residual','enstrophy_budget_relative_residual','stretch_parallel_at_Wmax','viscosity_parallel_at_Wmax']
            ep={key:abs(observed[key]-row[key])/max(abs(row[key]),1) for key in keys}
            assert max(ep.values())<1e-10
            for key in ['width_at_global_peak','width_at_central_roi_peak']:
                assert abs(observed[key]['minimum_chord_cells']-row[key]['minimum_chord_cells'])<1e-10
            _,rates,_=flow.rhs(h); assert abs(rates[3]-row['enstrophy_production'])/max(abs(row['enstrophy_production']),1)<1e-10
            endpoints[rid]=ep
        velocities.append(u); del h,wh,w,mag
    differences.append({'t':i*.01,'velocity_relative_L2_to_half':float(np.linalg.norm(velocities[0]-velocities[1])/np.linalg.norm(velocities[1]))})
    del velocities;gc.collect()
    if i%10==0: print('Checked both fields through '+str(i*.01),flush=True)
assert mean_drift<1e-10
comparisons={}
for key in ['Wmax','I','energy','enstrophy','max_to_mean_spin']:
    a=np.array([z[key] for z in records[IDS[0]]['series']]); b=np.array([z[key] for z in records[IDS[1]]['series']])
    comparisons[key]={'relative_curve_l2':float(np.linalg.norm(a-b)/np.linalg.norm(b)),'relative_max_error':float(np.max(abs(a-b))/np.max(abs(b))),'relative_endpoint_difference':float(abs(a[-1]-b[-1])/abs(b[-1]))}
a=np.array([z['Wmax'] for z in records[IDS[0]]['series']]);b=np.array([z['Wmax'] for z in records[IDS[1]]['series']])
first=next(i for i in range(1,41) if max(abs(a[:i+1]-b[:i+1]))/max(b[:i+1])>.01)
assert first==31
assert all(sha(STUDY/k)==v for k,v in hashes.items())
result={'status':'passed saved-result verification; timestep and spatial screens failed','checked_utc':datetime.now(timezone.utc).isoformat(),'new_simulations':0,'through':.4,'new_half_step_fields_checked':41,'ordinary_fields_rechecked':41,'saved_fields_checked':82,'source_hashes':hashes,'result_hashes':{rid:sha(STUDY/'runs'/rid/'result.json') for rid in IDS},'same_initial_arrays':True,'initial_hashes':initial,'prescribed_initial_array_exact_match':True,'mean_max_drift':mean_drift,'max_relative_remeasurement_error':max_error,'endpoint_diagnostic_errors':endpoints,'first_saved_timestep_screen_failure':first*.01,'timestep_screen_threshold':.01,'comparisons':comparisons,'velocity_comparisons':differences,'fields':checks,'scope':'Physical-space W, energy, enstrophy, mean and peak location independently remeasured at all saved times using an independently written Fourier curl. Frozen diagnostics check endpoint widths, spectra and instantaneous production. Recorded stage integrals are retained without trajectory reruns. The raw strain has nonsmooth periodic joins and initial maxima outside the central tube; spatial convergence remains unestablished.'}
(HERE/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ['fields','velocity_comparisons','endpoint_diagnostic_errors']}))
