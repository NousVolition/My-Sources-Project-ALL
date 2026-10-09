"""Remeasure the new 128-grid seed and its separation from saved baseline fields."""
from pathlib import Path
from datetime import datetime, timezone
import sys,json,hashlib,gc
import numpy as np
HERE=Path(__file__).resolve().parent
STUDY=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'matched-stretch-study'
sys.path.insert(0,str(STUDY))
from numerics import Flow
RID='perturb-n128-a0.001-s303'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    r=json.loads((STUDY/'runs'/RID/'result.json').read_text())
    protocol=json.loads((STUDY/'protocol.json').read_text())
    assert r['status']=='complete' and len(r['series'])==41 and r['series'][-1]['t']==.4
    assert r['source_hashes']==protocol['source_hashes']
    assert all(sha(STUDY/k)==v for k,v in r['source_hashes'].items())
    assert r['job']['n']==128 and r['job']['amplitude']==.001 and r['job']['seed']==303
    n=128;flow=Flow(n,r['job']['nu'],workers=1)
    base=STUDY/'runs/baseline-n128-base'
    reference=json.loads((base/'result.json').read_text());assert reference['status']=='complete'
    h0=np.load(base/'field-000.npy',allow_pickle=False)
    initial=np.load(STUDY/'runs'/RID/'field-000.npy',allow_pickle=False)
    assert np.array_equal(initial,flow.perturb(h0,.001,303))
    assert np.array_equal(h0,flow.initial())
    u0=np.fft.irfftn(h0,s=(n,)*3,axes=(-3,-2,-1))
    norm0=float(np.sum(u0*u0));mean0=u0.mean(axis=(1,2,3))
    del initial,h0,u0
    errors_max=np.zeros(4);drift=0.;fields=[];reference_fields=[]
    for i,row in enumerate(r['series']):
        assert abs(row['t']-.01*i)<1e-12
        file=STUDY/'runs'/RID/f'field-{i:03d}.npy';ref=base/f'field-{i:03d}.npy'
        h=np.load(file,allow_pickle=False);hb=np.load(ref,allow_pickle=False)
        assert np.isfinite(h).all() and np.isfinite(hb).all()
        u=np.fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1))
        omega=np.fft.irfftn(flow.curl(h),s=(n,)*3,axes=(-3,-2,-1))
        du=np.fft.irfftn(h-hb,s=(n,)*3,axes=(-3,-2,-1))
        measured=np.array([np.linalg.norm(omega,axis=0).max(),.5*np.sum(u*u)*flow.dx**3,.5*np.sum(omega*omega)*flow.dx**3,np.sqrt(np.sum(du*du)/norm0)])
        expected=np.array([row[k] for k in ['Wmax','energy','enstrophy','D_relative_initial_baseline_l2']])
        errors=abs(measured-expected)/np.maximum(abs(expected),1e-15);assert errors.max()<1e-10
        errors_max=np.maximum(errors_max,errors);drift=max(drift,float(np.max(abs(u.mean(axis=(1,2,3))-mean0))))
        fields.append({'file':file.name,'sha256':sha(file)});reference_fields.append({'file':ref.name,'sha256':sha(ref)})
        del hb,u,omega,du
        if i<40:del h
        if i%10==0:print('Matched seed: checked saved time '+str(row['t']),flush=True)
    assert drift<1e-10 and r['max_speed_cfl']<.5
    assert max(abs(z['energy_budget_relative_residual']) for z in r['series'])<.005
    assert np.all(np.diff([z['I'] for z in r['series']])>=0)
    end=flow.observe(h,.4,r['accum'],r['series'][0]);last=r['series'][-1]
    keys=['Wmax','energy','enstrophy','enstrophy_production','high_band_enstrophy_fraction','divergence_max','energy_budget_relative_residual','enstrophy_budget_relative_residual']
    endpoint={k:abs(end[k]-last[k])/max(abs(last[k]),1) for k in keys};assert max(endpoint.values())<1e-10
    for k in ['width_at_global_peak','width_at_central_roi_peak']:
        assert np.isclose(end[k]['minimum_chord_cells'],last[k]['minimum_chord_cells'],atol=1e-10,rtol=0)
    _,rates,_=flow.rhs(h)
    prod=abs(rates[3]-last['enstrophy_production'])/max(abs(last['enstrophy_production']),1);assert prod<1e-10
    d0=r['series'][0]['D_relative_initial_baseline_l2'];d1=last['D_relative_initial_baseline_l2']
    result={'status':'passed','checked_utc':datetime.now(timezone.utc).isoformat(),'new_simulations':0,
            'run':RID,'new_saved_fields_checked':41,'reference_fields_compared':41,
            'result_sha256':sha(STUDY/'runs'/RID/'result.json'),'baseline_result_sha256':sha(base/'result.json'),
            'source_hashes':r['source_hashes'],'initial_field_matches_prescribed_perturbation':True,
            'initial_relative_perturbation':d0,'mean_max_drift':drift,
            'max_relative_W_energy_enstrophy_D_errors':errors_max.tolist(),'endpoint_errors':endpoint,
            'direct_vs_rhs_enstrophy_production_error':prod,'fields':fields,'reference_fields':reference_fields,
            'endpoint_log_growth_rate':float(np.log(d1/d0)/.4),
            'rate_scope':'Finite-amplitude separation over the recorded window; no renormalization and no asymptotic Lyapunov claim.',
            'scope':'Independent physical-space W, energy, enstrophy and normalized difference from the saved baseline at all 41 output times. Endpoint diagnostics checked with frozen code; stage integrals retained without rerunning trajectories.'}
    assert all(sha(STUDY/k)==v for k,v in r['source_hashes'].items())
    (HERE/'matched-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Matched seed verification passed.',flush=True)
if __name__=='__main__':main()

