"""Read completed parameter runs, independently remeasure saved fields, never evolve."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;STUDY=ROOT.parent
sys.path.insert(0,str(STUDY))
from numerics import Flow
IDS=["angle_degrees-n64-v5.0"]
REUSED=["baseline-n64-base","angle_degrees-n64-v-5.0"]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    data={};checks=[];source_checks=[];baseline=np.load(STUDY/'runs/baseline-n64-base/field-000.npy')
    revision=json.loads((STUDY/'revisions/windows-write-retry.json').read_text())
    current={name:sha(STUDY/name) for name in revision['new']}
    assert current==revision['new'] and revision['old']['numerics.py']==current['numerics.py']
    for rid in IDS+REUSED:
        path=STUDY/'runs'/rid/'result.json';r=json.loads(path.read_text());assert r['status']=='complete' and len(r['series'])==41 and r['series'][-1]['t']==.4
        assert r['source_hashes']==current or r['source_hashes']==revision['old']
        source_checks.append({'id':rid,'matches_current_sources':r['source_hashes']==current,'source_hashes':r['source_hashes']})
        data[rid]=r
        if rid not in IDS:continue
        n=r['job']['n'];f=Flow(n,r['job']['nu'],workers=1);fields=[];max_errors=np.zeros(3);mean0=np.array(r['series'][0]['mean_velocity']);max_mean=0
        for i,row in enumerate(r['series']):
            assert abs(row['t']-.01*i)<1e-12
            fp=path.parent/f'field-{i:03d}.npy';h=np.load(fp,allow_pickle=False);assert np.isfinite(h).all()
            u=np.fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1));omega=np.fft.irfftn(f.curl(h),s=(n,)*3,axes=(-3,-2,-1))
            measured=[np.linalg.norm(omega,axis=0).max(),.5*np.sum(u*u)*f.dx**3,.5*np.sum(omega*omega)*f.dx**3]
            expected=[row[k] for k in ('Wmax','energy','enstrophy')];errors=abs(np.array(measured)-expected)/np.maximum(abs(np.array(expected)),1)
            assert errors.max()<1e-11;max_errors=np.maximum(max_errors,errors);max_mean=max(max_mean,float(np.max(abs(np.mean(u,axis=(1,2,3))-mean0))))
            fields.append({'file':fp.name,'sha256':sha(fp)})
        initial=np.load(path.parent/'field-000.npy');kwargs={k:r['job'][k] for k in ['radius','spin_factor','strain_factor','angle_degrees'] if k in r['job']}
        assert np.array_equal(initial,f.initial(**kwargs))
        if r['job']['family']=='viscosity':assert np.array_equal(initial,baseline)
        end=f.observe(h,.4,r['accum'],r['series'][0]);last=r['series'][-1]
        endpoint_keys=['Wmax','energy','enstrophy','enstrophy_production','high_band_enstrophy_fraction','divergence_max','energy_budget_relative_residual','enstrophy_budget_relative_residual']
        errs={k:abs(end[k]-last[k])/max(abs(last[k]),1) for k in endpoint_keys};assert max(errs.values())<1e-10
        _,rates,_=f.rhs(h);production_error=abs(rates[3]-last['enstrophy_production'])/max(abs(last['enstrophy_production']),1);assert production_error<1e-10
        assert np.all(np.diff([x['I'] for x in r['series']])>=0) and r['max_speed_cfl']<.5 and max_mean<1e-10
        assert max(abs(x['energy_budget_relative_residual']) for x in r['series'])<.005
        checks.append({'id':rid,'fields':fields,'result_sha256':sha(path),'initial_field_matches_prescribed_formula':True,'mean_velocity_max_drift':max_mean,'max_relative_W_energy_enstrophy_errors':max_errors.tolist(),'endpoint_errors':errs,'direct_vs_rhs_enstrophy_production_error':production_error})
        print('VERIFIED '+rid+' 41 saved fields',flush=True)
    result={'new_completed_runs':IDS,'reused_completed_runs':REUSED,'runs':data}
    (ROOT/'measurements.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    (ROOT/'verification.json').write_text(json.dumps({'status':'passed','saved_fields_checked':41,'numerical_solver_hash_unchanged':True,'source_checks':source_checks,'historical_metadata_revision':revision,'new_simulations':0,'runs':checks,'note':'I and integrated budget rates retain recorded per-step accumulations; no time evolution was repeated. Records accept only the documented pre/post metadata-write-retry versions. The numerical solver is identical; per-run source hashes are recorded.'},indent=2)+'\n')
    print('Completed plus-five-degree tilt verified.')
if __name__=='__main__':main()
