"""Remeasure saved viscosity-sweep fields without advancing any trajectory.

Usage: python verify.py /path/to/matched-stretch-study
Requires the original local saved fields, NumPy, SciPy, and frozen numerics.py.
"""
from pathlib import Path
from datetime import datetime, timezone
import sys,json,hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
STUDY=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(STUDY))
from numerics import Flow
RID='viscosity-n128-nu0.01'
def sha(p):
    d=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):d.update(b)
    return d.hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))

def main():
    folder=STUDY/'runs'/RID;r=read(folder/'result.json');p=read(STUDY/'protocol.json')
    assert r['status']=='complete' and len(r['series'])==41 and r['last_output']==40
    assert r['source_hashes']==p['source_hashes'] and all(sha(STUDY/n)==v for n,v in p['source_hashes'].items())
    assert r['job']=={'id':RID,'n':128,'family':'viscosity','nu':.01}
    assert r['dt']==.01/205 and r['end_time']==.4 and p['external_force']==0 and p['domain_side']==6
    f=Flow(128,.01,workers=1)
    h0=np.load(folder/'field-000.npy',allow_pickle=False)
    assert np.array_equal(h0,np.load(STUDY/'initial-n128.npy',allow_pickle=False))
    assert np.array_equal(h0,f.initial())
    mean0=h0[:,0,0,0].real/128**3
    fields=[];maxima={};spectra=[]
    for i,row in enumerate(r['series']):
        assert abs(row['t']-.01*i)<1e-12 and not row.get('hard_failure')
        path=folder/f'field-{i:03d}.npy';h=np.load(path,allow_pickle=False)
        assert np.isfinite(h).all()
        u=np.fft.irfftn(h,s=(128,)*3,axes=(-3,-2,-1))
        wh=f.curl(h);w=np.fft.irfftn(wh,s=(128,)*3,axes=(-3,-2,-1))
        mag=np.sqrt(np.sum(w*w,axis=0))
        independent={'Wmax':float(mag.max()),'energy':float(.5*np.sum(u*u)*f.dx**3),'enstrophy':float(.5*np.sum(w*w)*f.dx**3)}
        drift=float(np.max(abs(u.mean(axis=(1,2,3))-mean0)));assert drift<1e-10
        div=np.fft.irfftn(1j*sum(k*c for k,c in zip(f.k,h)),s=(128,)*3,axes=(-3,-2,-1))
        scaled=float(np.max(abs(div))*f.dx/row['peak_speed']);assert scaled<1e-10
        # Independent physical-space viscous losses and stretching production.
        gradloss=0.;production=0.
        for a in range(3):
            for b in range(3):
                grad=np.fft.irfftn(1j*f.k[b]*h[a],s=(128,)*3,axes=(-3,-2,-1))
                production+=float(np.sum(w[a]*w[b]*grad)*f.dx**3)
            gradw=np.fft.irfftn(1j*f.k[a]*wh,s=(128,)*3,axes=(-3,-2,-1))
            gradloss+=float(np.sum(gradw*gradw)*f.dx**3)
        independent.update(enstrophy_production=production,viscous_energy_loss_rate=2*.01*independent['enstrophy'],viscous_enstrophy_loss_rate=.01*gradloss)
        for k,v in independent.items():
            err=abs(v-row[k])/max(abs(row[k]),1.)
            maxima[k]=max(maxima.get(k,0),err);assert err<1e-10,(i,k,err)
        # Frozen observation routine independently re-reads spectral fractions and chord widths.
        accum=np.array([row['I'],row['integrated_viscous_energy_loss'],0.,row['integrated_enstrophy_production'],row['integrated_viscous_enstrophy_loss']])
        re=f.observe(h,row['t'],accum,r['series'][0])
        for k in ('high_band_energy_fraction','high_band_enstrophy_fraction','energy_budget_relative_residual','enstrophy_budget_relative_residual'):
            err=abs(re[k]-row[k]);maxima[k]=max(maxima.get(k,0),err);assert err<1e-10,(i,k,err)
        for k in ('width_at_global_peak','width_at_central_roi_peak'):
            assert np.isclose(re[k]['minimum_chord_cells'],row[k]['minimum_chord_cells'],rtol=0,atol=1e-10)
        assert abs(row['energy_budget_relative_residual'])<.005
        maxima['mean_drift']=max(maxima.get('mean_drift',0),drift)
        maxima['scaled_divergence']=max(maxima.get('scaled_divergence',0),scaled)
        fields.append({'t':row['t'],'file':path.name,'sha256':sha(path),'physical_measurements':independent})
        if i in (0,10,20,30,40):
            kmag=np.sqrt(f.k2);edges=np.arange(0,np.ceil(kmag.max())+2)
            energy=.5*(h.conj()*h).real.sum(axis=0)*f.weights[0]*f.scale
            ens=.5*(wh.conj()*wh).real.sum(axis=0)*f.weights[0]*f.scale
            spectra.append({'t':row['t'],'k_edges':edges.tolist(),'energy':np.histogram(kmag,bins=edges,weights=energy)[0].tolist(),'enstrophy':np.histogram(kmag,bins=edges,weights=ens)[0].tolist()})
        print(f'Remeasured {i+1}/41 at t={row["t"]:.2f}',flush=True)
    assert r['max_speed_cfl']<.5 and all(sha(STUDY/n)==v for n,v in p['source_hashes'].items())
    report={'status':'passed','checked_utc':datetime.now(timezone.utc).isoformat(),'run':RID,'through':.4,'fields_checked':41,'new_simulations':0,'initial_field_prescribed':True,'initial_mean':mean0.tolist(),'source_hashes':p['source_hashes'],'result_sha256':sha(folder/'result.json'),'protocol_sha256':sha(STUDY/'protocol.json'),'max_errors':maxima,'fields':fields,'spectra':spectra,'scope':'Independent NumPy physical-space FFT measurements at every output; frozen width/spectral algorithm remeasurement. Saved Heun stage integrals checked against endpoint budgets; trajectories were not repeated.'}
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Verification passed',flush=True)
if __name__=='__main__':main()
