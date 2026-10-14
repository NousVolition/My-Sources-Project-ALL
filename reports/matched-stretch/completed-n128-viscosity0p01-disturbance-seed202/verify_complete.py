"""Independent completed-pair audit; no integration. Usage: python verify_complete.py PILOT_DIRECTORY [OUTPUT_JSON]."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,hashlib
import numpy as np
PILOT=Path(sys.argv[1]).resolve()
OUTPUT=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else Path(__file__).with_name('verification.json')
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
p=read(PILOT/'protocol.json');STUDY=Path(p['original_study']);REF=Path(p['reference_directory'])
sys.path.insert(0,str(PILOT/'source'))
from numerics import Flow
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
def real(x):return np.fft.irfftn(x,s=(128,)*3,axes=(-3,-2,-1))
maxima={}
def check(key,got,want,tol=1e-10,absolute=False):
 err=abs(float(got)-float(want))/(1 if absolute else max(abs(float(want)),1e-20))
 maxima[key]=max(maxima.get(key,0.),err)
 assert err<tol,(key,got,want,err)
def main():
 r=read(PILOT/'result.json');b=read(REF/'result.json');old=read(PILOT/'saved-field-verification.json')
 assert r['status']==b['status']=='complete' and len(r['series'])==len(b['series'])==old['field_pairs']==41
 assert old['status']=='passed' and old['through']==.4 and r['last_output']==b['last_output']==40
 assert p['grid']==128 and p['viscosity']==.01 and p['amplitude']==.001 and p['seed']==202
 assert p['dt']==r['dt']==b['dt']==.01/205 and p['end_time']==.4 and p['output_dt']==.01
 assert p['workers']==1 and p['external_force']==0 and p['domain_side']==6
 assert sha(PILOT/'protocol.json')==r['protocol_sha256']
 for k,v in p['pilot_source_hashes'].items():assert sha(PILOT/k)==v
 for k,v in p['original_source_hashes'].items():assert sha(STUDY/k)==v
 assert sha(PILOT/'source/numerics.py')==p['original_source_hashes']['numerics.py']
 f=Flow(128,.01,workers=1);initial_b=np.load(REF/'field-000.npy',allow_pickle=False);initial_d=np.load(PILOT/'fields/field-000.npy',allow_pickle=False)
 assert np.array_equal(initial_b,f.initial()) and np.array_equal(initial_b,np.load(STUDY/'initial-n128.npy',allow_pickle=False))
 assert np.array_equal(initial_d,f.perturb(initial_b,.001,202)) and np.array_equal(initial_d,np.load(PILOT/'initial-disturbed.npy',allow_pickle=False))
 assert sha(REF/'field-000.npy')==p['reference_initial_sha256'] and sha(PILOT/'initial-disturbed.npy')==p['initial_disturbed_sha256']
 norm0=float(np.sqrt(np.sum(real(initial_b)**2)*f.dx**3));check('initial_reference_norm',norm0,p['reference_initial_l2'])
 lowp=STUDY/'runs/perturb-n128-a0.001-s202';lowb=STUDY/'runs/baseline-n128-base'
 lr=read(lowp/'result.json');lb=read(lowb/'result.json')
 assert lr['status']==lb['status']=='complete' and len(lr['series'])==len(lb['series'])==41
 assert lr['job']['nu']==lb['job']['nu']==.001
 assert np.array_equal(initial_d,np.load(lowp/'field-000.npy',allow_pickle=False))
 assert np.array_equal(initial_b,np.load(lowb/'field-000.npy',allow_pickle=False))
 mean0=np.array(p['mean_velocity']);records=[];spectra=[];reference_rows=b['series']
 for i,(row,br) in enumerate(zip(r['series'],reference_rows)):
  assert abs(row['t']-i*.01)<1e-12 and row['t']==br['t']==lr['series'][i]['t']==lb['series'][i]['t']
  assert not row.get('hard_failure') and not br.get('hard_failure')
  hd=np.load(PILOT/f'fields/field-{i:03d}.npy',allow_pickle=False);hb=np.load(REF/f'field-{i:03d}.npy',allow_pickle=False)
  delta=real(hd-hb);physical_D=float(np.sqrt(np.sum(delta*delta)*f.dx**3)/norm0)
  check('difference_norm',physical_D,row['D_relative_initial_reference_l2'])
  check('amplification',physical_D/.001,row['amplification'])
  transfer=0.;parts={}
  for role,h,rr,first,h0 in [('disturbed',hd,row,r['series'][0],initial_d),('reference',hb,br,b['series'][0],initial_b)]:
   assert np.isfinite(h).all()
   u=real(h);wh=f.curl(h);w=real(wh)
   measured={'Wmax':float(np.sqrt(np.sum(w*w,axis=0)).max()),'energy':float(.5*np.sum(u*u)*f.dx**3),'enstrophy':float(.5*np.sum(w*w)*f.dx**3)}
   drift=float(np.max(abs(u.mean(axis=(1,2,3))-mean0)));assert drift<1e-10
   div=real(1j*sum(k*c for k,c in zip(f.k,h)));scaled=float(np.max(abs(div))*f.dx/rr['peak_speed']);assert scaled<1e-10
   maxima[role+'_mean_drift']=max(maxima.get(role+'_mean_drift',0),drift);maxima[role+'_scaled_divergence']=max(maxima.get(role+'_scaled_divergence',0),scaled)
   production=0.;gradloss=0.
   for ai in range(3):
    for bj in range(3):
     grad=real(1j*f.k[bj]*h[ai]);production+=float(np.sum(w[ai]*w[bj]*grad)*f.dx**3)
     if role=='reference':transfer-=float(np.sum(delta[ai]*delta[bj]*grad)*f.dx**3)
    gradw=real(1j*f.k[ai]*wh);gradloss+=float(np.sum(gradw*gradw)*f.dx**3)
   measured.update(enstrophy_production=production,viscous_energy_loss_rate=.02*measured['enstrophy'],viscous_enstrophy_loss_rate=.01*gradloss)
   for k,v in measured.items():check(role+'_'+k,v,rr[k])
   accum=[rr['I'],rr['integrated_viscous_energy_loss'],0.,rr['integrated_enstrophy_production'],rr['integrated_viscous_enstrophy_loss']]
   observed=f.observe(h,rr['t'],np.array(accum),first)
   for k in ['high_band_energy_fraction','high_band_enstrophy_fraction','energy_budget_relative_residual','enstrophy_budget_relative_residual']:check(role+'_'+k,observed[k],rr[k],absolute=True)
   for k in ['width_at_global_peak','width_at_central_roi_peak']:
    check(role+'_'+k,observed[k]['minimum_chord_cells'],rr[k]['minimum_chord_cells'],absolute=True)
   assert abs(rr['energy_budget_relative_residual'])<p['hard_stops']['energy_budget_relative_residual']
   measured.update(mean_drift=drift,scaled_divergence=scaled,distance_from_initial_l2=float(np.sqrt(f.inner(h-h0,h-h0)/f.inner(h0,h0))))
   parts[role]=measured
   if i in [0,10,20,30,40]:
    kmag=np.sqrt(f.k2);edges=np.arange(0,np.ceil(kmag.max())+2)
    ew=.5*np.sum((h.conj()*h).real,axis=0)*f.weights[0]*f.scale
    zw=.5*np.sum((wh.conj()*wh).real,axis=0)*f.weights[0]*f.scale
    spectra.append({'role':role,'t':row['t'],'k_edges':edges.tolist(),'energy':np.histogram(kmag,bins=edges,weights=ew)[0].tolist(),'enstrophy':np.histogram(kmag,bins=edges,weights=zw)[0].tolist()})
  dw=real(f.curl(hd-hb));loss=float(.01*np.sum(dw*dw)*f.dx**3);net=transfer-loss
  rd=f.rhs(hd)[0];rb=f.rhs(hb)[0];rhs_net=f.inner(hd-hb,rd-rb)
  scale=max(abs(net),abs(loss),abs(transfer),1e-20)
  identity=abs(net-rhs_net)/scale;assert identity<1e-8
  maxima['difference_transfer_identity']=max(maxima.get('difference_transfer_identity',0),identity)
  for k,v in [('difference_energy_transfer_from_reference',transfer),('difference_energy_viscous_loss_rate',loss),('difference_energy_net_rate',net)]:
   err=abs(v-row[k])/scale;maxima[k]=max(maxima.get(k,0),err);assert err<1e-8
  de=float(.5*np.sum(delta*delta)*f.dx**3);check('difference_energy',de,row['difference_energy'])
  assert np.max(abs(delta.mean(axis=(1,2,3))))<1e-10
  lowh=np.load(lowp/f'field-{i:03d}.npy',allow_pickle=False);lowref=np.load(lowb/f'field-{i:03d}.npy',allow_pickle=False)
  lowD=float(np.sqrt(np.sum(real(lowh-lowref)**2)*f.dx**3)/norm0)
  check('lower_viscosity_difference_norm',lowD,lr['series'][i]['D_relative_initial_baseline_l2'])
  field_hashes={'disturbed':sha(PILOT/f'fields/field-{i:03d}.npy'),'reference':sha(REF/f'field-{i:03d}.npy')}
  assert field_hashes['disturbed']==old['fields'][i]['disturbed_sha256']
  assert field_hashes['reference']==row['reference_field_sha256']==old['fields'][i]['reference_sha256']
  records.append({'t':row['t'],'field_hashes':field_hashes,'physical_measurements':parts,'D':physical_D,'amplification':physical_D/.001,'transfer':transfer,'viscous_loss':loss,'net_rate':net,'rhs_net_rate':rhs_net,'difference_energy':de,'lower_viscosity_D':lowD,'lower_viscosity_amplification':lowD/.001})
  print(f'Extended paired audit {i+1}/41 at {row["t"]:.2f}',flush=True)
 assert r['max_speed_cfl']<.5 and b['max_speed_cfl']<.5
 for k,v in p['pilot_source_hashes'].items():assert sha(PILOT/k)==v
 for k,v in p['original_source_hashes'].items():assert sha(STUDY/k)==v
 ts=np.array([z['t'] for z in records]);rates=np.array([z['net_rate'] for z in records])
 trapezoid=float(np.sum(.5*(rates[:-1]+rates[1:])*np.diff(ts)))
 dechange=records[-1]['difference_energy']-records[0]['difference_energy']
 report={'status':'passed','checked_utc':datetime.now(timezone.utc).isoformat(),'through':.4,'field_pairs':41,'fields_remeasured':82,'lower_viscosity_norm_pairs_checked':41,'new_simulations':0,'initial_reference_prescribed':True,'initial_disturbance_seed202_prescribed':True,'reference_initial_l2':norm0,'max_errors':maxima,'max_speed_cfl':{'disturbed':r['max_speed_cfl'],'reference':b['max_speed_cfl']},'result_sha256':sha(PILOT/'result.json'),'reference_result_sha256':sha(REF/'result.json'),'lower_viscosity_result_sha256':sha(lowp/'result.json'),'lower_viscosity_reference_result_sha256':sha(lowb/'result.json'),'protocol_sha256':sha(PILOT/'protocol.json'),'pilot_source_hashes':p['pilot_source_hashes'],'original_source_hashes':p['original_source_hashes'],'records':records,'spectra':spectra,'saved_output_difference_budget':{'quadrature':'trapezoid at 0.01 saved spacing; not a Heun-stage integral','difference_energy_change':dechange,'integrated_net_rate':trapezoid,'relative_gap':(dechange-trapezoid)/max(abs(dechange),1e-20)},'scope':'Independent NumPy physical-space norms, means, divergence, stretching and losses; frozen spectral/width routine remeasurement. Transfer checked against physical gradient and RHS difference. Existing Heun-stage individual budgets rechecked; no trajectories repeated.'}
 OUTPUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
 print('Extended verification passed',flush=True)
if __name__=='__main__':main()

