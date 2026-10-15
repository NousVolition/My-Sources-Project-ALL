"""Build reports from completed, independently checked saved results; no integration.
Usage: python -B build_completed_batch.py
Requires the original local study/pilot directories and verification/base-inputs beside this file.
"""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, shutil, os, ast, re, posixpath
from urllib.parse import urlsplit, unquote
import numpy as np
P=Path(__file__).resolve().parent
S=Path('C:/Users/ziman/Documents/Codex/2026-10-05/ca/work/matched-stretch-study')
H=P.parents[1]/'viscosity-disturbance0p1'
D=P.parents[1]/'viscosity-disturbance0p01'
Q=P.parent/'disturbance01-20261014'
ROOT=P/'publication'
PAIR='reports/matched-stretch/completed-n128-viscosity0p1-disturbance-seed202'
MAT='reports/matched-stretch/completed-n128-parameter-batch-20261015'
JOBS=['viscosity-n128-nu0.002','viscosity-n128-nu0.0005','viscosity-n128-nu0.0','radius-n128-v0.18','radius-n128-v0.22','spin_factor-n128-v0.9','spin_factor-n128-v1.1']
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def put(path,text):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8',newline='')
def js(path,data):put(path,json.dumps(data,indent=2)+'\n')
def cp(src,dst):dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
def complete(r):assert r['status']=='complete' and r['last_output']==40 and len(r['series'])==41
base=read(P/'base-inputs.json');protocol=read(S/'protocol.json')
for n,h in protocol['source_hashes'].items():assert sha(S/n)==h
os.environ.setdefault('MPLCONFIGDIR',str(P/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,out,name):fig.savefig(out/name,dpi=155);plt.close(fig)
def finish_axes(axes):
 for ax in np.asarray(axes).flat:ax.set_xlabel('Saved model time');ax.grid(alpha=.2);ax.legend(fontsize=8)

def pair_report():
 out=ROOT/PAIR;out.mkdir(parents=True,exist_ok=True)
 p=read(H/'protocol.json');r=read(H/'result.json');v=read(P/'verification-pair.json');lowv=read(Q/'verification.json')
 complete(r);assert v['status']=='passed' and v['field_pairs']==41 and v['result_sha256']==sha(H/'result.json')
 low=read(D/'result.json');lb=read(S/'runs/viscosity-n128-nu0.01/result.json');ll=read(S/'runs/perturb-n128-a0.001-s202/result.json');llb=read(S/'runs/baseline-n128-base/result.json')
 for z in [low,lb,ll,llb]:complete(z)
 assert sha(D/'result.json')==lowv['result_sha256'] and sha(S/'runs/viscosity-n128-nu0.01/result.json')==lowv['reference_result_sha256']
 for n,h in p['source_hashes'].items():assert sha(H/n)==h
 assert sha(H/'protocol.json')==v['protocol_sha256']
 for name in ['result.json','protocol.json','initial-audit.json','publication-authorization.json','saved-field-verification.json','summary.json','report.html','comparison.png']:
  cp(H/name,out/name)
 for n in ['run_pair.py','diagnostics.py','build_report.py','verify_saved.py','source/numerics.py']:cp(H/n,out/'reproduction/pilot'/n)
 for src,dst in [(P/'verification-pair.json','verification.json'),(P/'verify_pair.py','verify_pair.py'),(D/'result.json','viscosity0p01-result.json'),(S/'runs/viscosity-n128-nu0.01/result.json','viscosity0p01-reference-result.json'),(S/'runs/perturb-n128-a0.001-s202/result.json','viscosity0p001-result.json'),(S/'runs/baseline-n128-base/result.json','viscosity0p001-reference-result.json'),(Q/'verification.json','previous-viscosity0p01-verification.json'),(S/'run_suite.py','reproduction/original-run_suite.py'),(Path(__file__),'build_completed_batch.py')]:cp(src,out/dst)
 rows=r['series'];refs=[x['reference_diagnostics'] for x in rows];t=np.array([x['t'] for x in rows]);last=rows[-1];br=refs[-1]
 assert list(t)==[x['t'] for x in low['series']]==[x['t'] for x in ll['series']]
 amps={'.1':[x['amplification'] for x in rows],'.01':[x['amplification'] for x in low['series']],'.001':[x['lower_viscosity_amplification'] for x in v['records']]}
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for nu,a in amps.items():axes[0].semilogy(t,a,label='viscosity '+nu)
 axes[0].axhline(1,color='gray',lw=1);axes[0].set(title='Net change since the start',ylabel='D(t) / D(0), logarithmic scale')
 for key,label in [('difference_energy_transfer_from_reference','Transfer from reference'),('difference_energy_viscous_loss_rate','Viscous removal'),('difference_energy_net_rate','Net instantaneous rate')]:axes[1].plot(t,[x[key] for x in rows],label=label)
 axes[1].axhline(0,color='gray',lw=1);axes[1].set(title='Viscosity 0.1: signed instantaneous rate',ylabel='Difference energy per model time')
 finish_axes(axes);fig.suptitle('Identical seed and saved times; each disturbance uses its matching reference');save(fig,out,'disturbance-and-rates.png')
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for role,rr,color in [('reference',refs,'#176b9c'),('disturbed',rows,'#bb4d28')]:
  for key,style,label in [('energy_budget_relative_residual','-','E'),('enstrophy_budget_relative_residual','--','Z')]:axs[0,0].semilogy(t,np.maximum(np.abs([x[key] for x in rr]),1e-16),style,color=color,label=role+' '+label)
  for key,style,label in [('high_band_energy_fraction','-','E'),('high_band_enstrophy_fraction','--','Z')]:axs[0,1].semilogy(t,[x[key] for x in rr],style,color=color,label=role+' '+label)
  for key,style,label in [('width_at_global_peak','-','global'),('width_at_central_roi_peak','--','central')]:axs[1,0].plot(t,[x[key]['minimum_chord_cells'] for x in rr],style,color=color,label=role+' '+label)
  for key,style in [('energy','-'),('enstrophy','--')]:axs[1,1].plot(t,[x[key]/rr[0][key] for x in rr],style,color=color,label=role+' '+key)
 axs[1,0].axhline(6,color='gray',lw=1)
 for ax,title,ylabel in zip(axs.flat,['Individual Heun-stage budget residuals','High-band fractions','Minimum half-peak chord; screen at 6 cells','Physical totals relative to initial'],['Absolute relative residual','Fraction of total','Grid cells','Ratio to initial']):ax.set(title=title,ylabel=ylabel)
 finish_axes(axs);fig.suptitle('Viscosity 0.1: budget consistency does not establish spatial resolution');save(fig,out,'budgets-resolution.png')
 fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for z in v['spectra']:
  if z['t'] not in [0,.2,.4]:continue
  k=(np.array(z['k_edges'][:-1])+np.array(z['k_edges'][1:]))/2
  for ax,key in zip(axs,['energy','enstrophy']):
   y=np.array(z[key]);mask=y>1e-22;ax.semilogy(k[mask],y[mask],'-' if z['role']=='reference' else '--',color={0:'#258060',.2:'#226ab5',.4:'#bd542d'}[z['t']],label=f"{z['role']}, t={z['t']:.2f}")
 for ax,title in zip(axs,['Shell energy','Shell enstrophy']):ax.set(title=title,xlabel='Angular wavenumber magnitude',ylabel='Shell-integrated amount');ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.suptitle('Reference solid, disturbed dashed; unit-width angular-wavenumber bins');save(fig,out,'spectra.png')
 fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
 for role in ['reference','disturbed']:axs[0].plot(t,[x['physical_measurements'][role]['distance_from_initial_l2'] for x in v['records']],label=role)
 axs[0].set(title='Distance from each own initial field',ylabel='L2 distance / own initial L2 norm')
 de=np.array([x['difference_energy'] for x in rows]);rate=np.array([x['difference_energy_net_rate'] for x in rows]);coarse=np.r_[0,np.cumsum(.5*(rate[1:]+rate[:-1])*np.diff(t))]
 axs[1].plot(t,de-de[0],label='Measured difference-energy change');axs[1].plot(t,[x['integrated_difference_net'] for x in rows],'--',label='Every accepted-step trapezoidal integral');axs[1].plot(t,coarse,':',label='Coarse saved-output trapezoid')
 axs[1].set(title='Integrated difference-energy budget',ylabel='Model difference energy');finish_axes(axs);save(fig,out,'return-and-difference-budget.png')
 fails={role:sum(not x['resolved_screen'] for x in rr) for role,rr in [('reference',refs),('disturbed',rows)]}
 summary={'status':'complete_verified','through':.4,'paired_outputs':41,'independent_fields_checked':82,'amplification':last['amplification'],'percent_above_start':100*(last['amplification']-1),'instantaneous_net_rate':last['difference_energy_net_rate'],'transfer':last['difference_energy_transfer_from_reference'],'viscous_loss':last['difference_energy_viscous_loss_rate'],'difference_budget_residual_relative_to_initial_energy':last['difference_energy_budget_relative_residual'],'spatial_screen_failures':fails,'t0p20_comparison':{nu:a[20] for nu,a in amps.items()},'lower_viscosities_complete':True,'limitations':p['limitations']}
 js(out/'analysis.json',summary)
 table='\n'.join(f"| {t[i]:.2f} | {amps['.1'][i]:.9g} | {amps['.01'][i]:.9g} | {amps['.001'][i]:.9g} | {rows[i]['difference_energy_net_rate']:+.9g} |" for i in [0,1,2,10,20,30,40])
 text=f'''# Completed viscosity 0.1 disturbance comparison — seed 202

[Matched-stretch reports](../README.md) · [Completed viscosity 0.01 comparison](../completed-n128-viscosity0p01-disturbance-seed202/README.md)

The disturbance reached **{last['amplification']:.9f} times its starting L2 size** at saved model time **0.40**, or **{100*(last['amplification']-1):.6f}% above its starting size**. Its instantaneous difference-energy rate was **positive ({last['difference_energy_net_rate']:+.9f})** at that state. Both trajectories are complete; all **41 reference and 41 disturbed saved fields** passed independent verification.

At the requested matching time **0.20**, the disturbance ratios are **{amps['.1'][20]:.9f}** at viscosity 0.1, **{amps['.01'][20]:.9f}** at 0.01, and **{amps['.001'][20]:.9f}** at 0.001. The viscosity-0.1 rate at 0.20 is **{rows[20]['difference_energy_net_rate']:+.9f}**. It has less separation than the lower-viscosity cases, while being larger than its own starting disturbance and growing instantaneously. Both lower-viscosity comparisons are complete through 0.40.

## Matched experiment and normalization

Grid 128, domain side 6, the supplied reference initial field, fixed mean, zero force, frozen Heun method, identical 0.1% seed-202 disturbance, one FFT worker and time step **{p['dt']:.17g}** were preserved. There are 205 steps per 0.01 output interval, ending at 0.40. This authorized pair contains two trajectories at viscosity 0.1; it does not reuse a reference with a different viscosity.

D(t) = ||u_disturbed(t) − u_reference(t)||₂ / ||u_reference(0)||₂. The fixed initial reference norm is **{v['reference_initial_l2']:.14g}**, D(0)=0.001 and D(0.40)={last['D_relative_initial_reference_l2']:.12g}. Each lower-viscosity comparison uses its own reference at the identical saved time and the same initial normalization. No missing times were interpolated.

| Saved time | D(t)/D(0), ν=0.1 | ν=0.01 | ν=0.001 | Instantaneous difference-energy rate, ν=0.1 |
|---|---:|---:|---:|---:|
{table}

At t=0.01 the disturbance remains below its starting size while its instantaneous rate is already positive. By t=0.02 it is above its starting size. These saved observations distinguish net change from current growth; they do not locate an exact reversal between outputs. The earlier one-step preflight was not evidence of sustained shrinking.

![Disturbance ratios and instantaneous rates](disturbance-and-rates.png)

## Transfer and budget checks

For δu=u_disturbed−U and Eδ=½||δu||₂², the zero-force periodic identity is dEδ/dt = −∫δuᵢδuⱼ∂ⱼUᵢ dx − ν||∇δu||₂². At 0.40, transfer is **{last['difference_energy_transfer_from_reference']:.9f}**, viscous removal is **{last['difference_energy_viscous_loss_rate']:.9f}**, and the net rate is **{last['difference_energy_net_rate']:+.9f}**. Independent physical gradients agree with the spectral RHS difference; maximum scaled identity discrepancy is {v['max_errors']['difference_transfer_identity']:.4g}.

The runner accumulated the difference-energy rate by **trapezoidal integration at every accepted time step**, using rates at the old and accepted new states. This is not a Heun predictor-stage integral. The final measured difference-energy change is {de[-1]-de[0]:.12g}, and the recorded integral is {last['integrated_difference_net']:.12g}. Their residual divided by initial difference energy is **{last['difference_energy_budget_relative_residual']:.9g}**. Saved-field verification rechecks the identity, endpoint energies, recorded integral and residual; it does not reconstruct unsaved intermediate states. Coarse quadrature at 0.01 saved spacing is shown separately and is not a time-integration error estimate.

Individual energy/enstrophy budgets retain their recorded Heun-stage integrals. Final relative energy residuals are {br['energy_budget_relative_residual']:.6g} / {last['energy_budget_relative_residual']:.6g}, and enstrophy residuals {br['enstrophy_budget_relative_residual']:.6g} / {last['enstrophy_budget_relative_residual']:.6g} (reference / disturbed). Small residuals do not prove convergence.

![Budget, spectral-fraction and width checks](budgets-resolution.png)

## Resolution and interpretation limits

**Raw strain is nonsmooth across periodic joins, and initial maxima lie outside the central tube.** Reference / disturbed fields fail **{fails['reference']} / {fails['disturbed']} of 41** spatial screens. At 0.40, minimum global half-peak chords are **{br['width_at_global_peak']['minimum_chord_cells']:.6g} / {last['width_at_global_peak']['minimum_chord_cells']:.6g} cells**, below the six-cell screen; central widths are {br['width_at_central_roi_peak']['minimum_chord_cells']:.6g} / {last['width_at_central_roi_peak']['minimum_chord_cells']:.6g} cells. High-band enstrophy fractions are {br['high_band_enstrophy_fraction']:.6%} / {last['high_band_enstrophy_fraction']:.6%}.

![Saved shell spectra](spectra.png)

This is one finite-amplitude seed, one grid and one time step. No paired finer-grid or half-step test was authorized. Spatial/time convergence is unestablished; this does not establish a physical stability threshold, infinitesimal stability or an asymptotic Lyapunov exponent. Peak vorticity is a separate diagnostic, not disturbance growth. Initial-return distances below do not establish recurrence or periodicity and do not search all state pairs.

![Initial-return and integrated difference-energy diagnostics](return-and-difference-budget.png)

## Verification and data

[Extended verification](verification.json) checks all 82 paired fields, prescribed initial states, exact seed normalization, means, divergence, physical norms, stretching, viscous losses, transfer identity, endpoint budgets, frozen hashes, spectra and widths. All 41 matching viscosity-0.001 norms were remeasured. Viscosity-0.01 comparisons reuse its [previous completed verification](previous-viscosity0p01-verification.json), with unchanged result hashes. The [existing verifier](saved-field-verification.json) also independently remeasured all 82 fields. No trajectory was repeated.

Spectral fractions and chord widths are remeasured using the frozen observation routine: consistency checks, not a second discretization. [Verification code](verify_pair.py) runs as `python -B verify_pair.py PILOT_DIRECTORY OUTPUT_JSON` in the original study environment; original large saved field arrays are required and are not uploaded. Per-field SHA-256 values are retained. [Publication builder](build_completed_batch.py) requires the original local paths and adjacent verification inputs.

- [Paired raw diagnostics, including reference rows](result.json), [analysis summary](analysis.json), [existing report summary](summary.json), [standalone report](report.html)
- [Viscosity 0.01 disturbed](viscosity0p01-result.json) and [reference](viscosity0p01-reference-result.json); [viscosity 0.001 disturbed](viscosity0p001-result.json) and [reference](viscosity0p001-reference-result.json)
- [Frozen protocol](protocol.json), [initial audit](initial-audit.json), [later publication authorization](publication-authorization.json)
- [Frozen runner](reproduction/pilot/run_pair.py), [diagnostics](reproduction/pilot/diagnostics.py), [report builder](reproduction/pilot/build_report.py), [original verifier](reproduction/pilot/verify_saved.py), [numerics](reproduction/pilot/source/numerics.py), [original matrix controller](reproduction/original-run_suite.py)

The frozen protocol's historical local-only restriction remains byte-for-byte intact. The later explicit authorization supersedes only publication scope. Historical intake reviews remain unchanged; the original matrix continues toward its unchanged endpoint and final report.
'''
 put(out/'README.md',text);return summary

def matrix_report():
 out=ROOT/MAT;out.mkdir(parents=True,exist_ok=True);runs={};ver={}
 for job in JOBS:
  r=read(S/'runs'/job/'result.json');v=read(P/f'verification-{job}.json');complete(r)
  assert v['status']=='passed' and v['fields_checked']==41 and v['result_sha256']==sha(S/'runs'/job/'result.json')
  runs[job]=r;ver[job]=v;cp(S/'runs'/job/'result.json',out/'diagnostics'/f'{job}.json');cp(P/f'verification-{job}.json',out/'verification'/f'{job}.json')
 controls=['baseline-n128-base','baseline-n128-half','baseline-n64-base','baseline-n64-half','viscosity-n128-nu0.01']+[j.replace('n128','n64') for j in JOBS]
 for job in controls:
  r=read(S/'runs'/job/'result.json');complete(r);runs[job]=r;cp(S/'runs'/job/'result.json',out/'diagnostics'/f'{job}.json')
 for src,dst in [(S/'protocol.json','protocol.json'),(S/'numerics.py','reproduction/numerics.py'),(S/'run_suite.py','reproduction/run_suite.py'),(P/'verify_matrix.py','verify_matrix.py'),(Path(__file__),'build_completed_batch.py')]:cp(src,out/dst)
 t=np.array([x['t'] for x in runs[JOBS[0]]['series']]);assert all(list(t)==[x['t'] for x in r['series']] for r in runs.values())
 labels={j:j.replace('-n128','').replace('viscosity-nu','nu=').replace('-v','=') for j in JOBS};colors=dict(zip(JOBS,plt.cm.tab10(np.arange(7))))
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for j in JOBS+['baseline-n128-base']:
  rr=runs[j]['series'];kw={'label':labels.get(j,'baseline nu=0.001'),'color':colors.get(j,'black'),'lw':1.4}
  for ax,key in zip(axs.flat,['Wmax','energy','enstrophy','I']):ax.plot(t,[x[key]/rr[0][key] if key!='I' else x[key] for x in rr],**kw)
 for ax,title,yl in zip(axs.flat,['Peak vorticity','Kinetic energy','Enstrophy','Integrated peak vorticity'],['W/W(0)','E/E(0)','Z/Z(0)','I(t)']):ax.set(title=title,ylabel=yl)
 finish_axes(axs);fig.suptitle('Completed one-factor sweeps at grid 128; each curve uses its own initial state');save(fig,out,'parameter-curves.png')
 fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained');ix=np.arange(7)
 for ax,key,title in zip(axs,['Wmax','energy','enstrophy'],['Peak vorticity','Energy','Enstrophy']):
  for n,offset,col in [(64,-.18,'#176b9c'),(128,.18,'#bb4d28')]:
   vals=[]
   for j in JOBS:
    rr=runs[j.replace('n128',f'n{n}')]['series'];vals.append(rr[-1][key]/rr[0][key])
   ax.bar(ix+offset,vals,.36,label=f'grid {n}',color=col)
  ax.set_xticks(ix,[labels[j] for j in JOBS],rotation=65,ha='right');ax.set(title=title+' at t=0.40',ylabel='Final / own initial');ax.grid(axis='y',alpha=.2);ax.legend()
 fig.suptitle('Existing matched 64/128 controls; disagreement is not a convergence proof');save(fig,out,'grid-endpoints.png')
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for j in JOBS:
  rr=runs[j]['series'];label=labels[j];col=colors[j]
  for ax,key in zip(axs[0],['energy_budget_relative_residual','enstrophy_budget_relative_residual']):ax.semilogy(t,np.maximum(np.abs([x[key] for x in rr]),1e-16),label=label,color=col)
  axs[1,0].semilogy(t,[x['high_band_enstrophy_fraction'] for x in rr],label=label,color=col)
  axs[1,1].plot(t,[x['width_at_global_peak']['minimum_chord_cells'] for x in rr],label=label,color=col)
 axs[1,0].axhline(.01,color='gray',lw=1);axs[1,1].axhline(6,color='gray',lw=1)
 for ax,title,yl in zip(axs.flat,['Energy budget','Enstrophy budget','High-band enstrophy; screen at 1%','Global minimum chord; screen at 6 cells'],['Absolute relative residual','Absolute relative residual','Fraction','Grid cells']):ax.set(title=title,ylabel=yl)
 finish_axes(axs);save(fig,out,'budgets-resolution.png')
 fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
 for j in JOBS:
  sp=ver[j]['spectra'][-1];k=(np.array(sp['k_edges'][:-1])+np.array(sp['k_edges'][1:]))/2
  for ax,key in zip(axs,['energy','enstrophy']):
   y=np.array(sp[key]);mask=y>1e-22;ax.semilogy(k[mask],y[mask],label=labels[j],color=colors[j])
 for ax,key in zip(axs,['energy','enstrophy']):ax.set(title='Final shell '+key,xlabel='Angular wavenumber magnitude',ylabel='Shell-integrated amount');ax.grid(alpha=.2);ax.legend(fontsize=8)
 save(fig,out,'spectra.png')
 fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
 for j in JOBS:ax.plot(t,[x['distance_from_initial_l2'] for x in ver[j]['fields']],label=labels[j],color=colors[j])
 ax.set(title='Distance from each trajectory’s own initial field',ylabel='L2 distance / own initial L2 norm');finish_axes([ax]);save(fig,out,'initial-return.png')
 summary=[]
 for j in JOBS:
  rr=runs[j]['series'];a=rr[-1];summary.append({'job':j,'status':'complete_verified','saved_outputs':41,'dt':runs[j]['dt'],'nu':runs[j]['job']['nu'],'W_ratio':a['Wmax']/rr[0]['Wmax'],'energy_ratio':a['energy']/rr[0]['energy'],'enstrophy_ratio':a['enstrophy']/rr[0]['enstrophy'],'I':a['I'],'spatial_screen_failures':sum(not x['resolved_screen'] for x in rr),'final_global_chord_cells':a['width_at_global_peak']['minimum_chord_cells'],'final_central_chord_cells':a['width_at_central_roi_peak']['minimum_chord_cells'],'final_high_band_enstrophy_fraction':a['high_band_enstrophy_fraction'],'max_speed_cfl':runs[j]['max_speed_cfl'],'max_absolute_energy_budget_residual':max(abs(x['energy_budget_relative_residual']) for x in rr),'max_absolute_enstrophy_budget_residual':max(abs(x['enstrophy_budget_relative_residual']) for x in rr)})
 state=read(S/'task-state.json');js(out/'analysis.json',{'status':'complete_verified_batch','new_fields_checked':287,'new_simulations':0,'runs':summary,'original_matrix_snapshot':state,'comparisons':'Existing completed matching saved times only; controls reused, not rerun or newly independently reverified.'})
 table='\n'.join(f"| [{x['job']}](diagnostics/{x['job']}.json) | {x['dt']:.7g} | {x['W_ratio']:.7g} | {x['energy_ratio']:.7g} | {x['enstrophy_ratio']:.7g} | {x['spatial_screen_failures']}/41 |" for x in summary)
 limits='\n'.join(f"| {labels[x['job']]} | {x['final_global_chord_cells']:.5g} | {x['final_central_chord_cells']:.5g} | {x['final_high_band_enstrophy_fraction']:.5%} | {x['max_absolute_energy_budget_residual']:.3g} | {x['max_absolute_enstrophy_budget_residual']:.3g} |" for x in summary)
 links='\n'.join(f"- [{j} independent verification](verification/{j}.json)" for j in JOBS)
 text=f'''# Completed grid-128 parameter batch — October 15, 2026

[Matched-stretch reports](../README.md) · [Completed viscosity 0.01 sweep](../completed-n128-viscosity0p01/README.md)

Seven additional original-matrix runs reached **0.40** with **41 saved fields each**. All **287 fields** were independently remeasured. The frozen supplied initial-field family, per-run viscosity, mean, domain, Heun method and zero force were preserved. No new simulations or horizon extensions were introduced.

The batch contains viscosities 0.002, 0.0005 and 0; radius 0.18 and 0.22; and spin factors 0.9 and 1.1. Radius/spin cases retain viscosity 0.001. These are the already authorized one-factor changes to the supplied field, not perturbation-separation experiments. Their means and normalization are checked against each run's prescribed initial state. Each curve is divided by its own initial value; absolute quantities remain in raw diagnostics.

| Run | Preserved dt | W(0.40)/W(0) | E(0.40)/E(0) | Z(0.40)/Z(0) | Failed spatial screens |
|---|---:|---:|---:|---:|---:|
{table}

![Parameter curves](parameter-curves.png)

## Existing controls and numerical limits

The matching completed grid-64 parameter controls and completed grid-128 baseline (viscosity 0.001) are included at identical saved times. Completed baseline half-step records and the previous viscosity-0.01 sweep are also supplied as context. Each retains its own recorded viscosity, initial field and time step. The 64/128 endpoint comparison below does not establish spatial convergence; baseline half-step checks do not establish time convergence of every parameter trajectory.

![Matching grid controls at the endpoint](grid-endpoints.png)

**Raw strain is nonsmooth across periodic joins, and the initial maxima lie outside the central tube.** The published spatial screens require high-band energy below 0.1%, high-band enstrophy below 1%, and minimum half-peak chord at least six cells. These are screening criteria, not proof of resolution. The table distinguishes spatial limitations from budget consistency.

| Run | Final global width, cells | Final central width, cells | Final high-band Z | Maximum absolute E budget residual | Maximum absolute Z budget residual |
|---|---:|---:|---:|---:|---:|
{limits}

![Budgets and resolution checks](budgets-resolution.png)

Individual budgets use recorded Heun-stage integrals, checked against saved endpoint energies and enstrophies. Independent physical-space checks cover norms, means, divergence, enstrophy production and viscous losses. All declared numerical hard stops pass; neither small residuals nor a viscosity change proves a physical stability transition. At zero viscosity the viscous terms are exactly zero; unresolved transfer to small scales remains a limitation.

![Final shell spectra](spectra.png)

The initial-return diagnostic measures distance from each trajectory's own initial velocity field. It does not establish recurrence, periodicity or a search over all state pairs. These unperturbed sweeps do not supply disturbance-growth or Lyapunov measurements.

![Initial-return diagnostics](initial-return.png)

## Verification and reproducibility

The [verifier](verify_matrix.py) runs as `python -B verify_matrix.py STUDY_DIRECTORY JOB_ID OUTPUT_JSON`. It checks all 41 saved fields per new run, the exact prescribed initial field, unchanged source/protocol hashes, mean, divergence, physical norms, stretching and losses, saved budget endpoints, spectra, widths and declared hard stops. Width and spectral remeasurement uses the frozen routine and is a consistency check, not a second discretization. Original large Fourier fields are required locally and are not embedded; their SHA-256 values are retained below.

{links}

[Raw diagnostics](diagnostics/) include the seven new runs and existing completed controls. Existing controls are reused without rerunning or claiming a new independent full-field verification. [Analysis](analysis.json) records summary values and the current completion snapshot. [Frozen protocol](protocol.json), [numerical source](reproduction/numerics.py), [matrix controller](reproduction/run_suite.py) and [publication builder](build_completed_batch.py) retain provenance.

At this report's snapshot, the original matrix has **{len(state['completed'])}/48 complete**, **{len(state['active'])} active**, **{len(state['pending'])} queued**, and **{len(state['failures'])} recorded failures**. Active jobs: {', '.join(state['active'])}. This is a completed batch report; the original matrix's final report awaits its remaining jobs. Historical intake reviews remain unchanged.
'''
 put(out/'README.md',text);return summary

pair=pair_report();matrix=matrix_report()
for path,text in base['files'].items():
 notes=[]
 for rel,label,detail in [(PAIR,'Viscosity 0.1 disturbed/reference pair',f"41 independently verified paired outputs through 0.40; D/D(0)={pair['amplification']:.6f}, instantaneous rate positive"),(MAT,'Seven completed grid-128 parameter runs','287 saved fields independently remeasured; original matrix continues')]:
  link=posixpath.relpath(rel+'/README.md',posixpath.dirname(path) or '.');assert rel not in text
  notes.append(f"**Completed October 15, 2026:** [{label}]({link}) — {detail}. Raw-strain, initial-peak and resolution limits remain.\n\n")
 newline='\r\n' if '\r\n' in text else '\n'
 head,sep,tail=text.partition(newline*2);assert sep;put(ROOT/path,head+sep+''.join(notes).replace('\n',newline)+tail)
files=[]
for f in ROOT.rglob('*'):
 if not f.is_file():continue
 data=f.read_bytes();path=f.relative_to(ROOT).as_posix()
 if f.suffix=='.py':ast.parse(data)
 if f.suffix=='.json':json.loads(data)
 files.append({'path':path,'local':str(f.resolve()),'size':len(data),'sha256':hashlib.sha256(data).hexdigest(),'blob_sha':hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()})
paths={x['path'] for x in base['base']['blobs']}|{x['path'] for x in files};checks=0
for rel in [PAIR,MAT]:
 for url in re.findall(r'\]\(([^)]+)\)',(ROOT/rel/'README.md').read_text(encoding='utf-8')):
  if urlsplit(url).scheme or url.startswith('#'):continue
  target=posixpath.normpath(posixpath.join(rel,unquote(url.split('#')[0])))
  assert target in paths or any(x.startswith(target+'/') for x in paths),(url,target)
  checks+=1
js(P/'payload.json',files);js(P/'quality-verification.json',{'status':'passed','python_json_parse':'passed','links_checked':checks,'payload_files':len(files),'paired_fields_checked':82,'matrix_fields_checked':287,'new_simulations':0})
print(json.dumps({'payload_files':len(files),'bytes':sum(x['size'] for x in files),'links_checked':checks,'pair':pair,'matrix':matrix},indent=2))
