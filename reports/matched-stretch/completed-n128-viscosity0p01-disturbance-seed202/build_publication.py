"""Build the completed publication from immutable saved diagnostics; no integration.
Usage: python -B build_publication.py PILOT_DIRECTORY
Requires verification.json, base-inputs.json and verify_complete.py beside this script.
"""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,shutil,os,hashlib,ast,re,posixpath
from urllib.parse import urlsplit,unquote
import numpy as np
P=Path(__file__).resolve().parent
D=Path(sys.argv[1]).resolve()
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
p=read(D/'protocol.json');S=Path(p['original_study']);B=Path(p['reference_directory'])
v=read(P/'verification.json');assert v['status']=='passed' and v['field_pairs']==41
r=read(D/'result.json');b=read(B/'result.json')
low=read(S/'runs/perturb-n128-a0.001-s202/result.json')
lb=read(S/'runs/baseline-n128-base/result.json')
assert all(x['status']=='complete' and len(x['series'])==41 for x in [r,b,low,lb])
assert hashlib.sha256((D/'result.json').read_bytes()).hexdigest()==v['result_sha256']
assert hashlib.sha256((B/'result.json').read_bytes()).hexdigest()==v['reference_result_sha256']
base=read(P/'base-inputs.json')
rel='reports/matched-stretch/completed-n128-viscosity0p01-disturbance-seed202'
root=P/'publication';out=root/rel;out.mkdir(parents=True,exist_ok=True)
def copy(src,dst):
 dst=out/dst;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
for name in ['result.json','protocol.json','initial-audit.json','publication-authorization.json','saved-field-verification.json','summary.json','report.html','comparison.png']:
 copy(D/name,name)
for src,dst in [(B/'result.json','reference-result.json'),(S/'runs/perturb-n128-a0.001-s202/result.json','lower-viscosity-result.json'),(S/'runs/baseline-n128-base/result.json','lower-viscosity-reference-result.json'),(P/'verification.json','verification.json'),(P/'verify_complete.py','verify_complete.py'),(Path(__file__),'build_publication.py')]:copy(src,dst)
for name in ['run_pair.py','diagnostics.py','build_report.py','verify_saved.py','source/numerics.py']:
 copy(D/name,'reproduction/pilot/'+name)
copy(S/'run_suite.py','reproduction/original-run_suite.py')
os.environ.setdefault('MPLCONFIGDIR',str(P/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
t=np.array([x['t'] for x in r['series']]);colors={'disturbed':'#bb4d28','reference':'#176b9c'}
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for role,result in [('disturbed',r),('reference',b)]:
 rows=result['series'];color=colors[role]
 for key,style,label in [('energy_budget_relative_residual','-','Energy'),('enstrophy_budget_relative_residual','--','Enstrophy')]:
  axs[0,0].semilogy(t,np.maximum(np.abs([x[key] for x in rows]),1e-16),style,color=color,label=role+' '+label)
 for key,style,label in [('high_band_energy_fraction','-','Energy'),('high_band_enstrophy_fraction','--','Enstrophy')]:
  axs[0,1].semilogy(t,[x[key] for x in rows],style,color=color,label=role+' '+label)
 for key,style,label in [('width_at_global_peak','-','Global peak'),('width_at_central_roi_peak','--','Central peak')]:
  axs[1,0].plot(t,[x[key]['minimum_chord_cells'] for x in rows],style,color=color,label=role+' '+label)
 axs[1,1].plot(t,[x['energy']/rows[0]['energy'] for x in rows],color=color,label=role+' E/E(0)')
 axs[1,1].plot(t,[x['enstrophy']/rows[0]['enstrophy'] for x in rows],'--',color=color,label=role+' Z/Z(0)')
axs[0,0].set(title='Individual Heun-stage budget residuals',ylabel='Absolute relative residual')
axs[0,1].set(title='High-band spectral fractions',ylabel='Fraction of total')
axs[1,0].set(title='Minimum measured half-peak chord',ylabel='Grid cells')
axs[1,1].set(title='Energy and enstrophy evolve differently',ylabel='Ratio to initial value')
for ax in axs.flat:ax.set(xlabel='Saved model time');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Viscosity 0.01: small budget residuals do not establish spatial resolution')
fig.savefig(out/'budgets-resolution.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
timecolors={0:'#258060',.2:'#226ab5',.4:'#bd542d'}
for s in v['spectra']:
 if s['t'] not in timecolors:continue
 k=(np.array(s['k_edges'][:-1])+np.array(s['k_edges'][1:]))/2
 for ax,key in zip(axs,['energy','enstrophy']):
  y=np.array(s[key]);mask=y>1e-22
  ax.semilogy(k[mask],y[mask],'-' if s['role']=='reference' else '--',color=timecolors[s['t']],label=f"{s['role']}, t={s['t']:.2f}")
for ax,title in zip(axs,['Shell energy','Shell enstrophy']):
 ax.set(title=title,xlabel='Angular wavenumber magnitude',ylabel='Shell-integrated amount');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Saved spectra: reference solid, disturbed dashed; unit-width angular-wavenumber bins')
fig.savefig(out/'spectra.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
for role in ['reference','disturbed']:
 axs[0].plot(t,[x['physical_measurements'][role]['distance_from_initial_l2'] for x in v['records']],label=role,color=colors[role])
axs[0].set(title='Distance from each trajectory\'s own initial field',ylabel='L2 distance / own initial L2 norm')
energy=np.array([x['difference_energy'] for x in v['records']]);rates=np.array([x['net_rate'] for x in v['records']])
integral=np.r_[0,np.cumsum(.5*(rates[1:]+rates[:-1])*np.diff(t))]
axs[1].plot(t,energy-energy[0],label='Measured difference-energy change')
axs[1].plot(t,integral,'--',label='Saved-rate trapezoidal integral')
axs[1].set(title='Difference-energy budget at saved-output spacing',ylabel='Model difference energy')
for ax in axs:ax.set(xlabel='Saved model time');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Initial-return diagnostic and coarse budget quadrature; no recurrence or stage-accuracy claim')
fig.savefig(out/'return-and-difference-budget.png',dpi=160);plt.close(fig)
last=r['series'][-1];br=b['series'][-1];vl=v['records'][-1]
fails={role:sum(not row['resolved_screen'] for row in rr['series']) for role,rr in [('disturbed',r),('reference',b)]}
summary={'status':'complete_verified','verified_at':v['checked_utc'],'time':.4,'paired_saved_outputs':41,'amplification':last['amplification'],'percent_above_start':(last['amplification']-1)*100,'instantaneous_difference_energy_net_rate':last['difference_energy_net_rate'],'lower_viscosity_amplification':vl['lower_viscosity_amplification'],'spatial_screen_failures':fails,'saved_output_difference_budget':v['saved_output_difference_budget'],'scope':'One 0.1% seed 202, grid 128, fixed Heun step. No asymptotic or physical stability threshold claim.','limits':p['limitations']}
(out/'analysis.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
gap=v['saved_output_difference_budget']['relative_gap']
table='\n'.join(f"| {x['t']:.2f} | {x['amplification']:.8g} | {x['lower_viscosity_amplification']:.8g} | {x['net_rate']:+.8g} |" for x in v['records'] if round(x['t']*100) in [0,1,10,20,30,40])
readme=f"""# Completed viscosity 0.01 disturbance comparison — seed 202

[Matched-stretch reports](../README.md) · [Completed undisturbed reference](../completed-n128-viscosity0p01/README.md)

The 0.1% disturbance grew to **{last['amplification']:.8f} times its starting L2 size** by saved model time **0.40**, or **{(last['amplification']-1)*100:.4f}% above its starting size**. Its instantaneous difference-energy rate at 0.40 was **positive ({last['difference_energy_net_rate']:+.8f})**: it was still growing at that saved state. Both trajectories are complete, with all 41 paired outputs independently checked.

The matching completed viscosity 0.001 case reached **{vl['lower_viscosity_amplification']:.8f} times its starting disturbance** at the same time. This comparison records less finite-time separation at viscosity 0.01 for this seed; it does not establish a physical stability transition, infinitesimal stability, or an asymptotic Lyapunov exponent.

## What was held fixed

Grid 128, domain side 6, the supplied initial reference field, fixed mean velocity, zero force, frozen Heun method, one FFT worker and time step {p['dt']:.17g} were preserved. There are 205 steps per 0.01 output interval and 41 outputs from 0 to 0.40. The disturbance is the identical 0.1%, seed-202 initial perturbation used in the lower-viscosity comparison. The completed viscosity-0.01 reference was reused; this pilot added only one disturbed trajectory.

D(t) = ||u_disturbed(t) − u_reference(t)||₂ / ||u_reference(0)||₂, with the fixed initial reference norm **{v['reference_initial_l2']:.14g}**. D(0) = 0.001; D(0.40) = **{last['D_relative_initial_reference_l2']:.12g}**. Each lower-viscosity difference uses its own matching reference, with the same initial normalization. All comparisons below use identical saved times; both lower-viscosity trajectories are complete.

| Saved time | D(t)/D(0), ν=0.01 | D(t)/D(0), ν=0.001 | Instantaneous difference-energy rate, ν=0.01 |
|---|---:|---:|---:|
{table}

The instantaneous rate is negative at t=0 and positive at the first evolved saved output, t=0.01. No reversal time is interpolated between those outputs. Change since the start and the instantaneous sign are different measurements.

![Disturbance, transfer and peak-vorticity comparisons](comparison.png)

## Difference-energy transfer and budgets

For δu = u_disturbed − u_reference and Eδ = ½||δu||₂², the zero-force periodic identity is dEδ/dt = −∫δuᵢδuⱼ ∂ⱼUᵢ dx − ν||∇δu||₂². At 0.40, transfer from the reference flow is **{last['difference_energy_transfer_from_reference']:.8f}**, viscous removal is **{last['difference_energy_viscous_loss_rate']:.8f}**, and their difference is **{last['difference_energy_net_rate']:+.8f}**. Independent physical-gradient evaluation agrees with the spectral RHS difference; the maximum scaled identity error is {v['max_errors']['difference_transfer_identity']:.3g}.

The individual reference and disturbed energy/enstrophy budgets use their recorded Heun-stage integrals. Final absolute relative residuals are {abs(br['energy_budget_relative_residual']):.6g} / {abs(last['energy_budget_relative_residual']):.6g} for energy and {abs(br['enstrophy_budget_relative_residual']):.6g} / {abs(last['enstrophy_budget_relative_residual']):.6g} for enstrophy (reference / disturbed).

A separate trapezoidal integral of the saved difference-energy rates at spacing 0.01 gives {v['saved_output_difference_budget']['integrated_net_rate']:.9g}, versus the measured difference-energy change {v['saved_output_difference_budget']['difference_energy_change']:.9g}; their signed relative gap is {gap:.6g}. **This is coarse saved-output quadrature, not a Heun-stage difference-energy integral.** The pilot did not retain stage-by-stage difference-energy rates, so this gap is not a time-integration error estimate.

![Budgets, spectral fractions, widths and physical totals](budgets-resolution.png)

## Resolution and interpretation limits

The raw strain is **not smooth across periodic joins**, and **the initial maxima lie outside the central tube**. The reference and disturbed trajectories each fail {fails['reference']} and {fails['disturbed']} of 41 spatial screens, respectively. At 0.40 their high-band enstrophy fractions are {br['high_band_enstrophy_fraction']:.6%} and {last['high_band_enstrophy_fraction']:.6%}; global minimum half-peak widths are {br['width_at_global_peak']['minimum_chord_cells']:.6g} and {last['width_at_global_peak']['minimum_chord_cells']:.6g} cells. Central-region widths are {br['width_at_central_roi_peak']['minimum_chord_cells']:.6g} and {last['width_at_central_roi_peak']['minimum_chord_cells']:.6g} cells. These observations do not demonstrate resolved spatial dynamics.

![Saved energy and enstrophy spectra](spectra.png)

This is one finite-amplitude seed at one grid and time step. No paired finer-grid or half-step control was authorized, so time and spatial convergence of this separation are unestablished. Passing finite-field, mean, divergence, CFL and energy-budget hard stops, and small budget residuals, do not cure those limitations. Peak vorticity is a separate physical diagnostic and is not used as a substitute for disturbance growth.

The initial-return diagnostic below measures distance from each trajectory's own initial velocity field. It does not establish recurrence or periodicity, and is not a search over arbitrary state pairs.

![Initial-return distances and coarse difference-energy budget](return-and-difference-budget.png)

## Verification and reproducibility

[Extended verification](verification.json) independently remeasured all **41 reference and 41 disturbed fields**, with all 41 lower-viscosity difference norms also checked. [Verification code](verify_complete.py) checks the prescribed reference, exact seed-202 perturbation, fixed normalization and mean, physical L2/energy/enstrophy/peak-vorticity norms, divergence, independent stretching and viscous terms, the transfer identity, individual budget endpoints, frozen hashes, spectra and widths. The spectral/width remeasurement uses the frozen routine; it is a consistency check, not a second discretization. The [existing verifier's record](saved-field-verification.json) supplies another physical-space norm check and hashes of all paired fields. No trajectories were repeated during verification.

The extended verifier runs as `python -B verify_complete.py PILOT_DIRECTORY OUTPUT_JSON` in the study environment and requires the original saved Fourier fields and original reference paths recorded in the protocol. Large field arrays are not embedded in this report. Per-field SHA-256 values are retained in the verification records. The copied code documents the original environment and paths; copying a report folder alone is not sufficient to reproduce the simulation.

- [Disturbed raw diagnostics](result.json), [reference diagnostics](reference-result.json), [lower-viscosity disturbed diagnostics](lower-viscosity-result.json), [lower-viscosity reference diagnostics](lower-viscosity-reference-result.json)
- [Frozen protocol](protocol.json), [initial audit](initial-audit.json), [later publication authorization](publication-authorization.json)
- [Analysis summary](analysis.json), [existing report summary](summary.json), [standalone HTML report](report.html), [publication builder](build_publication.py)
- [Frozen pilot runner](reproduction/pilot/run_pair.py), [diagnostics](reproduction/pilot/diagnostics.py), [report builder](reproduction/pilot/build_report.py), [existing verifier](reproduction/pilot/verify_saved.py), [numerical source](reproduction/pilot/source/numerics.py), [original matrix controller](reproduction/original-run_suite.py)

The frozen protocol's earlier local-only publication wording is preserved byte-for-byte. The later explicit publication authorization supersedes only that restriction, without altering any numerical setting or source. Historical intake reviews remain unchanged. The original 48-job matrix and separate viscosity-0.1 comparison continue under their existing scope.
"""
(out/'README.md').write_text(readme,encoding='utf-8')
link=rel+'/README.md'
for path,text in base['files'].items():
 target=root/path;target.parent.mkdir(parents=True,exist_ok=True)
 relative=posixpath.relpath(link,posixpath.dirname(path) or '.')
 note=f"**Completed viscosity 0.01 disturbance comparison:** [Verified seed-202 result]({relative}) — 41 paired outputs through 0.40; disturbance {last['amplification']:.6f}× its starting size and still growing instantaneously. Single-seed finite-time result with known resolution limits.\n\n"
 head,sep,tail=text.partition('\n\n');assert sep and link not in text
 target.write_text(head+sep+note+tail,encoding='utf-8')
files=[]
for file in root.rglob('*'):
 if not file.is_file():continue
 data=file.read_bytes();path=file.relative_to(root).as_posix()
 if file.suffix=='.py':ast.parse(data)
 if file.suffix=='.json':json.loads(data)
 files.append({'path':path,'local':str(file.resolve()),'size':len(data),'sha256':hashlib.sha256(data).hexdigest(),'blob_sha':hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()})
paths={x['path'] for x in base['base']['blobs']}|{x['path'] for x in files};checks=0
for url in re.findall(r'\]\(([^)]+)\)',readme):
 if urlsplit(url).scheme or url.startswith('#'):continue
 target=posixpath.normpath(posixpath.join(rel,unquote(url.split('#')[0])))
 assert target in paths,(url,target)
 checks+=1
(P/'payload.json').write_text(json.dumps(files,indent=2)+'\n',encoding='utf-8')
(P/'quality-verification.json').write_text(json.dumps({'status':'passed','python_and_json_parse':'passed','new_report_links_checked':checks,'payload_files':len(files),'all_41_paired_fields_verified':True,'new_simulations':0},indent=2)+'\n',encoding='utf-8')
print(json.dumps({'payload_files':len(files),'bytes':sum(x['size'] for x in files),'links_checked':checks,'analysis':summary},indent=2))

