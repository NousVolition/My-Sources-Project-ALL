"""Build the completed grid-256 baseline report from saved fields/diagnostics.
Run beside base-inputs.json and verification-baseline-n256-base.json.
No trajectory is advanced. Original local study files are required.
"""
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import urlsplit,unquote
import json,hashlib,shutil,os,ast,re,posixpath
import numpy as np
P=Path(__file__).resolve().parent
S=Path('C:/Users/ziman/Documents/Codex/2026-10-05/ca/work/matched-stretch-study')
ROOT=P/'publication';DEST='reports/matched-stretch/completed-n256-baseline'
JOB='baseline-n256-base';CONTROLS=['baseline-n64-base','baseline-n64-half','baseline-n128-base','baseline-n128-half']
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,text):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text,encoding='utf-8',newline='')
def js(p,x):put(p,json.dumps(x,indent=2)+'\n')
def cp(a,b):b.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(a,b)
base=read(P/'base-inputs.json');protocol=read(S/'protocol.json');v=read(P/f'verification-{JOB}.json')
assert v['status']=='passed' and v['fields_checked']==41 and v['through']==.4
for n,h in protocol['source_hashes'].items():assert sha(S/n)==h
assert v['protocol_sha256']==sha(S/'protocol.json')
runs={j:read(S/'runs'/j/'result.json') for j in [JOB]+CONTROLS}
for j,r in runs.items():assert r['status']=='complete' and r['last_output']==40 and len(r['series'])==41
assert v['result_sha256']==sha(S/'runs'/JOB/'result.json')
out=ROOT/DEST;out.mkdir(parents=True,exist_ok=True)
for j in runs:cp(S/'runs'/j/'result.json',out/'diagnostics'/f'{j}.json')
for src,dst in [(S/'protocol.json','protocol.json'),(S/'numerics.py','reproduction/numerics.py'),(S/'run_suite.py','reproduction/run_suite.py'),(P/'verify_matrix.py','verify_matrix.py'),(Path(__file__),'build_publication.py'),(P/f'verification-{JOB}.json',f'verification/{JOB}.json')]:cp(src,out/dst)
t=np.array([x['t'] for x in runs[JOB]['series']]);assert all(list(t)==[x['t'] for x in r['series']] for r in runs.values())
os.environ.setdefault('MPLCONFIGDIR',str(P/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
colors={64:'#176b9c',128:'#bb4d28',256:'#6a4195'}
def save(fig,name):fig.savefig(out/name,dpi=155);plt.close(fig)
def axes_finish(axes):
 for ax in np.asarray(axes).flat:ax.set_xlabel('Saved model time');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for n in [64,128,256]:
 rr=runs[f'baseline-n{n}-base']['series']
 for ax,key in zip(axs.flat,['Wmax','energy','enstrophy','I']):ax.plot(t,[x[key]/rr[0][key] if key!='I' else x[key] for x in rr],label=f'grid {n}',color=colors[n])
for ax,title,yl in zip(axs.flat,['Peak vorticity','Kinetic energy','Enstrophy','Integrated peak vorticity'],['W/W(0)','E/E(0)','Z/Z(0)','I(t)']):ax.set(title=title,ylabel=yl)
axes_finish(axs);fig.suptitle('Completed baseline grids; each normalized curve uses its own initial value');save(fig,'grid-curves.png')
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for n in [64,128]:
 a=runs[f'baseline-n{n}-base']['series'];b=runs[f'baseline-n{n}-half']['series']
 for ax,key in zip(axs.flat,['Wmax','energy','enstrophy','I']):
  vals=np.array([abs(y[key]-x[key])/max(abs(x[key]),1e-300) if x[key] else 0. for x,y in zip(a,b)])
  ax.semilogy(t,np.maximum(vals,1e-16),label=f'grid {n}',color=colors[n])
for ax,title in zip(axs.flat,['Peak vorticity','Energy','Enstrophy','Integrated peak vorticity']):ax.set(title=title,ylabel='Absolute half-step difference / base value')
axes_finish(axs);fig.suptitle('Completed 64/128 time-step controls; grid-256 half-step run remains incomplete');save(fig,'half-step-controls.png')
rr=runs[JOB]['series'];first=rr[0];last=rr[-1]
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,key,label in [(axs[0,0],'energy_budget_relative_residual','energy'),(axs[0,1],'enstrophy_budget_relative_residual','enstrophy')]:ax.semilogy(t,np.maximum(np.abs([x[key] for x in rr]),1e-16),label=label,color=colors[256]);ax.set(title=label.capitalize()+' budget',ylabel='Absolute relative residual')
for key,label in [('high_band_energy_fraction','energy'),('high_band_enstrophy_fraction','enstrophy')]:axs[1,0].semilogy(t,[x[key] for x in rr],label=label)
axs[1,0].axhline(.001,color='gray',ls=':',label='energy screen: 0.1%');axs[1,0].axhline(.01,color='gray',ls='--',label='enstrophy screen: 1%');axs[1,0].set(title='High-band fractions',ylabel='Fraction')
for key,label in [('width_at_global_peak','global peak'),('width_at_central_roi_peak','central peak')]:axs[1,1].plot(t,[x[key]['minimum_chord_cells'] for x in rr],label=label)
axs[1,1].axhline(6,color='gray',ls='--',label='six-cell screen');axs[1,1].set(title='Minimum half-peak chords',ylabel='Grid cells');axes_finish(axs);save(fig,'budgets-resolution.png')
fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
for sp in v['spectra']:
 k=(np.array(sp['k_edges'][:-1])+np.array(sp['k_edges'][1:]))/2
 for ax,key in zip(axs,['energy','enstrophy']):
  y=np.array(sp[key]);mask=y>1e-22;ax.semilogy(k[mask],y[mask],label=f"t={sp['t']:.2f}")
for ax,key in zip(axs,['energy','enstrophy']):ax.set(title='Grid-256 shell '+key,xlabel='Angular wavenumber magnitude',ylabel='Shell-integrated amount');ax.grid(alpha=.2);ax.legend(fontsize=8)
save(fig,'spectra.png')
fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained');ax.plot(t,[x['distance_from_initial_l2'] for x in v['fields']],label='grid 256',color=colors[256]);ax.set(title='Distance from the grid-256 initial field',ylabel='L2 distance / own initial L2 norm');axes_finish([ax]);save(fig,'initial-return.png')
summary=[]
for j,r in runs.items():
 a=r['series'][0];b=r['series'][-1]
 summary.append({'job':j,'n':r['job']['n'],'dt':r['dt'],'nu':r['job']['nu'],'new_independent_field_verification':j==JOB,'W_ratio':b['Wmax']/a['Wmax'],'energy_ratio':b['energy']/a['energy'],'enstrophy_ratio':b['enstrophy']/a['enstrophy'],'I':b['I'],'initial_mean':a['mean_velocity'],'initial_Wmax':a['Wmax'],'initial_energy':a['energy'],'initial_enstrophy':a['enstrophy'],'initial_peak_location':a['peak_location'],'spatial_screen_failures':sum(not x['resolved_screen'] for x in r['series']),'final_high_band_energy_fraction':b['high_band_energy_fraction'],'final_high_band_enstrophy_fraction':b['high_band_enstrophy_fraction'],'final_global_width_cells':b['width_at_global_peak']['minimum_chord_cells'],'final_central_width_cells':b['width_at_central_roi_peak']['minimum_chord_cells'],'max_absolute_energy_budget_residual':max(abs(x['energy_budget_relative_residual']) for x in r['series']),'max_absolute_enstrophy_budget_residual':max(abs(x['enstrophy_budget_relative_residual']) for x in r['series'])})
state=read(S/'task-state.json');new=summary[0]
comparison={}
for n in [64,128]:
 a=runs[f'baseline-n{n}-base']['series'];b=runs[f'baseline-n{n}-half']['series'];comparison[str(n)]={k:{'endpoint_relative_difference':abs(b[-1][k]-a[-1][k])/abs(a[-1][k]),'maximum_relative_difference':max(abs(y[k]-x[k])/abs(x[k]) if x[k] else 0. for x,y in zip(a,b))} for k in ['Wmax','energy','enstrophy','I']}
js(out/'analysis.json',{'status':'complete_verified_batch','new_fields_checked':41,'new_simulations':0,'runs':summary,'completed_half_step_comparisons':comparison,'original_matrix_snapshot':state,'comparisons':'Only completed controls at identical saved times; existing control fields are not rerun or newly independently reverified.'})
table='\n'.join(f"| [{x['job']}](diagnostics/{x['job']}.json) | {x['dt']:.7g} | {x['W_ratio']:.7g} | {x['energy_ratio']:.7g} | {x['enstrophy_ratio']:.7g} | {x['spatial_screen_failures']}/41 |" for x in summary)
means='\n'.join(f"| {n} | "+', '.join(f'{z:.12g}' for z in runs[f'baseline-n{n}-base']['series'][0]['mean_velocity'])+' |' for n in [64,128,256])
put(out/'README.md',f'''# Completed grid-256 baseline — October 15, 2026

[Matched-stretch reports](../README.md)

The original **baseline-n256-base** reached **0.40** with **41 saved fields**, all independently remeasured. The supplied initial-field formula and original grid-specific starting array, viscosity **0.001**, domain side **6**, retained mean, Heun method, zero force and endpoint were preserved. Its time step remains **0.01/409 = {runs[JOB]['dt']:.17g}**. No trajectory was repeated or extended.

At 0.40, peak vorticity is **{new['W_ratio']:.7g} times** its starting value; energy is **{new['energy_ratio']:.7g} times** its starting value, and enstrophy is **{new['enstrophy_ratio']:.7g} times** its starting value. These are finite-time numerical observations. **38 of 41 spatial screens fail**, so the larger grid does not establish spatial convergence or resolved physical growth.

## Completed grid and time-step controls

All table entries are complete through 0.40 and compared at the same 41 saved times. Every normalized curve uses its own initial value; absolute values remain in the diagnostics. The grid-64 and grid-128 results and their half-step controls are reused from existing completed runs. Only the newly completed grid-256 baseline receives a new full-field audit in this publication.

| Run | Preserved dt | W(0.40)/W(0) | E(0.40)/E(0) | Z(0.40)/Z(0) | Failed spatial screens |
|---|---:|---:|---:|---:|---:|
{table}

![Completed grid curves](grid-curves.png)

The existing **grid-256 half-step control is still running**. Its endpoint comparison is unavailable and is not inferred from the completed lower-grid controls. Grid sensitivity and time-step sensitivity are separate questions; agreement at one grid or in one quantity does not establish convergence for the other quantities.

![Completed half-step controls](half-step-controls.png)

**The supplied raw strain is nonsmooth across periodic joins, and initial maxima lie outside the central tube.** Sampling this initial-field formula on different grids gives different discrete means. Those means are preserved per run; they are not silently equalized or subtracted in these comparisons.

| Grid | Preserved initial mean velocity (x, y, z) |
|---:|---|
{means}

## Budgets, spectra and width limits

For grid 256, the final global minimum half-peak chord is **{new['final_global_width_cells']:.5g} cells**, the central chord is **{new['final_central_width_cells']:.5g} cells**, and the high-band energy/enstrophy fractions are **{new['final_high_band_energy_fraction']:.5%} / {new['final_high_band_enstrophy_fraction']:.5%}**. The declared spatial screens require energy below 0.1%, enstrophy below 1%, and minimum chord at least six cells. Passing such screens would itself not prove resolution.

The maximum absolute relative energy/enstrophy budget residuals are **{new['max_absolute_energy_budget_residual']:.5g} / {new['max_absolute_enstrophy_budget_residual']:.5g}**. Saved Heun-stage integrals are checked against endpoint budgets; physical norms, stretching production and viscous losses are independently remeasured. All declared numerical hard stops pass. Small budget residuals do not repair spatial under-resolution or prove a physical stability threshold.

![Budgets and spatial screens](budgets-resolution.png)

![Remeasured shell spectra at five saved times](spectra.png)

The distance below uses each saved field and this run's own initial field. It is an initial-return diagnostic, not a search for recurrence among every pair of states or evidence of periodicity. This unperturbed baseline does not provide a disturbed/reference separation rate or an asymptotic Lyapunov exponent.

![Initial-return diagnostic](initial-return.png)

## Verification and provenance

[Independent verification](verification/{JOB}.json) records all 41 field SHA-256 hashes, exact initial-field matching, unchanged source/protocol hashes, conserved mean, scaled divergence, physical energy/enstrophy/vorticity, stretching production, viscous losses, budget consistency, spectra, widths and hard-stop checks. Physical-space FFT measurements use NumPy independently of the simulation's SciPy transforms. Spectral/width observation remeasurement also uses the frozen routine and is a consistency check, not a second discretization. Complete original Fourier arrays are required locally and are not embedded in this report.

Run the [verifier](verify_matrix.py) as `python -B verify_matrix.py STUDY_DIRECTORY baseline-n256-base OUTPUT_JSON`. [Raw diagnostics](diagnostics/), [analysis](analysis.json), [frozen protocol](protocol.json), [numerical source](reproduction/numerics.py), [controller](reproduction/run_suite.py) and [publication builder](build_publication.py) provide provenance. No active numerical source was modified.

The original matrix snapshot is **{len(state['completed'])}/48 complete**, **{len(state['active'])} active**, **{len(state['pending'])} queued**, and **{len(state['failures'])} recorded failures**. Active: {', '.join(state['active'])}. Its final report awaits the remaining authorized jobs. Historical intake reviews remain unchanged.
''')
for path,text in base['files'].items():
 assert DEST not in text
 link=posixpath.relpath(DEST+'/README.md',posixpath.dirname(path) or '.')
 note=f'**Completed October 15, 2026:** [Grid-256 baseline]({link}) — all 41 saved fields independently remeasured through 0.40; 38 spatial screens fail and the original half-step control continues. Raw-strain, initial-peak and resolution limits remain.\n\n'
 nl='\r\n' if '\r\n' in text else '\n';head,sep,tail=text.partition(nl*2);assert sep;put(ROOT/path,head+sep+note.replace('\n',nl)+tail)
payload=[]
for f in ROOT.rglob('*'):
 if not f.is_file():continue
 data=f.read_bytes();rel=f.relative_to(ROOT).as_posix()
 if f.suffix=='.py':ast.parse(data)
 if f.suffix=='.json':json.loads(data)
 payload.append({'path':rel,'local':str(f.resolve()),'size':len(data),'sha256':hashlib.sha256(data).hexdigest(),'blob_sha':hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()})
paths={x['path'] for x in base['base']['blobs']}|{x['path'] for x in payload};count=0
for url in re.findall(r'\]\(([^)]+)\)',(out/'README.md').read_text(encoding='utf-8')):
 if urlsplit(url).scheme or url.startswith('#'):continue
 target=posixpath.normpath(posixpath.join(DEST,unquote(url.split('#')[0])));assert target in paths or any(x.startswith(target+'/') for x in paths),(url,target);count+=1
js(P/'payload.json',payload);js(P/'quality-verification.json',{'status':'passed','python_json_parse':'passed','relative_links_checked':count,'payload_files':len(payload),'new_fields_checked':41,'new_simulations':0})
print(json.dumps({'payload_files':len(payload),'bytes':sum(x['size'] for x in payload),'links':count,'new_run':new},indent=2))
