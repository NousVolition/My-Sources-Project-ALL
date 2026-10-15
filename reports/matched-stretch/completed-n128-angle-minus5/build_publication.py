"""Build reports from completed, independently checked saved results; no integration.
Usage: python -B build_publication.py
Requires the original local study directory and verification/base-inputs beside this file.
"""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, shutil, os, ast, re, posixpath
from urllib.parse import urlsplit, unquote
import numpy as np
P=Path(__file__).resolve().parent
S=Path('C:/Users/ziman/Documents/Codex/2026-10-05/ca/work/matched-stretch-study')
ROOT=P/'publication'
MAT='reports/matched-stretch/completed-n128-angle-minus5'
JOBS=['angle_degrees-n128-v-5.0']
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

def matrix_report():
 out=ROOT/MAT;out.mkdir(parents=True,exist_ok=True);runs={};ver={}
 for job in JOBS:
  r=read(S/'runs'/job/'result.json');v=read(P/f'verification-{job}.json');complete(r)
  assert v['status']=='passed' and v['fields_checked']==41 and v['result_sha256']==sha(S/'runs'/job/'result.json')
  runs[job]=r;ver[job]=v;cp(S/'runs'/job/'result.json',out/'diagnostics'/f'{job}.json');cp(P/f'verification-{job}.json',out/'verification'/f'{job}.json')
 controls=['baseline-n128-base','baseline-n128-half','baseline-n64-base','baseline-n64-half','angle_degrees-n64-v-5.0']
 for job in controls:
  r=read(S/'runs'/job/'result.json');complete(r);runs[job]=r;cp(S/'runs'/job/'result.json',out/'diagnostics'/f'{job}.json')
 for src,dst in [(S/'protocol.json','protocol.json'),(S/'numerics.py','reproduction/numerics.py'),(S/'run_suite.py','reproduction/run_suite.py'),(P/'verify_matrix.py','verify_matrix.py'),(Path(__file__),'build_publication.py')]:cp(src,out/dst)
 t=np.array([x['t'] for x in runs[JOBS[0]]['series']]);assert all(list(t)==[x['t'] for x in r['series']] for r in runs.values())
 labels={j:'axis tilt -5 degrees' for j in JOBS};colors=dict(zip(JOBS,plt.cm.tab10(np.arange(len(JOBS)))))
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for j in JOBS+['baseline-n128-base']:
  rr=runs[j]['series'];kw={'label':labels.get(j,'baseline axis tilt 0 degrees'),'color':colors.get(j,'black'),'lw':1.4}
  for ax,key in zip(axs.flat,['Wmax','energy','enstrophy','I']):ax.plot(t,[x[key]/rr[0][key] if key!='I' else x[key] for x in rr],**kw)
 for ax,title,yl in zip(axs.flat,['Peak vorticity','Kinetic energy','Enstrophy','Integrated peak vorticity'],['W/W(0)','E/E(0)','Z/Z(0)','I(t)']):ax.set(title=title,ylabel=yl)
 finish_axes(axs);fig.suptitle('Completed strain-axis controls at grid 128; each curve uses its own initial state');save(fig,out,'parameter-curves.png')
 fig,axs=plt.subplots(1,3,figsize=(11,4.5),layout='constrained');ix=np.arange(len(JOBS))
 for ax,key,title in zip(axs,['Wmax','energy','enstrophy'],['Peak vorticity','Energy','Enstrophy']):
  for n,offset,col in [(64,-.18,'#176b9c'),(128,.18,'#bb4d28')]:
   vals=[]
   for j in JOBS:
    rr=runs[j.replace('n128',f'n{n}')]['series'];vals.append(rr[-1][key]/rr[0][key])
   ax.bar(ix+offset,vals,.36,label=f'grid {n}',color=col)
  ax.set_xticks(ix,[labels[j] for j in JOBS],rotation=0,ha='right');ax.set(title=title+' at t=0.40',ylabel='Final / own initial');ax.grid(axis='y',alpha=.2);ax.legend()
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
 state=read(S/'task-state.json');js(out/'analysis.json',{'status':'complete_verified_batch','new_fields_checked':41,'new_simulations':0,'runs':summary,'original_matrix_snapshot':state,'comparisons':'Existing completed matching saved times only; controls reused, not rerun or newly independently reverified.'})
 table='\n'.join(f"| [{x['job']}](diagnostics/{x['job']}.json) | {x['dt']:.7g} | {x['W_ratio']:.7g} | {x['energy_ratio']:.7g} | {x['enstrophy_ratio']:.7g} | {x['spatial_screen_failures']}/41 |" for x in summary)
 limits='\n'.join(f"| {labels[x['job']]} | {x['final_global_chord_cells']:.5g} | {x['final_central_chord_cells']:.5g} | {x['final_high_band_enstrophy_fraction']:.5%} | {x['max_absolute_energy_budget_residual']:.3g} | {x['max_absolute_enstrophy_budget_residual']:.3g} |" for x in summary)
 links='\n'.join(f"- [{j} independent verification](verification/{j}.json)" for j in JOBS)
 text=f'''# Completed grid-128 strain-axis tilt -5 degrees — October 15, 2026

[Matched-stretch reports](../README.md) · [Previous completed parameter batch](../completed-n128-parameter-batch-20261015/README.md)

The original **angle_degrees-n128-v-5.0** run reached **0.40** with **41 saved fields**, all independently remeasured. The frozen supplied initial-field family, per-run viscosity, mean, domain, Heun method and zero force were preserved. No new simulations or horizon extensions were introduced.

This authorized one-factor case tilts the supplied raw strain axis by -5 degrees and retains viscosity 0.001, domain side 6, grid 128 and zero force. It preserves its original Heun step of 0.01/285 = 3.508771929824561e-05. The exact prescribed initial field was checked using the frozen initial-field routine with angle_degrees=-5.0. Its initial mean is {ver[JOBS[0]]['initial_mean']}; the zero-angle baseline has a different prescribed initial mean. Each trajectory conserves its own initial mean. Peak vorticity, energy and enstrophy curves are normalized by their own initial values; I(t) remains absolute. Absolute quantities remain in the raw diagnostics. This is a parameter-control experiment, not a disturbed/reference separation test.

| Run | Preserved dt | W(0.40)/W(0) | E(0.40)/E(0) | Z(0.40)/Z(0) | Failed spatial screens |
|---|---:|---:|---:|---:|---:|
{table}

![Parameter curves](parameter-curves.png)

## Existing controls and numerical limits

The matching completed -5-degree grid-64 control and completed zero-angle grid-128 baseline (both viscosity 0.001) are included at identical saved times. Completed baseline grid-64 and half-step records are also supplied as context. Each retains its own recorded viscosity, initial field and time step. The 64/128 endpoint comparison below does not establish spatial convergence; baseline half-step checks do not establish time convergence of this tilted trajectory.

![Matching grid controls at the endpoint](grid-endpoints.png)

**Raw strain is nonsmooth across periodic joins, and the initial maxima lie outside the central tube.** The published spatial screens require high-band energy below 0.1%, high-band enstrophy below 1%, and minimum half-peak chord at least six cells. These are screening criteria, not proof of resolution. The table distinguishes spatial limitations from budget consistency.

| Run | Final global width, cells | Final central width, cells | Final high-band Z | Maximum absolute E budget residual | Maximum absolute Z budget residual |
|---|---:|---:|---:|---:|---:|
{limits}

![Budgets and resolution checks](budgets-resolution.png)

Individual budgets use recorded Heun-stage integrals, checked against saved endpoint energies and enstrophies. Independent physical-space checks cover norms, means, divergence, enstrophy production and viscous losses. All declared numerical hard stops pass; neither small residuals nor a viscosity change proves a physical stability transition. Unresolved transfer to small scales remains a limitation.

![Final shell spectra](spectra.png)

The initial-return diagnostic measures distance from each trajectory's own initial velocity field. It does not establish recurrence, periodicity or a search over all state pairs. This unperturbed strain-axis sweep does not supply disturbance-growth or Lyapunov measurements.

![Initial-return diagnostics](initial-return.png)

## Verification and reproducibility

The [verifier](verify_matrix.py) runs as `python -B verify_matrix.py STUDY_DIRECTORY JOB_ID OUTPUT_JSON`. It checks all 41 saved fields per new run, the exact prescribed initial field, unchanged source/protocol hashes, mean, divergence, physical norms, stretching and losses, saved budget endpoints, spectra, widths and declared hard stops. Width and spectral remeasurement uses the frozen routine and is a consistency check, not a second discretization. Original large Fourier fields are required locally and are not embedded; their SHA-256 values are retained below.

{links}

[Raw diagnostics](diagnostics/) include the newly completed strain-axis run and existing completed controls. Existing controls are reused without rerunning or claiming a new independent full-field verification. [Analysis](analysis.json) records summary values and the current completion snapshot. [Frozen protocol](protocol.json), [numerical source](reproduction/numerics.py), [matrix controller](reproduction/run_suite.py) and [publication builder](build_publication.py) retain provenance.

At this report's snapshot, the original matrix has **{len(state['completed'])}/48 complete**, **{len(state['active'])} active**, **{len(state['pending'])} queued**, and **{len(state['failures'])} recorded failures**. Active jobs: {', '.join(state['active'])}. This is a completed batch report; the original matrix's final report awaits its remaining jobs. Historical intake reviews remain unchanged.
'''
 put(out/'README.md',text);return summary

matrix=matrix_report()
for path,text in base['files'].items():
 notes=[]
 for rel,label,detail in [(MAT,'Grid-128 strain-axis tilt -5 degrees','all 41 saved fields independently remeasured through 0.40; original matrix continues')]:
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
for rel in [MAT]:
 for url in re.findall(r'\]\(([^)]+)\)',(ROOT/rel/'README.md').read_text(encoding='utf-8')):
  if urlsplit(url).scheme or url.startswith('#'):continue
  target=posixpath.normpath(posixpath.join(rel,unquote(url.split('#')[0])))
  assert target in paths or any(x.startswith(target+'/') for x in paths),(url,target)
  checks+=1
js(P/'payload.json',files);js(P/'quality-verification.json',{'status':'passed','python_json_parse':'passed','links_checked':checks,'payload_files':len(files),'matrix_fields_checked':41,'new_simulations':0})
print(json.dumps({'payload_files':len(files),'bytes':sum(x['size'] for x in files),'links_checked':checks,'matrix':matrix},indent=2))
