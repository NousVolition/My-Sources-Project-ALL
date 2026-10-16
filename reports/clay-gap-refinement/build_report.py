"""Create scientific plots and a standalone viewer from completed recordings."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,html,hashlib,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_one import load
R=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
def find(name):
 paths=[p for p in R.glob('*/data/*.json') if p.stem==name]+[p for p in (R/'baseline-data').glob('*.json') if p.stem==name]
 return load(paths[0]) if paths else None
def error_chart(selected,path,title):
 fig,axes=plt.subplots(1,2,figsize=(11,max(4.8,1+len(selected))),layout='constrained')
 for ax,key,label in zip(axes,['gap_arclength_rms','surface_gap_rms'],['Surface-length-weighted RMS gap','Original equal-edge-sampling RMS gap']):
  yy=np.arange(len(selected));values=[100*c[key]['maximum_scaled'] for _,c in selected];ax.barh(yy,values,color=['#946946' if c['a']=='shared144-one' else '#167d8d' for _,c in selected]);ax.set_yticks(yy,[l+f"\n{c['start']:g}–{c['end']:g}" for l,c in selected]);ax.invert_yaxis();ax.axvline(5,color='#5e5e5e',linestyle='--',label='Unchanged 5% target');ax.set_xlabel('Maximum gap disagreement (%)');ax.set_title(label);ax.set_xlim(0,max(8,max(values,default=0)*1.22));ax.legend(loc='lower right',frameon=False)
  for i,x in enumerate(values):ax.text(x+.08,i,(f'{x:.2f}%' if x>=.1 else f'{x:.4f}%'),va='center')
 fig.suptitle(title);fig.savefig(path);plt.close(fig)
def main():
 a=json.loads((R/'analysis.json').read_text());figdir=R/'figures';figdir.mkdir(exist_ok=True)
 cs=a['comparisons']
 labels=[('Original 144 / 168, one-way','shared144-one','extend168-one'),('Revised 72 / 96, one-way','length100-72-one','length100-96-one'),('Revised 144 / 168, one-way','length100-144-one','length100-168-one'),('Revised 144 / 168, reciprocal','length100-144-two','length100-168-two')]
 selected=[(label,next((c for c in cs if c['a']==x and c['b']==y),None)) for label,x,y in labels];selected=[(l,c) for l,c in selected if c]
 error_chart(selected,figdir/'mesh-gaps.png','Full saved histories; intervals shown for every comparison')
 harder_labels=[('Harder 72 / 96, one-way','harder72-one','harder96-one'),('Harder 96 / 120, one-way','harder96-one','harder120-one'),('Harder 96 / 120, reciprocal','harder96-two','harder120-two')]
 harder=[(l,next((c for c in cs if c['a']==x and c['b']==y),None)) for l,x,y in harder_labels];harder=[(l,c) for l,c in harder if c]
 error_chart(harder,figdir/'harder-gaps.png','Harsher forcing: α = 256, ρ = 1000, viscosity = 0.009375')
 selected+=harder
 fig,axes=plt.subplots(2,1,figsize=(9,6.4),layout='constrained')
 for name,color,style in [('length100-144-one','#c36c37','--'),('length100-168-one','#147c8d','-'),('length100-144-two','#c36c37','--'),('length100-168-two','#147c8d','-')]:
  s=find(name)
  if not s:continue
  r,_,t,_,d=s;ax=axes[1 if r['parameters']['feedback'] else 0];ax.plot(t,[v['gap_arclength_rms'] for v in d],style,color=color,label=f"{r['parameters']['angles']} × {r['parameters']['layers']}")
 for ax,label in zip(axes,['One-way Lorenz → material','Reciprocal Lorenz ↔ material']):ax.set_ylabel('Weighted RMS gap\n(model length units)');ax.set_xlabel('Model time (dimensionless)');ax.set_title(label);ax.legend(frameon=False);ax.grid(alpha=.2)
 axes[-1].set_xlabel('Model time (dimensionless)');fig.savefig(figdir/'gap-histories.png');plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
 for ax,(names,title) in zip(axes,[(['length100-72-one','length100-72-two'],'Capacity 0.01'),(['middle72-one','middle72-two'],'Capacity 100'),(['weak72-one','weak72-two'],'Capacity 10,000')]):
  for name,color in zip(names,['#c36c37','#147c8d']):
   if name=='weak72-two' and find('weak72-two-tighter'):name='weak72-two-tighter'
   s=find(name)
   if s:r,_,t,y,_=s;ax.plot(t,y[:,-7],color=color,label='Reciprocal' if r['parameters']['feedback'] else 'One-way')
  ax.set_title(title);ax.set_xlabel('Model time');ax.set_ylabel('Lorenz X (model units)');ax.legend(frameon=False);ax.grid(alpha=.2)
 fig.suptitle('Feedback strength changes the drive; default feedback is strongly damping');fig.savefig(figdir/'reciprocal-drive.png');plt.close(fig)
 lines=['# Hug: gap refinement and longer Lorenz stress tests','',
 'This follow-up preserves the published model and tests a separate proposed material revision. The revision adds a fixed reference-gradient jump energy to resist deformation concentrated into smaller and smaller regions. **It changes the constitutive model; it is not a correction that makes the old model exact.**','',
 'The nominal test uses coupling α = 128, Lorenz ρ = 500, viscosity 0.01875, preload 0.35, initial Lorenz state (−8, 8, 27), the same tilted ellipse and square boundary, and the same area and shape floors. All quantities are dimensionless.','',
 '[Raw-array release and restoration catalog](data-release.json). Numerical arrays are distributed as a verified GitHub release archive; code, summaries and the viewer remain in this source folder. Restore arrays before running the analysis or verification commands below.','',
 '## Measured mesh differences','',
 'Both statistics use the unchanged denominator max(0.01, |fine gap|), the unchanged 5% target, and all common saved times (normally every 0.005). No gap rescaling, clipping, fitted correction or remeshing was used. A passed gap test alone does not establish convergence of imprint or all state variables.','',
 'The measured gap is radial separation from the ellipse, not exact closest-point distance. The weighted RMS integrates its square along the current inner boundary with 32-point Gauss quadrature per edge. The original statistic retains 17 equally spaced samples per edge (including endpoints). The conservative whole-edge clearance used for collision protection is a separate quantity. Reported maxima are over saved times, not certified maxima over continuous time.','',
 '| Model / coupling | Meshes | Common interval | Weighted gap max | Original gap max | Full gap check |','|---|---|---|---:|---:|---|']
 for label,c in selected:lines.append(f"| {label} | {c['coarse_mesh'][0]}×{c['coarse_mesh'][1]} / {c['fine_mesh'][0]}×{c['fine_mesh'][1]} | {c['start']:g}–{c['end']:g} | {100*c['gap_arclength_rms']['maximum_scaled']:.4f}% | {100*c['surface_gap_rms']['maximum_scaled']:.4f}% | {'Pass' if c['full_mesh_check_passed'] else 'Unresolved'} |")
 fine=next(c for c in cs if c['a']=='length100-144-one' and c['b']=='length100-168-one');medium=next(c for c in cs if c['a']=='length100-72-one' and c['b']=='length100-96-one')
 lines+=['',f"Imprint differs by at most **{100*fine['imprint']['maximum_scaled']:.3f}%** on the revised 144/168 comparison through time 1.0. On 72/96 through 0.5 it still differs by **{100*medium['imprint']['maximum_scaled']:.3f}%**, although both gap statistics pass."]
 hard=next(c for c in cs if c['a']=='harder96-one' and c['b']=='harder120-one');lines+=['',f"Under the harsher forcing, the 96/120 imprint difference is still **{100*hard['imprint']['maximum_scaled']:.3f}%**. Passing the gap target does not remove this remaining uncertainty in plastic memory."]
 lines+=['','![Mesh comparisons](figures/mesh-gaps.png)','','![Harsher forcing](figures/harder-gaps.png)','','![Longer histories](figures/gap-histories.png)','','## What changed and why','',
 'The original squeeze diagnostic shows active same-boundary self-contact terms, many between nearby material neighbors, around the previous error peaks at time 0.125–0.130. This is evidence of a difficult localized deformation/contact regime, not proof that any one force term alone caused the old disagreement. Strict BDF and independent Radau restarts test time error separately; every stopped attempt remains in the archive.','',
 'The proposed reference energy is E_reg = (μ ℓ² / 2) Σ_e (L_e / d_e) ||F_left − F_right||². Here L_e is an interior edge length in the reference mesh, d_e is the normal distance between adjacent reference triangle centroids, and F is each triangle’s deformation gradient. The fixed material length ℓ = 0.1 is one quarter of the initial annulus thickness 0.4; it was declared before examining the revised mesh results. Additional ℓ = 0.05 and 0.2 runs expose sensitivity to this new modeling choice. The positive energy vanishes for every affine deformation and has an exact constant sparse Hessian. Plastic rest lengths and their fading memory retain the original evolution law. This term is elastic and can alter how the imprint relaxes.','',
 'The numerical backend eliminates diagonal internal-variable rows and the small Lorenz border by exact Schur complements. Independent full sparse real/complex solves check the result. This backend improvement leaves the equations unchanged. A half-unit segmented reference drive preserves every value on the published first interval exactly when extending the one-way simulation. Independent DOP853/Radau reference calculations check input precision.','',
 '## Actual runs','',f"{a['completed']} of {a['attempts']} new attempts completed their requested intervals. All new recorded histories passed geometry: **{a['all_geometry_passed']}**; separate material/Lorenz/combined energy audits passed: **{a['all_energy_passed']}**. The energy target is 0.1%. Stopped attempts are reported with their actual endpoints, never as completed tests.",'',
 '| Recording | Mesh | Direction | Actual / requested end | Energy error max | Status |','|---|---|---|---|---:|---|']
 for r in a['runs']:lines.append(f"| [{r['name']}]({r['path']}) | {r['parameters']['angles']}×{r['parameters']['layers']} | {'reciprocal' if r['parameters']['feedback'] else 'one-way'} | {r['end']:.8g} / {r['requested_end']:g} | {100*max(r['energy_errors'].values()):.5g}% | {'Complete' if r['complete'] else 'Stopped: '+r['failure']['message']} |")
 lines+=['','## Accuracy and uncertainty','',
 'Numerical uncertainty is represented by mesh differences, independent time-method/tolerance comparisons, recomputed energy balance, positive clearances and certified straight paths between accepted states. These are deterministic error indicators, not statistical confidence intervals. See [analysis.json](analysis.json), [exact execution settings](data/executed-protocol.json), [factor checks](data/factor-check.json), [affine/gradient/reference checks](data/regularization-check.json), and [all saved-state audits](data/audits/).','',
 '| Independent time comparison | Interval checked | Scaled full-state difference | Outcome |','|---|---|---:|---|']
 for c in a['time_controls']:lines.append(f"| {c['a']} / {c['b']} | {c['start']:g}–{c['end']:g} | {100*c.get('state_scaled',float('inf')):.6g}% | {'Pass' if c.get('passed') else 'Unresolved'} |")
 lines+=['','The looser C = 10,000 BDF/Radau comparison fails the 0.01% scaled-state target. The higher-precision Radau pair (rtol 2e−8 / 2e−9) is checked separately, and the viewer uses the tighter recording. The original looser result remains in the table.','',
 '| Sensitivity (changes, not error bars) | Comparison | Max weighted gap change | Max imprint change |','|---|---|---:|---:|']
 for c in a.get('sensitivities',[]):lines.append(f"| {c['kind']} | {c['a']} / {c['b']} | {100*c['gap_arclength_rms']['maximum_scaled']:.4g}% | {100*c['imprint']['maximum_scaled']:.4g}% |")
 lines+=['','The material-length sweep changes the equations, so its differences must not be treated as numerical error or used to infer which length is physically correct. The small initial-perturbation test uses the looser integrator; its observed change is comparable to that configuration’s numerical uncertainty and does not isolate physical sensitivity.']
 lines+=['',
 'A short checkpoint comparison tests only the stated interval and retains the parent’s accumulated uncertainty. It cannot validate earlier evolution or a later full run. The original-model peak controls assess time error locally; they do not repair the old 7.72% / 10.19% spatial gap difference. Agreement in gap also does not imply agreement in imprint: the revised 72/96 imprint comparison must be read from analysis.json.','',
 '## Both directions of influence','',
 'The Lorenz state drives radial material loading. The reciprocal rule adds λ A_dot / C to X_dot, with λ = α/10, A = Σ_i m_i |x_i| and capacity C. Material port power is −λ X A_dot; the Lorenz port receives the opposite power. Material deformation, plastic dissipation and its changing motion therefore affect the drive through this explicitly chosen scalar coupling. Native Lorenz power is tracked separately; the combined system is not assumed to conserve energy.','',
 'The default C = 0.01 causes very strong damping of X. Intermediate C = 100 and larger-capacity C = 10,000 tests explore moving-drive feedback. Nearby initial conditions test finite-horizon sensitivity; they do not establish a Lyapunov exponent or long-term statistical convergence. The harsher tests double α and ρ and halve viscosity. Their mesh accuracy must not be inferred from the nominal forcing tests.','',
 '![Reciprocal input histories](figures/reciprocal-drive.png)','','## Reproduce','',
 'Use Python 3.12 and the pinned requirements. Every recording includes source SHA-256 hashes, all parameters, solver tolerances, budgets, actual endpoint, work counters and raw states. Solvers refuse to overwrite an attempt. BLAS libraries use one thread. [baseline.json](baseline.json) identifies the immutable previous GitHub commit; baseline sources and data are preserved locally. Checkpoint prefixes preserve exact raw states and explicitly retain the original full-run audit scope.','',
 '```sh','python -m pip install -r requirements.txt','python restore_data.py','python restore_data.py --verify-only','python verify.py','python check_factor.py','python check_regularization.py','python check_derivatives.py','python check_quadrature.py',
 'python run.py my-new-case --n 168 --layers 15 --end 1 --seconds 10800',
 'python run.py my-two-way-case --n 168 --layers 15 --end 3 --two --seconds 3600',
 'python replay.py solver/data/harder72-two.json /absolute/path/to/a-new-empty-location','```','',
 'The replay command creates a result-free copy and checks every saved array for bitwise equality. It requires a destination that does not yet exist. Failed/resource-limited runs can complete to different endpoints on different hardware; their declared resource caps are not physical failure thresholds. Replaying the first one-way interval requires substantial compute time.','',
 'The check commands preserve the stored publication evidence by default; `--write` explicitly replaces a check’s evidence file in a new build. To regenerate analysis, figures or the report, use a separate copy and run `python analyze.py` and `python build_report.py`. A rebuilt report needs its browser check and a new manifest; it is not the immutable original publication. Study-local Git attributes preserve the audited line endings. Bitwise replay was verified on the recorded runtime; other numerical libraries may require tolerance comparisons.','',
 '## Scope','',
 'This is a two-dimensional abstract yielding spring material with a proposed nonlocal elastic regularization, fixed obstacle and square boundary. It has not been calibrated against clay measurements. Lorenz is an adopted reduced dynamical drive, not a solved Navier–Stokes flow or a quantum calculation. These tests do not prove regularity of Navier–Stokes, establish real-clay accuracy, or demonstrate nuclear quantum effects. The original model’s unresolved results remain visible rather than being replaced by the revised model’s results.','',
 '[Open the standalone saved-state viewer](report.html). The viewer selects actual saved states; graphics are downsampled and rounded for display only. Numerical analysis uses complete double-precision raw recordings.','']
 (R/'README.md').write_text('\n'.join(lines),encoding='utf-8')
 # Exact protocol reconstructed from immutable metadata, not a suggested plan.
 protocol=[]
 for r in a['runs']:
  p=R/r['path'];v=json.loads(p.read_text());protocol.append({k:v[k] for k in ['name','parameters','initial','start_time','requested_end','completed_end','record_step','parent','failure','audit','source_hashes']})
 (R/'data/executed-protocol.json').write_text(json.dumps(dict(thread_profile={'numpy':1,'scipy':1},attempts=protocol),indent=2)+'\n')
 viewer=[]
 for title,one,two in [('Revised material','length100-168-one','length100-168-two'),('Original material','original168-long-one','original168-long-two'),('Revised: capacity 100','middle72-one','middle72-two'),('Revised: capacity 10,000','weak72-one','weak72-two')]:
  if two=='weak72-two' and find('weak72-two-tighter'):two='weak72-two-tighter'
  A,B=find(one),find(two)
  if A is None and title=='Revised material':A,B=find('length100-96-one'),find('length100-96-two')
  if A is None or B is None:continue
  limit=min(A[2][-1],B[2][-1]);common=[(i,int(np.argmin(abs(B[2]-t)))) for i,t in enumerate(A[2]) if t<=limit+1e-10 and min(abs(B[2]-t))<1e-10];ids=np.unique(np.linspace(0,len(common)-1,min(31,len(common))).astype(int));frames=[]
  for idx in ids:
   i,j=common[idx];frame={'time':float(A[2][i]),'panels':[]}
   for s,k in [(A,i),(B,j)]:
    r,z,t,y,d=s;n=len(z['initial_positions']);frame['panels'].append(dict(xy=np.round(y[k,:2*n].reshape(n,2),5).tolist(),lorenz=np.round(y[:k+1,-7:-4],6).tolist(),gap=d[k]['gap_arclength_rms'],imprint=d[k]['imprint']))
   frames.append(frame)
  viewer.append(dict(title=title,records=[A[0]['name'],B[0]['name']],n=A[0]['parameters']['angles'],layers=A[0]['parameters']['layers'],triangles=A[1]['triangles'].tolist(),outline=np.round(A[1]['outline'],6).tolist(),frames=frames))
 (R/'data/playback.json').write_text(json.dumps(viewer,separators=(',',':'))+'\n')
 template=(R/'report-template.html').read_text(encoding='utf-8');assert template.count('PLAYBACK_DATA')==1
 table=''.join(f"<tr><td>{html.escape(l)}</td><td>0–{c['end']:g}</td><td>{100*c['gap_arclength_rms']['maximum_scaled']:.4f}%</td><td>{100*c['surface_gap_rms']['maximum_scaled']:.4f}%</td></tr>" for l,c in selected)
 (R/'report.html').write_text(template.replace('PLAYBACK_DATA',json.dumps(viewer,separators=(',',':'))).replace('RESULT_ROWS',table),encoding='utf-8')
 print(json.dumps(dict(runs=len(a['runs']),viewer_cases=len(viewer),figures=4)))
if __name__=='__main__':main()
