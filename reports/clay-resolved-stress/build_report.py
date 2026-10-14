"""Make figures, result prose and recorded-state playback from verified raw data."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,argparse,html,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import load,sha,Mesh,Parameters
ROOT=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--inline-output',type=Path);args=ap.parse_args();a=json.loads((ROOT/'analysis.json').read_text());by={r['name']:r for r in a['runs']}
 def series(name):return load(ROOT/by[name]['path'] if name in by else ROOT/'parents'/(name+'.json'))
 text=[f"**{a['completed']} of {a['attempts']} new integrations completed** their requested interval. Inherited parent histories are not counted as new integrations.",'','| Attempt | Mesh | Reached / requested | Status |','|---|---:|---:|---|']
 for r in a['runs']:text.append(f"| {r['name']} | {r['parameters']['angles']} × {r['parameters']['layers']} | {r['end']:.9g} / {r['requested_end']:g} | {r['failure']['message'] if r['failure'] else 'completed'} |")
 text+=['','### Mesh accuracy','','| Comparison | Common saved interval | Maximum surface-weighted gap difference | Maximum original gap difference | Full interval mesh check |','|---|---:|---:|---:|---|']
 for c in a['comparisons']:text.append(f"| {c['a']} / {c['b']} | [{c['start']:g}, {c['end']:g}] | {100*c['gap_arclength_rms']['maximum_scaled']:.5g}% | {100*c['surface_gap_rms']['maximum_scaled']:.5g}% | {c['full_mesh_check_passed']} |")
 text+=['','The original 48 / 72 disagreement is retained as historical evidence. A newly refined pair only passes when both full runs finish and both gap metrics remain below the same 5% target. Shorter or stopped histories cannot validate a longer requested interval. One-way drive histories also need 0.01% agreement.','', '| One-way comparison | Maximum scaled Lorenz-drive difference | Drive target met |','|---|---:|---|']
 for c in a['comparisons']:
  if c['a'].endswith('-one'):text.append(f"| {c['a']} / {c['b']} | {100*c['lorenz_scaled']:.5g}% | {c['drive_agreement_passed']} |")
 text+=['', '### Independent numerical accuracy','','| Control | Independently recomputed interval | Maximum scaled physical-state difference | 0.01% target |','|---|---:|---:|---|']
 for c in a['time_controls']:text.append(f"| {c['a']} / {c['b']} | [{c['start']:.9g}, {c['end']:g}] | {100*c['state_scaled']:.5g}% | {c['passed']} |")
 text+=['','State normalization is componentwise max(1, peak absolute reference component), excluding the four work bookkeeping variables. Inherited parent prefixes are excluded. Controls include a different time integrator, tighter tolerances, an equivalent force assembly and the independent-input formulation. They are not all isolated time-step experiments. A stopped comparison or a short control cannot certify a later interval.','', '### Both directions of influence','','| Mesh and common end | Material motion change with feedback | Lorenz Z fluctuation change with feedback |','|---|---:|---:|']
 for p in a['pairs']:text.append(f"| {p['mesh']} points, {p['end']:g} | {p['motion_change_percent']:+.5g}% | {p['z_std_change_percent']:+.5g}% |")
 text+=['','Material motion RMS comes from accumulated viscous dissipation divided by viscosity, reference area and duration. Lorenz Z fluctuations use time-weighted saved samples. The chosen feedback law introduces strong damping; these changes do not establish general chaos suppression.','', '### Geometry and energy','']
 if a['runs']:text+=[f"All geometry audits passed: **{a['all_geometry_passed']}**. All material, Lorenz and combined energy audits passed: **{a['all_energy_passed']}**. Maximum scaled full-history balance residual: **{max(max(r['energy_errors'].values()) for r in a['runs']):.6g}**. Minimum certified area ratio **{min(r['minimum_area'] for r in a['runs']):.6g}**, deformation quality **{min(r['minimum_quality'] for r in a['runs']):.6g}**, relative cell altitude **{min(r['minimum_altitude'] for r in a['runs']):.6g}**."]
 text+=['','![Same-forcing gap histories](figures/gap-histories.png)','', '![Gap disagreement under refinement](figures/mesh-errors.png)','', '![Effects of reciprocal coupling](figures/reciprocal.png)']
 readme=(ROOT/'README.md').read_text(encoding='utf-8');lo=readme.index('<!-- RESULTS -->')+len('<!-- RESULTS -->');hi=readme.index('<!-- END RESULTS -->');(ROOT/'README.md').write_text(readme[:lo]+'\n\n'+'\n'.join(text)+'\n\n'+readme[hi:],encoding='utf-8')
 protocol=[json.loads((ROOT/r['path']).read_text()) for r in a['runs']];protocol=[{k:r[k] for k in ['name','parameters','initial','requested_end','start_time','record_step','parent','audit','source_hashes']} for r in protocol]
 profiles={r['name']:dict(numpy_blas_threads=1 if r['name'].startswith('workflow-smoke') or 'indexed_contact.py' in r['source_hashes'] else 24,scipy_blas_threads=1) for r in protocol}
 (ROOT/'data/executed-protocol.json').write_text(json.dumps(dict(attempts=protocol,runtime_profiles=profiles,report_sha256=sha(__file__)),indent=2)+'\n')
 plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10});figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
 fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
 for ax,mode in zip(axes,['one','two']):
  names=['extend72-'+mode,'fine96-'+mode,('extend120-two' if mode=='two' and 'extend120-two' in by else 'fine120-'+mode),('shared144-one' if mode=='one' and 'shared144-one' in by else ('extend144-two' if mode=='two' and 'extend144-two' in by else 'fine144-'+mode)),('shared168-long-two' if mode=='two' and 'shared168-long-two' in by else ('extend168-one' if mode=='one' and 'extend168-one' in by else 'shared168-'+mode))]
  for name in names:
   if name not in by:continue
   r,z,t,y,d=series(name);line,=ax.plot(t,[v['gap_arclength_rms'] for v in d],label=f"{r['parameters']['angles']} points"+(' · common input' if r['audit'].get('split_one_way') else '')+(' (stopped)' if r['failure'] else ''));ax.plot(t[-1],d[-1]['gap_arclength_rms'],'x' if r['failure'] else 'o',color=line.get_color())
  ax.set(title='One-way Lorenz drive' if mode=='one' else 'Reciprocal Lorenz drive',xlabel='Model time',ylabel='Surface-weighted RMS gap (model units)');ax.legend()
 fig.suptitle('Strongest forcing: actual reached times');fig.savefig(figdir/'gap-histories.png',dpi=160);plt.close(fig)
 fig,axes=plt.subplots(2,1,figsize=(11,9),layout='constrained')
 for ax,mode in zip(axes,['one','two']):
  comps=[c for c in a['comparisons'] if c['a'].endswith('-'+mode)];labels=[]
  for c in comps:
   tag='common input' if c['a'].startswith('shared') else ('indexed' if c['a'].startswith('indexed') else 'ordinary')
   labels.append(f"{c['coarse_mesh'][0]} → {c['fine_mesh'][0]} · {tag} · [0, {c['end']:g}]")
  weighted=[100*c['gap_arclength_rms']['maximum_scaled'] for c in comps];original=[100*c['surface_gap_rms']['maximum_scaled'] for c in comps];y=np.arange(len(comps))
  ax.barh(y-.18,weighted,height=.36,label='Surface-length weighted');ax.barh(y+.18,original,height=.36,label='Original statistic');maximum=max([*weighted,*original],default=0);ax.set_xlim(0,max(.01,maximum*1.3));ax.set(yticks=y,yticklabels=labels,xlabel='Maximum scaled gap disagreement (%) · target 5%',title='One-way drive' if mode=='one' else 'Reciprocal drive');ax.invert_yaxis()
  if maximum>=5:ax.axvline(5,color='.4',ls='--',lw=1,label='5% target')
  for positions,values in [(y-.18,weighted),(y+.18,original)]:
   for pos,value in zip(positions,values):ax.annotate(f'{value:.4g}%',(value,pos),xytext=(4,0),textcoords='offset points',ha='left',va='center',fontsize=9)
  ax.legend(loc='lower right')
 fig.suptitle('Mesh comparisons: actual common recorded intervals')
 fig.savefig(figdir/'mesh-errors.png',dpi=160);plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 for ax,key,title in zip(axes,['motion_change_percent','z_std_change_percent'],['Material motion','Lorenz Z fluctuations']):
  ax.barh([f"{p['mesh']} points [0, {p['end']:g}]" for p in a['pairs']],[p[key] for p in a['pairs']],color='#007f88');ax.axvline(0,color='.5',lw=.8);ax.set(title=title,xlabel='Change with feedback (%)')
 fig.savefig(figdir/'reciprocal.png',dpi=160);plt.close(fig)
 available=[n for n in [72,96,120,144] if ('extend72-one' if n==72 else f'fine{n}-one') in by and ('extend72-two' if n==72 else f'fine{n}-two') in by]
 if 'shared168-one' in by and 'shared168-two' in by:available.append(168)
 wanted=sorted(set([min(available),max(available)]),reverse=True) if available else []
 cases=[];records={}
 for n in wanted:
  one,two=('extend72-one','extend72-two') if n==72 else (('extend168-one' if 'extend168-one' in by else 'shared168-one','shared168-two') if n==168 else (f'fine{n}-one',f'fine{n}-two'));left,right=series(one),series(two);common=np.intersect1d(left[2].round(9),right[2].round(9));count=3 if n>144 else (5 if n>120 else 7);times=common[np.unique(np.round(np.linspace(0,len(common)-1,min(count,len(common)))).astype(int))]
  label=f'Strongest forcing · {n} × {left[0]["parameters"]["layers"]}'
  if left[0]['failure'] or right[0]['failure']:label+=f' · shared time {times[-1]:g} · incomplete'
  for name,(r,z,t,y,d) in [(one,left),(two,right)]:
   m=Mesh(Parameters(**r['parameters']));frames=[]
   for clock in times:
    idx=int(np.argmin(abs(t-clock)));assert abs(t[idx]-clock)<1e-8;x,rest,l=m.unpack(y[idx]);v=d[idx]
    frames.append(dict(time=float(t[idx]),positions=x.round(7).tolist(),lorenz=l.round(5).tolist(),imprint=v['imprint'],overlap=v['penetration'],clearance=v['certified_clearance'],selfClearance=v['self_clearance'],wallClearance=v['wall_clearance'],port=float(y[idx,-1]),area=v['min_area_ratio'],quality=v['min_deformation_quality'],altitude=v['min_altitude_ratio'],gap=v['gap_arclength_rms']))
   ids=np.unique(np.round(np.linspace(0,len(t)-1,min(151,len(t)))).astype(int));records[name]=dict(frames=frames,triangles=m.tri.tolist(),initial=m.initial.round(7).tolist(),angles=m.p.angles,layers=m.p.layers,outline=z['outline'][::2].round(7).tolist(),trace=np.column_stack([t[ids],y[ids,-7],y[ids,-5]]).round(5).tolist())
  cases.append(dict(key='ellipse' if n==wanted[0] else f'mesh{n}',label=label,one=one,two=two,initialIndex=len(times)-1))
 template=(ROOT/'molding-template.html').read_text(encoding='utf-8').replace('contact-moldable-hug','resolved-material-hug').replace('contact-hug-v3','resolved-hug-v1').replace('The hug under stronger forcing','The hug: finest-mesh stress').replace('data.initial.slice(0,data.angles)','records[k].initial.slice(0,records[k].angles)').replace('data.triangles','records[k].triangles')
 template=template.replace("['Imprint',f=>f.imprint.toFixed(3)]","['Surface-weighted gap',f=>f.gap.toFixed(4)],['Smallest cell height / original',f=>(100*f.altitude).toFixed(1)+'%'],['Imprint',f=>f.imprint.toFixed(3)]")
 template=template.replace("const w=svg.parentElement.clientWidth;","svg.previousElementSibling.textContent=(k?'Two-way':'One-way')+' · '+records[k].angles+' × '+records[k].layers;const w=svg.parentElement.clientWidth;")
 template=template.replace('Recorded material motion · model units · guarded barrier contact','Recorded shared interval · 2-D model units')
 for field in ['clearance','selfClearance','wallClearance']:template=template.replace(f'f.{field}.toFixed(4)',f'(f.{field}>0&&f.{field}<.0001?f.{field}.toExponential(2):f.{field}.toFixed(4))')
 fragment=template.replace('__MOLDING_DATA__',json.dumps(dict(cases=cases,records=records),separators=(',',':')));assert len(fragment.encode())<1000000
 (ROOT/'molding-fragment.html').write_text(fragment,encoding='utf-8')
 if args.inline_output:args.inline_output.write_text(fragment,encoding='utf-8')
 table='<table class="table"><tr><th>Attempt</th><th>Reached / requested</th><th>Status</th></tr>'
 for r in a['runs']:table+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['name'],f"{r['end']:.9g} / {r['requested_end']:g}",r['failure']['message'] if r['failure'] else 'completed'])+'</tr>'
 table+='</table>'
 verdict=[]
 for mode,title in [('one','One-way'),('two','Reciprocal')]:
  comparisons=[c for c in a['comparisons'] if c['a'].endswith('-'+mode)]
  if not comparisons:continue
  c=max(comparisons,key=lambda c:(c['fine_mesh'][0],c['end']));status='passed' if c['full_mesh_check_passed'] else 'not passed'
  verdict.append(f"{title} {c['coarse_mesh'][0]} / {c['fine_mesh'][0]}: {100*c['gap_arclength_rms']['maximum_scaled']:.5g}% surface-weighted and {100*c['surface_gap_rms']['maximum_scaled']:.5g}% original gap disagreement over recorded time 0–{c['end']:g}; full-interval check {status}.")
 if any(not c['passed'] for c in a['time_controls']):verdict.append('An independent numerical repeat also failed its temporal accuracy target; see the retained control results.')
 summary='<p>'+html.escape(' '.join(verdict))+'</p>'
 report='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>The hug: finest-mesh stress</title><style>'+(ROOT/'report-style.css').read_text()+'</style><h1>Strongest Lorenz forcing: results</h1><p><a href="README.md">Results and reproduction</a> · <a href="analysis.json">Numerical audit</a> · <a href="checksums.json">File checksums</a></p>'+summary+fragment
 for name,title in [('gap-histories','Same-forcing gap histories'),('mesh-errors','Spatial agreement'),('reciprocal','Influence in both directions')]:report+='<h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+title+'">'
 report+='<h2>Every new integration</h2><div class="table-responsive">'+table+'</div><p>Only actual recorded states and reached intervals are shown. The 2-D material remains dimensionless and uncalibrated.</p></html>'
 (ROOT/'report.html').write_text(report,encoding='utf-8');print(json.dumps(dict(cases=len(cases),fragment_bytes=len(fragment.encode()),attempts=len(a['runs']))))
if __name__=='__main__':main()
