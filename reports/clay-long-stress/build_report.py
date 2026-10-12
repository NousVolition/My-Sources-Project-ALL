"""Generate the report and playback from actual recorded states only."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,argparse,html,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import load,Mesh,Parameters,sha
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inline-output',type=Path);args=ap.parse_args()
    a=json.loads((ROOT/'analysis.json').read_text());by={r['name']:r for r in a['runs']};figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
    def series(name):
        p=ROOT/by[name]['path'] if name in by else ROOT/'parents'/(name+'.json')
        return load(p)
    pct=lambda x:f'{100*x:.4g}%'
    text=[f"**{a['completed']} of {a['attempts']} new integrations completed** their requested interval. Parent checkpoints are inherited evidence and are not counted as new runs.",'','### Coarse-start resolution','']
    for mode,result in a['resolution'].items():
        text.append(f"- `{mode}`: accepted = {result['accepted']}; selected mesh = {result['selected_mesh']}. The selected record is `{result['selected']}`.")
    early=a['resolution']['one']['checks'][0]['errors'];last=a['resolution']['one']['checks'][-1]['errors']
    text+=['',f"The 24 × 3 to 48 × 5 gap disagreement remains {pct(early['gap_arclength_rms']['maximum_scaled'])} for surface-length weighting ({pct(early['surface_gap_rms']['maximum_scaled'])} for the original statistic). The workflow therefore rejects that pair and returns the 72 × 7 solution after the 48 × 5 / 72 × 7 check: maximum disagreements {pct(last['gap_arclength_rms']['maximum_scaled'])} and {pct(last['surface_gap_rms']['maximum_scaled'])}, respectively. These are mesh agreements, not measured physical errors. The native fixed 24-point solution has not been made accurate.",'', '| Resolution probe | Final surface-weighted gap | Difference from 72 × 7 at time 0.2 |','|---|---:|---:|']
    for r in a['axis_probes']:text.append(f"| {r['name']} | {r['gap_arclength_rms']['coarse_final']:.8g} | {pct(r['gap_arclength_rms']['final_scaled'])} |")
    text+=['','### Requested and completed long tests','','| Attempt | New interval starts | Reached / requested | Status |','|---|---:|---:|---|']
    for r in a['runs']:
        if r['name'].startswith('long'):text.append(f"| {r['name']} | {r['start']:g} | {r['end']:.9g} / {r['requested_end']:g} | {r['failure']['message'] if r['failure'] else 'completed'} |")
    text+=['','| Fine-mesh comparison | Shared recorded end | Surface-weighted maximum difference | Original-statistic maximum difference | Full requested mesh check passed |','|---|---:|---:|---:|---|']
    for c in a['comparisons']:
        if c['a'].startswith('long48') and c['b'].startswith('long72'):
            text.append(f"| {c['a']} / {c['b']} | {c['end']:g} | {pct(c['gap_arclength_rms']['maximum_scaled'])} | {pct(c['surface_gap_rms']['maximum_scaled'])} | {c['full_mesh_check_passed']} |")
    text+=['','A partial fine run cannot validate the whole requested interval, even if its shared samples agree. The long 24-point trajectories test solver behavior and reciprocal sensitivity; they remain coarse and are not a substitute for a passed fine-mesh comparison.','', '### Independent integration controls', '', '| Control pair | Independently recomputed interval | Scaled state difference | Control target met |','|---|---:|---:|---|']
    for c in a['time_controls']:
        text.append(f"| {c['a']} / {c['b']} | [{c['start']:g}, {c['end']:g}] | {pct(c['state_scaled'])} | {c['passed']} |")
    text+=['','The state target is 0.01%, normalized componentwise by max(1, peak absolute reference component) and excluding the four work bookkeeping variables. Shared parent samples are excluded from these independent-control figures. A matched checkpoint is not counted as independent agreement.','', '### Effects in both directions', '', '| Mesh and shared window | Material motion change with feedback | Lorenz Z fluctuation change with feedback |','|---|---:|---:|']
    for p in a['pairs']:text.append(f"| {p['mesh']} points [{p['start']:g}, {p['end']:g}] | {p['motion_change_percent']:+.4g}% | {p['z_std_change_percent']:+.4g}% |")
    text+=['','Material motion RMS comes from integrated viscous dissipation, and Lorenz fluctuation is the time-weighted standard deviation of saved Z values. Initial transients, mesh error and the chosen reciprocal coupling affect these numbers. They do not prove chaos suppression or better physical clay behavior.','', '### Geometry, energy and limits', '']
    if a['runs']:
        text.append(f"All saved histories passed the fresh geometry audit: **{a['all_geometry_passed']}**. All energy budgets met their target: **{a['all_energy_passed']}**. Minimum certified area ratio: **{min(r['minimum_area'] for r in a['runs']):.6g}**; deformation quality: **{min(r['minimum_quality'] for r in a['runs']):.6g}**; relative cell altitude: **{min(r['minimum_altitude'] for r in a['runs']):.6g}**. Maximum scaled full-history energy residual: **{max(max(r['energy_errors'].values()) for r in a['runs']):.6g}**.")
    for r in a['runs']:
        if r['failure']:text.append(f"- `{r['name']}` stopped at {r['end']:.9g}: {r['failure']['message']}. Its states and failure record are retained.")
    text+=['','No crossing was detected in accepted recorded states or certified straight paths. This is not a guarantee about implicit dense-output curves, arbitrary geometries, or untested forcing. The underlying material remains a research model, with no measured-clay calibration or physical uncertainty distribution.','', '![Angular and radial resolution](figures/resolution.png)','', '![Long forcing histories](figures/long-gap.png)']
    readme=(ROOT/'README.md').read_text(encoding='utf-8');lo=readme.index('<!-- RESULTS -->')+len('<!-- RESULTS -->');hi=readme.index('<!-- END RESULTS -->')
    (ROOT/'README.md').write_text(readme[:lo]+'\n\n'+'\n'.join(text)+'\n\n'+readme[hi:],encoding='utf-8')
    protocol=[json.loads((ROOT/r['path']).read_text()) for r in a['runs']]
    protocol=[{k:r[k] for k in ['name','parameters','initial','requested_end','start_time','parent','audit','source_hashes']} for r in protocol]
    (ROOT/'data/executed-protocol.json').write_text(json.dumps(dict(attempts=protocol,report_sha256=sha(__file__)),indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained');labels=[];values=[]
    for name in ['mesh24-one','axis24x5','axis24x7','axis48x3','mesh48-one','mesh72-one']:
        if name.startswith('axis') and name not in by:continue
        r,z,t,y,d=series(name);labels.append(f"{r['parameters']['angles']} × {r['parameters']['layers']}");values.append(d[-1]['gap_arclength_rms'])
    ax.bar(labels,values,color='#007f88');ax.set(xlabel='Angular points × radial layers',ylabel='Surface-weighted RMS gap (model units)',title='Same parameter values and forcing, model time 0.2')
    for i,v in enumerate(values):ax.text(i,v+.0015,f'{v:.4f}',ha='center',fontsize=9)
    ax.set_ylim(0,max(values)*1.17);fig.savefig(figdir/'resolution.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for ax,mode in zip(axes,['one','two']):
        for n in [24,48,72]:
            name=f'long{n}-{mode}'
            if name not in by:continue
            r,z,t,y,d=series(name);line,=ax.plot(t,[v['gap_arclength_rms'] for v in d],label=f'{n} points'+(' (stopped)' if r['failure'] else ''))
            ax.plot(t[-1],d[-1]['gap_arclength_rms'],'x' if r['failure'] else 'o',color=line.get_color())
        ax.set(title='One-way drive' if mode=='one' else 'Reciprocal drive',xlabel='Model time',ylabel='Surface-weighted RMS gap');ax.legend()
    fig.suptitle('Hardest forcing: only times actually reached are plotted');fig.savefig(figdir/'long-gap.png',dpi=160);plt.close(fig)
    pairs=[p for p in a['pairs'] if p['start']>0];fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,key,title in zip(axes,['motion_change_percent','z_std_change_percent'],['Material motion','Lorenz Z fluctuations']):
        ax.barh([f"{p['mesh']} points [{p['start']:g}, {p['end']:g}]" for p in pairs],[p[key] for p in pairs],color='#007f88');ax.axvline(0,color='.5',lw=.8);ax.set(title=title,xlabel='Change with feedback (%)')
    fig.savefig(figdir/'reciprocal.png',dpi=160);plt.close(fig)
    cases=[];records={}
    wanted=[('ellipse','Coarse start → accepted 72 × 7 / 48 × 5','mesh72-one','mesh48-two')]
    wanted += [(f'long{n}',f'Hardest forcing · {n}-point mesh',f'long{n}-one',f'long{n}-two') for n in [24,48,72] if f'long{n}-one' in by and f'long{n}-two' in by]
    for key,label,one,two in wanted:
        left,right=series(one),series(two);common=np.intersect1d(left[2].round(9),right[2].round(9));count=11 if max(left[0]['parameters']['angles'],right[0]['parameters']['angles'])>=48 else 21
        times=common[np.unique(np.round(np.linspace(0,len(common)-1,min(count,len(common)))).astype(int))]
        stopped=left[0]['failure'] or right[0]['failure']
        if stopped:label+=f' · shared time {times[-1]:g} · incomplete'
        for name,(r,z,t,y,d) in [(one,left),(two,right)]:
            m=Mesh(Parameters(**r['parameters']));frames=[]
            for clock in times:
                i=int(np.argmin(abs(t-clock)));assert abs(t[i]-clock)<1e-8
                x,rest,l=m.unpack(y[i]);v=d[i]
                frames.append(dict(time=float(t[i]),positions=x.round(7).tolist(),lorenz=l.round(5).tolist(),imprint=v['imprint'],overlap=v['penetration'],clearance=v['certified_clearance'],selfClearance=v['self_clearance'],wallClearance=v['wall_clearance'],port=float(y[i,-1]),area=v['min_area_ratio'],quality=v['min_deformation_quality'],altitude=v['min_altitude_ratio'],gap=v['gap_arclength_rms']))
            ids=np.unique(np.round(np.linspace(0,len(t)-1,min(401,len(t)))).astype(int))
            records[name]=dict(frames=frames,triangles=m.tri.tolist(),initial=m.initial.round(7).tolist(),angles=m.p.angles,layers=m.p.layers,outline=z['outline'][::2].round(7).tolist(),trace=np.column_stack([t[ids],y[ids,-7],y[ids,-5]]).round(5).tolist())
        cases.append(dict(key=key,label=label,one=one,two=two,initialIndex=len(times)-1))
    template=(ROOT/'molding-template.html').read_text(encoding='utf-8').replace('contact-moldable-hug','long-forced-hug').replace('contact-hug-v3','long-hug-v6').replace('The hug under stronger forcing','The hug: longer forcing').replace('data.initial.slice(0,data.angles)','records[k].initial.slice(0,records[k].angles)').replace('data.triangles','records[k].triangles')
    template=template.replace("['Imprint',f=>f.imprint.toFixed(3)]","['Surface-weighted gap',f=>f.gap.toFixed(4)],['Smallest cell height / original',f=>(100*f.altitude).toFixed(1)+'%'],['Imprint',f=>f.imprint.toFixed(3)]")
    template=template.replace("const w=svg.parentElement.clientWidth;", "svg.previousElementSibling.textContent=(k?'Two-way':'One-way')+' · '+records[k].angles+' × '+records[k].layers;const w=svg.parentElement.clientWidth;")
    template=template.replace('Recorded material motion · model units · guarded barrier contact','Recorded shared interval · 2-D model units')
    for field in ['clearance','selfClearance','wallClearance']:
        template=template.replace(f'f.{field}.toFixed(4)',f'(f.{field}>0&&f.{field}<.0001?f.{field}.toExponential(2):f.{field}.toFixed(4))')
    fragment=template.replace('__MOLDING_DATA__',json.dumps(dict(cases=cases,records=records),separators=(',',':')))
    assert len(fragment.encode())<1000000
    (ROOT/'molding-fragment.html').write_text(fragment,encoding='utf-8')
    if args.inline_output:args.inline_output.write_text(fragment,encoding='utf-8')
    table='<table class="table"><tr><th>Attempt</th><th>Reached / requested</th><th>Result</th></tr>'
    for r in a['runs']:
        table+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['name'],f"{r['end']:.7g} / {r['requested_end']:g}",r['failure']['message'] if r['failure'] else 'completed'])+'</tr>'
    table+='</table>'
    report='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>The hug: longer forcing</title><style>'+(ROOT/'report-style.css').read_text()+'</style><h1>The hug: refinement and longer forcing</h1><p><a href="README.md">Results and reproduction</a> · <a href="analysis.json">Numerical audit</a> · <a href="checksums.json">File checksums</a></p><p>The fixed 24-point error remains. Coarse-start calculations now select an actual finer mesh when its gap comparison passes. Long stress runs are shown only over times actually reached; partial runs remain incomplete.</p>'+fragment
    for name,title in [('resolution','Angular and radial refinement'),('long-gap','Long forcing histories'),('reciprocal','Reciprocal effects over the extension')]:report+='<h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+title+'">'
    report+='<h2>Every new attempt</h2><div class="table-responsive">'+table+'</div><p>Parent checkpoints and new extensions are distinguished in the report. This is an abstract material model, not measured clay.</p></html>'
    (ROOT/'report.html').write_text(report,encoding='utf-8');print(json.dumps(dict(cases=len(cases),fragment_bytes=len(fragment.encode()),attempts=len(a['runs']))))
if __name__=='__main__':main()
