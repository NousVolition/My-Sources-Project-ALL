"""Build figures and offline playback from computed states only."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from analyze_results import mesh
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';FIG=ROOT/'figures'
def main():
    a=json.loads((DATA/'analysis.json').read_text());rows=json.loads((DATA/'all-runs.json').read_text());lookup={r['name']:r for r in rows}
    FIG.mkdir(exist_ok=True);plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
    selected=[('Previous material','prior-extreme-one'),('Reinforced / same forcing','prior-extreme-new-one'),('Harder / alpha 32','harder32-one'),('Hardest / alpha 64','harder64-one')]
    # Baseline file and new same-forcing run intentionally use different names.
    selected=[(label,name) for label,name in selected if (DATA/(name+'.npz')).exists()]
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for label,name in selected:
        raw=np.load(DATA/(name+'.npz'));row=lookup.get(name)
        if row is None:
            from contact_hug import Mesh as OldMesh
            from cell_hug import Parameters,geometry
            old=json.loads((DATA/(name+'.json')).read_text());m=OldMesh(Parameters(**old['parameters']))
        else:m=mesh(row)
        from cell_hug import geometry
        geoms=[geometry(m,m.unpack(s)[0]) for s in raw['state']]
        for ax,k in zip(axes,['min_area_ratio','min_quality_ratio','min_altitude_ratio']):ax.plot(raw['time'],[g[k] for g in geoms],label=label)
    for ax,title in zip(axes,['Smallest area / original area','Smallest shape quality / original','Smallest altitude / original']):
        ax.set(xlabel='Model time',ylabel='Ratio',title=title);ax.set_ylim(bottom=0);ax.legend(fontsize=8)
    fig.savefig(FIG/'cell-quality.png',dpi=160);plt.close(fig)
    # Show the same forcing, saved time and axes for the old and new material.
    if 'prior-extreme-new-one' in lookup:
        old_raw=np.load(DATA/'prior-extreme-one.npz');new_raw=np.load(DATA/'prior-extreme-new-one.npz')
        from contact_hug import Mesh as OldMesh
        from cell_hug import Parameters,geometry
        oldrow=json.loads((DATA/'prior-extreme-one.json').read_text());oldmesh=OldMesh(Parameters(**oldrow['parameters']))
        index=int(np.argmin([geometry(oldmesh,oldmesh.unpack(s)[0])['min_area_ratio'] for s in old_raw['state']]))
        selected_time=float(old_raw['time'][index]);fig,axes=plt.subplots(1,2,figsize=(9,4.7),layout='constrained')
        for ax,raw,m,label,color in [(axes[0],old_raw,oldmesh,'Earlier thin cells','#a96325'),(axes[1],new_raw,mesh(lookup['prior-extreme-new-one']),'Reinforced material','#007f88')]:
            idx=int(np.argmin(abs(raw['time']-selected_time)));x=m.unpack(raw['state'][idx])[0];outline=raw['outline']
            ax.fill(outline[:,0],outline[:,1],color='.85');ax.plot([-2.2,2.2,2.2,-2.2,-2.2],[-2.2,-2.2,2.2,2.2,-2.2],color='.6',lw=.8)
            ax.add_collection(PolyCollection(x[m.tri],facecolors=color,edgecolors=color,alpha=.35,linewidths=.65))
            ax.set(xlim=(-2.25,2.25),ylim=(-2.25,2.25),aspect='equal',xlabel='x (model units)',ylabel='y (model units)',title=label)
        fig.suptitle(f'Matched previous extreme one-way drive, saved time {selected_time:g}');fig.savefig(FIG/'molding-comparison.png',dpi=160);plt.close(fig)
    pairs=[p for p in a['pairs'] if p['key'] in ['prior-extreme-new','harder32','harder64','amplitude4','amplitude16','amplitude32']]
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,key,title in zip(axes,['material_motion_rms','z_std'],['Material motion','Lorenz Z fluctuations']):
        ax.bar(range(len(pairs)),[p['changes_percent'][key] for p in pairs]);ax.axhline(0,color='.5',lw=.8)
        ax.set(xticks=range(len(pairs)),xticklabels=[p['key'].replace('prior-extreme-new','prior combined') for p in pairs],ylabel='Two-way change from one-way (%)',title=title);ax.tick_params(axis='x',labelrotation=30)
    fig.suptitle('Matched second-half windows; each case has its own recorded duration');fig.savefig(FIG/'reciprocal.png',dpi=160);plt.close(fig)
    if a['long_blocks']:
        fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        for seed in range(3):
            blocks=[b for b in a['long_blocks'] if b['key']=='long-s'+str(seed)]
            for ax,key in zip(axes,['material_motion_rms','z_std']):ax.plot([(b['start']+b['end'])/2 for b in blocks],[b['changes_percent'][key] for b in blocks],marker='o',label='Start '+str(seed))
        for ax,title in zip(axes,['Material motion','Lorenz Z fluctuations']):ax.axhline(0,color='.5',lw=.8);ax.set(xlabel='Center of six-unit block',ylabel='Two-way change (%)',title=title);ax.legend()
        fig.savefig(FIG/'long-windows.png',dpi=160);plt.close(fig)
    records={};cases=[]
    for key,label in [('prior-extreme-new','Same previous extreme'),('harder32','Harder drive · alpha 32'),('harder64','Hardest drive · alpha 64'),('amplitude4','Molding · alpha 4'),('long-s0','Sustained harder drive')]:
        names=[key+'-'+s for s in ['one','two']]
        if not all(n in lookup and lookup[n]['failure'] is None for n in names):continue
        for n in names:
            row=lookup[n];raw=np.load(DATA/(n+'.npz'));m=mesh(row)
            ids=np.unique(np.round(np.linspace(0,len(raw['time'])-1,41)).astype(int));frames=[]
            for idx in ids:
                state=raw['state'][idx];x,rest,l=m.unpack(state);d=m.diagnostics(state)
                frames.append(dict(time=float(raw['time'][idx]),positions=x.round(5).tolist(),lorenz=l.round(5).tolist(),imprint=d['imprint'],contact=d['contact_fraction'],overlap=d['penetration'],clearance=d['certified_clearance'],selfClearance=d['self_clearance'],wallClearance=d['wall_clearance'],port=float(state[-1]),area=d['min_area_ratio'],quality=d['min_quality_ratio'],altitude=d['min_altitude_ratio']))
            ids=np.unique(np.round(np.linspace(0,len(raw['time'])-1,301)).astype(int))
            records[n]=dict(frames=frames,outline=raw['outline'][::2].round(5).tolist(),trace=np.column_stack([raw['time'][ids],raw['state'][ids,-7],raw['state'][ids,-5]]).round(5).tolist())
        cases.append(dict(key='ellipse' if key=='prior-extreme-new' else key,label=label,one=names[0],two=names[1],initialIndex=int(np.argmin([f['area'] for f in records[names[0]]['frames']]))))
    assert cases,'No validated completed pair available'
    m=mesh(lookup[cases[0]['one']]);payload=dict(cases=cases,records=records,triangles=m.tri.tolist(),initial=m.initial.round(5).tolist(),angles=24)
    template=(ROOT/'molding-template.html').read_text(encoding='utf-8').replace('contact-moldable-hug','cell-moldable-hug').replace('contact-hug-v3','cell-hug-v4').replace('The hug under stronger forcing','The hug with cell protection').replace("['Imprint',f=>f.imprint.toFixed(3)]","['Cell area / original',f=>f.area.toFixed(3)],['Shape quality / original',f=>f.quality.toFixed(3)],['Cell altitude / original',f=>f.altitude.toFixed(3)],['Imprint',f=>f.imprint.toFixed(3)]").replace('guarded barrier contact','contact and cell-shape barriers')
    fragment=template.replace('__MOLDING_DATA__',json.dumps(payload,separators=(',',':')));assert len(fragment.encode())<1000000
    (ROOT/'molding-fragment.html').write_text(fragment,encoding='utf-8')
    inline=Path('C:/Users/ziman/.codex/visualizations/2026/10/10/01a12397-aa86-7133-8023-8ff85a4b8a8c/hug-cell-protection.html');inline.write_text(fragment,encoding='utf-8')
    table=[]
    for r in rows:
        v=r['audit'];met=r['metrics'];table.append('<tr>'+''.join('<td>'+html.escape(str(s))+'</td>' for s in [r['name'],r['completed_end'],r['requested_end'],'incomplete' if r['failure'] else 'completed',f"{v['minimum_swept_area_ratio']:.5g}",f"{v['minimum_swept_quality_ratio']:.5g}",f"{v['minimum_swept_altitude_ratio']:.5g}",f"{max(met[k] for k in ['material_budget_error','lorenz_budget_error','total_budget_error']):.3g}" if met else 'n/a'])+'</tr>')
    text=(ROOT/'README.md').read_text(encoding='utf-8')
    report='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Hug cell protection and reciprocal stress</title><style>:root{--background:#fff;--foreground:#202c33;--muted:#eef0f1;--muted-foreground:#5a6872;--border:#d0d6dc;--viz-series-1:#007f88;--viz-series-2:#a96325;--primary:#202c33;--primary-foreground:#fff}body{font:16px/1.55 system-ui;color:var(--foreground);max-width:1040px;margin:30px auto;padding:0 16px}button,input,select{font:inherit}select{max-width:100%}.viz-controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap}.form-label{display:block}.form-select{display:block}.form-range{display:block;width:100%}.btn{padding:8px 14px}.table{width:100%;border-collapse:collapse}.table th,.table td{text-align:left;padding:6px;border-bottom:1px solid var(--border)}.table-responsive{overflow:auto}img{max-width:100%;height:auto}pre{white-space:pre-wrap}.results{overflow:auto}.results table{min-width:800px}.tabular-nums{font-variant-numeric:tabular-nums}@media(prefers-color-scheme:dark){:root{--background:#192127;--foreground:#e6edf1;--muted:#29353c;--muted-foreground:#b2c1cc;--border:#52616b;--viz-series-1:#53c3c6;--viz-series-2:#f3b16e}body{background:var(--background)}}</style><h1>Cell protection and harder reciprocal forcing</h1><p><a href="README.md">Model, exact settings, results and reproduction instructions</a> · <a href="data/analysis.json">Audit and comparisons</a> · <a href="checksums.json">File hashes</a></p>'+fragment+'<h2>Cell deformation</h2><img src="figures/cell-quality.png" alt="Measured cell area, shape quality and altitude compared with the previous material"><h2>Reciprocal effects</h2><img src="figures/reciprocal.png" alt="Finite-window changes to material motion and Lorenz fluctuation">'+('<img src="figures/long-windows.png" alt="Variations across initial states and time blocks">' if a['long_blocks'] else '')+'<h2>All recorded attempts</h2><div class="results"><table class="table"><thead><tr>'+''.join('<th>'+s+'</th>' for s in ['Run','Reached time','Requested time','Status','Swept area ratio','Swept quality ratio','Swept altitude ratio','Budget error'])+'</tr></thead><tbody>'+''.join(table)+'</tbody></table></div><h2>Scope</h2><p>Dimensionless 2-D material initially surrounding a fixed object. Resistance to compression and poor cell shape is a chosen material law. Independent solver agreement tests numerical accuracy, not physical clay. Finite runs do not prove global stability or the Navier–Stokes problem.</p></html>'
    report=report.replace('<h2>Reciprocal effects</h2>','<img src="figures/molding-comparison.png" alt="Old and reinforced cells at the same forcing and saved time"><h2>Reciprocal effects</h2>')
    (ROOT/'report.html').write_text(report,encoding='utf-8');print('Playback bytes',len(fragment.encode()),'cases',len(cases))
if __name__=='__main__':main()
