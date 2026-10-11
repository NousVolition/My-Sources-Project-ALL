"""Render only saved numerical results; never generate new motion for display."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,html,argparse,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import mesh_for
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inline-output',type=Path);args=ap.parse_args()
    a=json.loads((ROOT/'analysis.json').read_text());rows=a['runs'];by={r['key']:r for r in rows};figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
    current=lambda n,mode:next((key for key in [f'solver/mesh{n}-{mode}'] if key in by),None)
    def load(key):return np.load(ROOT/by[key]['relative_path'].replace('.json','.npz'))
    def detail(key):return json.loads((ROOT/by[key]['relative_path']).read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    history=json.loads((ROOT/'historical/comparison.json').read_text())
    for ax,oldfield,newfield,title in zip(axes,['gap_legacy','gap_arclength'],['surface_gap_rms','gap_arclength_rms'],['Original gap statistic','Surface-length weighted gap']):
        old=[r for r in history if r['label'].endswith('one')]
        ax.plot([r['angles'] for r in old],[r[oldfield] for r in old],'o--',label='Previous published model',color='#a96325')
        new=[by[k] for n in [24,48,72,96] if (k:=current(n,'one')) and by[k]['end']>=.2-1e-10]
        ax.plot([r['parameters']['angles'] for r in new],[r['final'][newfield] for r in new],'o-',label='Revised material',color='#007f88')
        ax.set(title=title,xlabel='Angular mesh points',ylabel='RMS radial gap (model units)',xticks=[24,48,72,96]);ax.legend(fontsize=9)
    fig.suptitle('One-way strong drive, matched model time 0.2');fig.savefig(figdir/'mesh-gap.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,mode in zip(axes,['one','two']):
        for n in [24,48,72,96]:
            key=current(n,mode)
            if not key:continue
            raw=load(key);d=detail(key)['diagnostics'];ax.plot(raw['time'],[v['gap_arclength_rms'] for v in d],label=f'{n} points')
        ax.set(title='One-way' if mode=='one' else 'Reciprocal',xlabel='Model time',ylabel='Surface-weighted RMS gap');ax.legend()
    fig.savefig(figdir/'gap-history.png',dpi=160);plt.close(fig)
    checks=a['solver_comparisons'];fig,ax=plt.subplots(figsize=(10,max(4,.36*len(checks))),layout='constrained')
    for i,c in enumerate(checks):
        ax.scatter(max(1e-16,c['scaled_state_error']),i,color='#007f88' if c['passed'] else '#b15826',s=38)
    ax.axvline(1e-4,color='.4',ls='--',label='0.01% target');ax.set_xscale('log');ax.set(xlabel='Largest scaled state difference (0 plotted at 1e-16)',yticks=range(len(checks)),yticklabels=[c['a']+' / '+c['b']+(' (shared partial interval)' if not c['full_requested_interval'] else '') for c in checks]);ax.legend();fig.savefig(figdir/'accuracy.png',dpi=160);plt.close(fig)
    selected=[r for r in rows if r['key'].startswith('solver/stress') and r['failure'] is None]
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    for mode,style in [('one','o-'),('two','s--')]:
        values=sorted([r for r in selected if r['name'].endswith(mode)],key=lambda r:r['parameters']['alpha'])
        for ax,field,title in zip(axes,['minimum_area','minimum_distortion','minimum_altitude'],['Area / reference area','Deformation quality','Altitude / reference altitude']):
            ax.plot([r['parameters']['alpha'] for r in values],[r[field] for r in values],style,label=mode+'-way');ax.set(title=title,xlabel='Drive amplitude alpha',ylabel='Smallest certified bound',ylim=(0,1.05));ax.legend()
    fig.suptitle('Joint stress sweep: stronger coupling, higher Lorenz rho, lower viscosity')
    fig.savefig(figdir/'stress-geometry.png',dpi=160);plt.close(fig)
    pairs=[p for p in a['pairs'] if p['name'].startswith('solver/')];fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    for ax,key,title in zip(axes,['motion_change_percent','z_std_change_percent'],['Material motion','Lorenz Z fluctuations']):
        ax.barh([p['name'].removeprefix('solver/')+f" [{p['start']:g}, {p['end']:g}]" for p in pairs],[p[key] for p in pairs],color='#007f88');ax.axvline(0,color='.4',lw=.8);ax.set(title=title,xlabel='Two-way change relative to one-way (%)')
    fig.savefig(figdir/'reciprocal.png',dpi=160);plt.close(fig)
    cases=[];records={}
    wanted=[(f'mesh{n}',f'{n}-point mesh · '+('coarse preview' if n==24 else 'strong drive · time 0.2'),current(n,'one'),current(n,'two')) for n in [72,48,24]]
    wanted += [(f'stress{n}',f'24-point stress test · drive amplitude {n}',f'solver/stress{n}-one',f'solver/stress{n}-two') for n in [16,32,64,128]]
    wanted += [('window0','Sustained drive · six time units','solver/window0-one','solver/window0-two')]
    for key,label,one,two in wanted:
        if not one or not two or one not in by or two not in by:continue
        if by[one]['failure'] or by[two]['failure']:continue
        aa,bb=load(one),load(two);common=np.intersect1d(aa['time'].round(9),bb['time'].round(9));count=11 if by[one]['parameters']['angles']>=48 else 21
        sampled=common[np.unique(np.round(np.linspace(0,len(common)-1,min(count,len(common)))).astype(int))]
        for name,raw in [(one,aa),(two,bb)]:
            r=by[name];m=mesh_for(r['parameters']);d=detail(name)['diagnostics'];frames=[]
            for t in sampled:
                idx=int(np.argmin(abs(raw['time']-t)));assert abs(raw['time'][idx]-t)<1e-8
                state=raw['state'][idx];x,rest,l=m.unpack(state);v=d[idx]
                frames.append(dict(time=float(raw['time'][idx]),positions=x.round(5).tolist(),lorenz=l.round(5).tolist(),imprint=v['imprint'],overlap=v['penetration'],clearance=v['certified_clearance'],selfClearance=v['self_clearance'],wallClearance=v['wall_clearance'],port=float(state[-1]),area=v['min_area_ratio'],quality=v['min_deformation_quality'],altitude=v['min_altitude_ratio'],gap=v['gap_arclength_rms']))
            ids=np.unique(np.round(np.linspace(0,len(raw['time'])-1,min(401,len(raw['time'])))).astype(int))
            records[name]=dict(frames=frames,triangles=m.tri.tolist(),initial=m.initial.round(5).tolist(),angles=m.p.angles,outline=raw['outline'][::2].round(5).tolist(),trace=np.column_stack([raw['time'][ids],raw['state'][ids,-7],raw['state'][ids,-5]]).round(5).tolist())
        cases.append(dict(key='ellipse' if not cases else key,label=label,one=one,two=two,initialIndex=len(sampled)-1))
    assert cases
    payload=dict(cases=cases,records=records)
    template=(ROOT/'molding-template.html').read_text(encoding='utf-8').replace('contact-moldable-hug','mesh-accurate-hug').replace('contact-hug-v3','mesh-hug-v5').replace('The hug under stronger forcing','The hug: mesh and drive tests').replace('data.initial.slice(0,data.angles)','records[k].initial.slice(0,records[k].angles)').replace('data.triangles','records[k].triangles')
    template=template.replace("['Imprint',f=>f.imprint.toFixed(3)]","['Surface-weighted gap',f=>f.gap.toFixed(4)],['Smallest cell height / original',f=>(100*f.altitude).toFixed(1)+'%'],['Imprint',f=>f.imprint.toFixed(3)]")
    template=template.replace('Recorded material motion · model units · guarded barrier contact','Recorded 2-D material motion · model units')
    fragment=template.replace('__MOLDING_DATA__',json.dumps(payload,separators=(',',':')))
    (ROOT/'molding-fragment.html').write_text(fragment,encoding='utf-8')
    if args.inline_output:
        assert len(fragment.encode())<1000000;args.inline_output.write_text(fragment,encoding='utf-8')
    css=':root{--background:#fff;--foreground:#202c33;--muted:#eef0f1;--muted-foreground:#5a6872;--border:#d0d6dc;--viz-series-1:#007f88;--viz-series-2:#a96325}body{font:16px/1.55 system-ui;color:var(--foreground);max-width:1040px;margin:30px auto;padding:0 16px;background:var(--background)}a{color:var(--viz-series-1)}button,input,select{font:inherit}select{max-width:100%}.viz-controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap}.form-label{display:block}.form-select{display:block}.form-range{display:block;width:100%}.btn{padding:8px 14px}.table{width:100%;border-collapse:collapse}.table th,.table td{text-align:left;padding:6px;border-bottom:1px solid var(--border)}.table-responsive{overflow:auto}img{max-width:100%;height:auto}pre{white-space:pre-wrap}.tabular-nums{font-variant-numeric:tabular-nums}@media(prefers-color-scheme:dark){:root{--background:#192127;--foreground:#e6edf1;--muted:#29353c;--muted-foreground:#b2c1cc;--border:#52616b;--viz-series-1:#53c3c6;--viz-series-2:#f3b16e}}'
    table='<table class="table"><thead><tr>'+''.join('<th>'+s+'</th>' for s in ['Run','Reached / requested','Result','Minimum height','Energy error'])+'</tr></thead><tbody>'
    for r in rows:
        values=[r['key'],f"{r['end']:.5g} / {r['requested_end']:g}",'completed' if r['failure'] is None else r['failure']['message'],f"{r['minimum_altitude']:.4g}",f"{max(r['energy_errors'].values()):.3g}"]
        table+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in values)+'</tr>'
    table+='</tbody></table>'
    report='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>The hug: mesh and accuracy</title><style>'+css+'</style><h1>The hug: mesh and accuracy</h1><p><a href="README.md">Results, model and reproduction</a> · <a href="analysis.json">Numerical comparisons</a> · <a href="checksums.json">File hashes</a></p><p>The coarse 24-point mesh is a preview. Quantitative gap claims require a matched refinement check. These are computed motions of an abstract 2-D material, not measurements of real clay.</p>'+fragment
    for name,title in [('mesh-gap','Gap and mesh resolution'),('gap-history','Gap at every saved time'),('accuracy','Numerical agreement, including failed earlier settings'),('stress-geometry','Cell protection under stronger forcing'),('reciprocal','Effects in both directions')]:report+='<h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+title+'">'
    report+='<details><summary>Every saved attempt, including repair pilots and incomplete runs</summary><div class="table-responsive">'+table+'</div></details><p>Finite mesh agreement is not an exact-error bound. Contact sampling, integration, geometry and energy are checked separately. Longer chaotic trajectories and new parameter choices need their own checks.</p></html>'
    (ROOT/'report.html').write_text(report,encoding='utf-8');print(json.dumps(dict(cases=len(cases),bytes=len(fragment.encode()))))
if __name__=='__main__':main()
