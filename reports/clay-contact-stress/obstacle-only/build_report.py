from pathlib import Path
import base64,json,hashlib,html,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from contact_hug import *
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';FIG=ROOT/'figures'

def main():
    analysis=json.loads((DATA/'analysis.json').read_text());rows=json.loads((DATA/'all-runs.json').read_text());lookup={r['name']:r for r in rows}
    assert analysis['all_audits_passed'] and analysis['all_independent_checks_passed']
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    # Numerical contact failure before/after, and actual positive certified margins.
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    old=analysis['baseline'];i=np.arange(2)
    axes[0].bar(i-.18,[r['old_overlap'] for r in old],.36,label='Old penalty');axes[0].bar(i+.18,[r['new_overlap'] for r in old],.36,label='New barrier')
    axes[0].set(xticks=i,xticklabels=['One-way','Two-way'],ylabel='Maximum radial penetration (model units)',title='Matched former failure: alpha = 2')
    for j,r in enumerate(old):axes[0].text(j+.18,.001,'0',ha='center');axes[0].text(j-.18,r['old_overlap']+.001,f"{r['old_overlap']:.4f}",ha='center')
    axes[0].set_ylim(0,.033);axes[0].legend()
    chosen=['repair-seed0','amplitude4','amplitude8','amplitude16','rho160','adversarial']
    for direction,offset in [('one',-.18),('two',.18)]:
        axes[1].bar(np.arange(len(chosen))+offset,[lookup[k+'-'+direction]['audit']['minimum_swept_clearance'] for k in chosen],.36,label=direction+'-way')
    axes[1].set(xticks=np.arange(len(chosen)),xticklabels=['a2','a4','a8','a16','rho160','combined'],ylabel='Minimum swept clearance lower bound',title='All accepted motion stays outside')
    axes[1].legend();fig.savefig(FIG/'contact.png',dpi=170);plt.close(fig)
    # Each block is a finite-window comparison, never presented as a confidence interval.
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for seed in range(3):
        blocks=[b for b in analysis['long_blocks'] if b['seed']==seed]
        for ax,key in zip(axes,['mesh_motion_rms','z_std']):
            ax.plot([b['start']+10 for b in blocks],[b['changes_percent'][key] for b in blocks],marker='o',label='Initial state '+str(seed))
    for ax,title in zip(axes,['Material motion','Lorenz Z fluctuations']):
        ax.axhline(0,color='.5',linewidth=.8);ax.set(title=title,xlabel='Center of 20-unit time block',ylabel='Two-way change from one-way (%)');ax.legend()
    fig.suptitle('Long runs: reaction can change the effect with the initial state and time window')
    fig.savefig(FIG/'reciprocal.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    comparison={c['key']:c for c in analysis['comparisons']}
    keys=['repair-seed1','amplitude4','amplitude8','rho160','fast','adversarial']
    for ax,field,title in zip(axes,['two_forward','two_return'],['Drive → hug','Hug → drive']):
        ax.bar(np.arange(len(keys)),[comparison[k][field] for k in keys]);ax.set(xticks=np.arange(len(keys)),xticklabels=['a2','a4','a8','rho160','fast','combined'],ylabel='Gross exchanged work (model units)',title=title)
    fig.suptitle('Both signs of exchange in the two-way runs; durations differ by test')
    fig.savefig(FIG/'work.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(10,7),layout='constrained')
    joint=[]
    for direction in ['one','two']:
        row=lookup['amplitude8-'+direction];raw=np.load(DATA/(row['name']+'.npz'));m=Mesh(Parameters(**row['parameters']))
        joint.append(np.array([m.diagnostics(s)['surface_gap_rms'] for s in raw['state']]))
    near_time=float(raw['time'][np.argmin(joint[0]+joint[1])])
    for r,direction in enumerate(['one','two']):
        row=lookup['amplitude8-'+direction];raw=np.load(DATA/(row['name']+'.npz'));m=Mesh(Parameters(**row['parameters']))
        for ax,time in zip(axes[r],[0,near_time,16]):
            idx=np.argmin(abs(raw['time']-time));x=m.unpack(raw['state'][idx])[0];outline=raw['outline']
            ax.fill(outline[:,0],outline[:,1],color='.85',label='Rigid object')
            ax.plot([-2.2,2.2,2.2,-2.2,-2.2],[-2.2,-2.2,2.2,2.2,-2.2],color='.6',linewidth=.8)
            ax.add_collection(PolyCollection(x[m.tri],facecolors='#238b8e' if r else '#b67f35',alpha=.3,edgecolors='#238b8e' if r else '#b67f35',linewidths=.5))
            ax.set(xlim=(-2.25,2.25),ylim=(-2.25,2.25),aspect='equal',xlabel='x (model units)',ylabel='y (model units)',title=f"{direction}-way, t={raw['time'][idx]:g}")
    fig.suptitle('Computed filled material under alpha = 8');fig.savefig(FIG/'molding.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    steady=analysis['steady_mesh'];axes[0].plot([s['angles'] for s in steady],[s['final_imprint'] for s in steady],marker='o');axes[0].set(xlabel='Angles (radial layers also increase)',ylabel='Final imprint (dimensionless)',title='Steady forming, t = 12')
    subs=[lookup['lobed-sub'+str(s)] for s in [8,16,32]]
    axes[1].plot([8,16,32],[s['metrics']['final']['surface_gap_rms'] for s in subs],marker='o');axes[1].set(xlabel='Pieces per material edge',ylabel='RMS radial surface gap',title='Lobed conservative offset, t = 8')
    for offset,key,label in [(-.18,'mesh_motion_rms','Material motion'),(.18,'z_std','Z fluctuations')]:
        axes[2].bar(np.arange(2)+offset,[comparison[c]['changes_percent'][key] for c in ['amplitude4','mesh48']],.36,label=label)
    axes[2].set(xticks=[0,1],xticklabels=['24 angles / 3 layers','48 angles / 5 layers'],ylabel='Two-way change (%)',title='Matched alpha-4 interaction');axes[2].legend()
    fig.savefig(FIG/'refinement.png',dpi=170);plt.close(fig)
    # Interactive playback contains recorded states, not an independently animated curve.
    cases=[dict(key='ellipse',label='Former strong case · alpha 2',base='repair-seed1'),dict(key='strong',label='Harder drive · alpha 8',base='amplitude8'),dict(key='rho160',label='Higher Lorenz forcing · rho 160',base='rho160'),dict(key='lobed',label='Lobed object · alpha 4',base='lobed')]
    records={}
    for case in cases:
        for direction in ['one','two']:
            name=case['base']+'-'+direction;row=lookup[name];raw=np.load(DATA/(name+'.npz'));m=Mesh(Parameters(**row['parameters']))
            indices=np.unique(np.round(np.linspace(0,len(raw['time'])-1,51)).astype(int));frames=[]
            for idx in indices:
                t=float(raw['time'][idx]);state=raw['state'][idx];x,rest,l=m.unpack(state);d=m.diagnostics(state)
                frames.append(dict(time=t,positions=np.round(x,5).tolist(),lorenz=np.round(l,5).tolist(),imprint=d['imprint'],contact=d['contact_fraction'],overlap=d['penetration'],clearance=d['certified_clearance'],wallClearance=d['wall_clearance'],port=float(state[-1])))
            ids=np.unique(np.round(np.linspace(0,len(raw['time'])-1,301)).astype(int))
            records[name]=dict(frames=frames,outline=np.round(raw['outline'][::2],5).tolist(),trace=np.column_stack([raw['time'][ids],raw['state'][ids,-7],raw['state'][ids,-5]]).round(5).tolist())
        case['one']=case['base']+'-one';case['two']=case['base']+'-two';case['initialIndex']=0;del case['base']
    first=Mesh(Parameters());payload=dict(cases=cases,records=records,triangles=first.tri.tolist(),initial=first.initial.round(5).tolist(),angles=24)
    fragment=(ROOT/'molding-template.html').read_text(encoding='utf-8').replace('__MOLDING_DATA__',json.dumps(payload,separators=(',',':')))
    assert len(fragment.encode('utf-8'))<1000000
    (ROOT/'molding-fragment.html').write_text(fragment,encoding='utf-8')
    # Results section in source documentation is generated from completed runs.
    changes=comparison['amplitude8']['changes_percent']
    solver=json.loads((DATA/'solver-checks.json').read_text())
    worst=max(c.get('error',0) for c in solver['checks'])
    block_ranges={key:[min(b['changes_percent'][key] for b in analysis['long_blocks']),max(b['changes_percent'][key] for b in analysis['long_blocks'])] for key in ['mesh_motion_rms','z_std']}
    text=f"""{analysis['configurations']} final-protocol configurations completed. All have zero measured object and wall overlap, positive swept whole-edge clearance, positive swept cell area, and passing saved-frame edge-crossing checks. The smallest swept lower bound is {analysis['minimum_swept_clearance']:.6g} model units. The largest scaled energy-balance error across material, Lorenz and combined checks is {analysis['maximum_budget_error']:.6g}, below the 0.001 target.

The matched old alpha-2 failures had maximum radial penetration {old[0]['old_overlap']:.6g} (one-way) and {old[1]['old_overlap']:.6g} (two-way). Both are zero in the upgraded contact model. This is a changed contact law, not a claim that the old solution was corrected by time refinement alone.

The stress sweep reaches alpha 16 (four times the former alpha-4 stress amplitude), rho 160, viscosity 0.15, contact strength 100000, varied reservoir capacities, three object shapes and three initial states. Matched long runs reach time 64. The most extreme combined alpha-16/rho-160/viscosity-0.15 test runs to time 2; most paired comparisons run to time 16, alpha-16 tests to time 8. All durations and settings are explicit in `data/protocol.json`.

In the alpha-8, initial-state-1 comparison over time 4–16, two-way reaction changes material motion by {changes['mesh_motion_rms']:+.3f}% and Lorenz Z fluctuations by {changes['z_std']:+.3f}%. Across nine 20-unit blocks from three longer runs, motion changes range from {block_ranges['mesh_motion_rms'][0]:+.3f}% to {block_ranges['mesh_motion_rms'][1]:+.3f}%, and Z fluctuation changes range from {block_ranges['z_std'][0]:+.3f}% to {block_ranges['z_std'][1]:+.3f}%. These are observed finite-window ranges, not confidence intervals or universal damping guarantees.

Reduced motion does not mean tighter molding. Signed Lorenz forcing can push the material outward toward the square. In the alpha-8 pair, final RMS radial surface gaps are {lookup['amplitude8-one']['metrics']['final']['surface_gap_rms']:.4f} (one-way) and {lookup['amplitude8-two']['metrics']['final']['surface_gap_rms']:.4f} (two-way), much larger than the contact activation width. The playback shows this expansion as well as the transient approach to the object.

All independent and refinement checks pass; the largest reported short state discrepancy is {worst:.6g} on the scaled state metric (target 0.0001). The zero-port one/two control is identical. Whole raw states, independent solver recordings, all per-run audits, mesh/subdivision comparisons and the retained earlier attempts are included. See `data/analysis.json`, `data/solver-checks.json` and the report for exact values.
"""
    readme=(ROOT/'README.md').read_text(encoding='utf-8');readme=re.sub(r'<!-- RESULTS -->.*?<!-- /RESULTS -->','<!-- RESULTS -->\n'+text+'<!-- /RESULTS -->',readme,flags=re.S);(ROOT/'README.md').write_text(readme,encoding='utf-8')
    # All test rows are visible, including numerical reliability and actual steps.
    table=''.join(f"<tr><td>{html.escape(r['name'])}</td><td>{r['end']:g}</td><td>{r['parameters']['alpha']:g}</td><td>{r['parameters']['rho']:g}</td><td>{r['audit']['minimum_swept_clearance']:.5g}</td><td>{max(r['metrics'][k] for k in ['material_budget_error','lorenz_budget_error','total_budget_error']):.3g}</td><td>{r['audit']['minimum_step']:.3g}</td><td>{r['audit']['rejected_trials']}</td></tr>" for r in rows)
    paired=''.join(f"<tr><td>{c['key']}</td><td>{c['changes_percent']['mesh_motion_rms']:+.2f}%</td><td>{c['changes_percent']['z_std']:+.2f}%</td><td>{c['two_forward']:.5g}</td><td>{c['two_return']:.5g}</td></tr>" for c in analysis['comparisons'])
    images=''.join(f'<figure><img alt="{label}" src="data:image/png;base64,{base64.b64encode((FIG/name).read_bytes()).decode()}"><figcaption>{label}</figcaption></figure>' for name,label in [('contact.png','Old penetration and new positive clearance bounds, in model units.'),('reciprocal.png','Separate initial states and finite time blocks; these are not confidence intervals.'),('work.png','Forward and return work are integrated on accepted steps; test durations differ.'),('molding.png','Computed material triangles, rigid object and matched forcing rules.'),('refinement.png','Refinement changes an approximation; the finest mesh is not an exact answer.')])
    title='The hug: stronger contact and reciprocal stress tests'
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title><style>
body{{margin:0;padding:24px;font:16px/1.5 system-ui;background:#fafafa;color:#263943;--foreground:#263943;--border:#cad1d5;--muted:#e4e8e8;--muted-foreground:#52666e;--viz-series-1:#197f8f;--viz-series-2:#ac7136;--primary:#263943;--primary-foreground:#fff}}main{{max-width:1000px;margin:auto}}h1,h2,h3{{font-weight:500}}h2{{margin-top:30px}}a{{color:#126a7d}}img{{width:100%;height:auto}}figure{{margin:24px 0}}table{{border-collapse:collapse;width:100%}}th,td{{padding:7px;text-align:left;border-bottom:1px solid var(--border)}}.table-responsive{{overflow-x:auto}}.viz-controls,.viz-row{{display:flex;gap:16px;flex-wrap:wrap;align-items:center}}.form-label{{display:flex;flex-direction:column;gap:4px;flex:1;min-width:0}}.form-select{{max-width:100%;padding:8px;font:inherit}}.btn{{padding:8px 12px;font:inherit}}.btn-primary{{background:var(--primary);color:var(--primary-foreground);border:0;border-radius:6px}}.form-range{{width:100%}}.text-small{{font-size:12px}}.tabular-nums{{font-variant-numeric:tabular-nums}}@media(max-width:430px){{body{{padding:12px}}.viz-controls .form-label{{flex-basis:100%}}}}
@media(prefers-color-scheme:dark){{body{{background:#192025;color:#e0e7eb;--foreground:#e0e7eb;--border:#52626d;--muted:#35434c;--muted-foreground:#b0c0ca;--viz-series-1:#69c2ce;--viz-series-2:#edb16a;--primary:#e0e7eb;--primary-foreground:#192025}}a{{color:#8cd2df}}}}
</style></head><body><main><h1>{title}</h1><p>{analysis['configurations']} completed configurations. No object or square-wall penetration in the final protocol. Smallest certified swept clearance: {analysis['minimum_swept_clearance']:.6g}. Largest scaled energy error: {analysis['maximum_budget_error']:.6g}.</p><p>The contact spring has been replaced by a barrier and guarded adaptive steps. Whole segments are checked. The material still yields and keeps a fading rest-shape imprint. Default contact introduces a small stand-off gap and changes the contact response.</p><p><a href="README.md">Equations, geometric bounds, sources, reproduction and limitations</a> · <a href="../clay-hug/report.html">Earlier material results</a> · <a href="data/protocol.json">Frozen protocol</a> · <a href="data/analysis.json">Recomputed results</a></p>
{fragment}<h2>Contact and reciprocal behavior</h2>{images}<h2>Matched influence in both directions</h2><p>Statistics exclude the first 4 model-time units where the run is long enough; shorter stress/control runs exclude their first half. Positive changes mean an increase with feedback. Work columns refer to the full two-way run. Different trajectories alone do not establish stabilization.</p><div class="table-responsive"><table><thead><tr><th>Pair</th><th>Material motion change</th><th>Z fluctuation change</th><th>Drive → hug work</th><th>Hug → drive work</th></tr></thead><tbody>{paired}</tbody></table></div>
<h2>Numerical and geometric checks</h2><p>13 analytic/gradient/domain checks, independent DOP853 integrations, step/tolerance refinement, zero-port controls, meshes with 24/48/72 angles, and lobed contact subdivisions with 8/16/32 pieces. Proper crossings of nonadjacent mesh edges were checked at every saved frame. The swept certificate applies to obstacle contact and cell area, not full continuous self-contact.</p><div class="table-responsive"><table><thead><tr><th>Configuration</th><th>End time</th><th>alpha</th><th>rho</th><th>Swept clearance</th><th>Max energy error</th><th>Smallest half step</th><th>Rejected trials</th></tr></thead><tbody>{table}</tbody></table></div>
<h2>What these results establish</h2><p>The numerical hug can react back on the Lorenz drive through an explicit work-balanced port. The repaired boundary excludes the fixed objects under the tested forcing. Effects on variability depend on initial state and time window. Refinement agreement supports numerical reliability for the tested equations.</p><p>This remains a dimensionless 2-D prototype without friction, fracture, physical clay calibration or experimental measurements. Finite runs do not prove universal stability, global existence, or any Navier–Stokes theorem. Older inaccurate and interrupted attempts remain under <a href="data/initial-attempts/protocol.json">initial attempts</a> and <a href="data/low-barrier-attempts/protocol.json">weak-barrier attempts</a>.</p></main></body></html>'''
    (ROOT/'report.html').write_text(page,encoding='utf-8')
    (DATA/'report-build.json').write_text(json.dumps(dict(builder_sha256=hashlib.sha256((ROOT/'build_report.py').read_bytes()).hexdigest(),fragment_bytes=len(fragment.encode()),selected_cases=cases,block_ranges=block_ranges),indent=2)+'\n')
    print('Built report, figures and',len(fragment.encode()),'byte playback',flush=True)
if __name__=='__main__':main()
