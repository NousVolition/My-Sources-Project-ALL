"""Read saved study results, compare completed runs, and render the progress page."""
import os
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(__file__).resolve().parents[1]/'matplotlib-cache'))
import argparse
import hashlib
import html
import json
import math
import time
from pathlib import Path
import numpy as np
from scipy import fft
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from protocol import ROOT,CASES,GRIDS,atomic_json
from visual_data import viewer_html

OUTPUT=ROOT.parents[1]/'outputs'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
esc=lambda x:html.escape(str(x))

def table(headers,rows):
    return '<div class="table"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(c)+'</td>' for c in row)+'</tr>' for row in rows)+'</tbody></table></div>'

def difference(a,b):
    """Align by exact saved model time; normalize max error by reference peak."""
    keys=('Wmax','W_central_roi','I','energy','spin_ratio','origin_axial_strain','central_marker_parallel_stretch_spin_weighted')
    aa={round(r['t'],10):r for r in a['rows']};bb={round(r['t'],10):r for r in b['rows']}
    times=sorted(set(aa)&set(bb));out={}
    for k in keys:
        av=np.array([aa[t][k] for t in times]);bv=np.array([bb[t][k] for t in times])
        out[k]={'max_curve_difference_relative_to_reference_peak':float(np.max(abs(av-bv))/max(np.max(abs(bv)),1e-300)),
                'last_time':times[-1],'a_last':float(av[-1]),'b_last':float(bv[-1])}
    apeak=max(a['rows'],key=lambda r:r['Wmax']);bpeak=max(b['rows'],key=lambda r:r['Wmax'])
    out['largest_saved_peak_times']=[apeak['t'],bpeak['t']]
    # Samples are .02 apart: peak times are saved-sample estimates.
    return out

def field_difference(a,b):
    out=[]
    for t in (.1,.2,.3,.4):
        pa=ROOT/'runs'/a['job']['id']/f'field-t{t:.2f}.npz'
        pb=ROOT/'runs'/b['job']['id']/f'field-t{t:.2f}.npz'
        if not pa.exists() or not pb.exists():continue
        with np.load(pa,allow_pickle=False) as x,np.load(pb,allow_pickle=False) as y:
            ha=x['h'];hb=y['h'];na=a['job']['n'];nb=b['job']['n']
        if na>nb:ha,hb,na,nb=hb,ha,nb,na
        if na!=nb:
            ids=np.rint(fft.fftfreq(na)*na).astype(int)%nb
            rr=np.arange(na//2+1)
            lower=hb[:,ids[:,None,None],ids[None,:,None],rr[None,None,:]]*(na/nb)**3
        else:lower=hb
        modes=fft.fftfreq(na)*na;rz=fft.rfftfreq(na)*na
        keep=(abs(modes[:,None,None])<na/3)&(abs(modes[None,:,None])<na/3)&(rz[None,None,:]<na/3)
        lower=lower*keep;ha=ha*keep
        ua=fft.irfftn(ha,s=(na,)*3,axes=(-3,-2,-1),workers=2)
        ub=fft.irfftn(lower,s=(na,)*3,axes=(-3,-2,-1),workers=2)
        err=ua-ub
        out.append({'t':t,'lower_grid':na,'higher_grid':nb,
            'relative_velocity_L2_common_modes':float(np.sqrt(np.sum(err*err)/np.sum(ub*ub))),
            'relative_velocity_Linf_common_modes':float(np.sqrt(np.sum(err*err,axis=0)).max()/np.sqrt(np.sum(ub*ub,axis=0)).max()),
            'scope':'Physical fields compared after restriction to the coarser retained Fourier modes; fine-only modes are excluded from this error.'})
    return out

def build():
    state=read(ROOT/'task-state.json');proto=read(ROOT/'protocol.json');pre=read(ROOT/'preflight.json')
    records={p.parent.name:read(p) for p in (ROOT/'runs').glob('*/result.json')}
    completed={k:v for k,v in records.items() if v['status']=='complete'}
    cache_path=ROOT/'comparisons.json';cache=read(cache_path) if cache_path.exists() else {}
    pairs=[]
    for case in CASES:
        for method in ('fourier','fd4'):
            for n in (80,112,160):pairs.append(('timestep',f'{case}-{method}-n{n}-base',f'{case}-{method}-n{n}-half'))
            for lo,hi in zip(GRIDS,GRIDS[1:]):pairs.append(('grid',f'{case}-{method}-n{lo}-base',f'{case}-{method}-n{hi}-base'))
        for n in GRIDS:pairs.append(('method',f'{case}-fourier-n{n}-base',f'{case}-fd4-n{n}-base'))
    for method in ('fourier','fd4'):
        for n in GRIDS:
            for case in ('compressive','exodus'):
                pairs.append(('configuration',f'aligned-{method}-n{n}-base',f'{case}-{method}-n{n}-base'))
    for kind,aid,bid in pairs:
        if aid not in completed or bid not in completed:continue
        a,b=completed[aid],completed[bid]
        key=aid+'__'+bid
        fingerprint=hashlib.sha256(json.dumps([a,b],sort_keys=True).encode()).hexdigest()
        if key in cache and cache[key].get('fingerprint')==fingerprint:continue
        item={'kind':kind,'a':aid,'b':bid,'fingerprint':fingerprint,'diagnostics':difference(a,b)}
        if kind!='configuration':item['fields']=field_difference(a,b)
        cache[key]=item
    atomic_json(cache_path,cache)
    count=len(completed);nplanned=len(proto['jobs'])
    summaries=[]
    for job in proto['jobs']:
        r=records.get(job['id'])
        if r is None:summaries.append([job['case'],job['method'],job['n'],'half' if job['half'] else 'base','queued','—','—','—','—']);continue
        end=r['rows'][-1];peak=max(r['rows'],key=lambda x:x['Wmax'])
        summaries.append([job['case'],job['method'],job['n'],'half' if job['half'] else 'base',r['status'],f'{end["t"]:.2f}',f'{peak["Wmax"]:.3f} at {peak["t"]:.2f}',f'{end["Wmax"]:.3f}',f'{end["I"]:.4f}'])
    plotpath=OUTPUT/'adversarial-vortex-study.png'
    fig,axes=plt.subplots(3,3,figsize=(15,11),layout='constrained')
    for i,case in enumerate(CASES):
        for n,c in zip(GRIDS,plt.cm.viridis(np.linspace(.1,.9,len(GRIDS)))):
            for method,ls in [('fourier','-'),('fd4','--')]:
                r=records.get(f'{case}-{method}-n{n}-base')
                if r is None:continue
                rows=r['rows'];ts=[x['t'] for x in rows]
                for j,k in enumerate(('Wmax','I','central_marker_parallel_stretch_spin_weighted')):
                    axes[i,j].plot(ts,[x[k] for x in rows],ls,color=c,label=f'{n}³ {method}',lw=1.3)
        for j,label in enumerate(('Peak vorticity','Accumulated maximum spin','Central-marker weighted strain')):
            axes[i,j].set_title(f'{case}: {label}',fontsize=11);axes[i,j].set_xlim(0,.4)
            axes[i,j].set_xlabel('Model time');axes[i,j].grid(alpha=.2)
        if axes[i,0].lines:axes[i,0].legend(fontsize=7,ncol=2)
        else:axes[i,0].text(.5,.5,'Queued',ha='center',transform=axes[i,0].transAxes)
    fig.suptitle(f'New vortex study — {count}/{nplanned} runs complete; saved observations only',fontsize=14)
    fig.savefig(plotpath,dpi=150);plt.close(fig)
    comparison_rows=[]
    for c in cache.values():
        if c['kind']=='configuration':continue
        metrics=c['diagnostics']
        comparison_rows.append([c['kind'],c['a'],c['b'],f'{100*metrics["Wmax"]["max_curve_difference_relative_to_reference_peak"]:.3f}%',f'{100*metrics["I"]["max_curve_difference_relative_to_reference_peak"]:.3f}%'])
    failures=[f'{k}: {r.get("error")}' for k,r in records.items() if r['status']=='failed']
    warnings=[]
    for k,r in completed.items():
        flagged=[x for x in r['rows'] if x['resolution_flags']]
        if flagged:warnings.append([k,f'{flagged[0]["t"]:.2f}',f'{min(x["global_width"]["minimum_chord_cells"] for x in r["rows"]):.2f}',f'{100*max(x["high_band_enstrophy_fraction"] for x in r["rows"]):.2f}%'])
    initial_table=[]
    for r in pre['initial_rows']:
        if r['method']=='fourier' and r['n']==48:
            initial_table.append([r['case'],f'{r["diagnostics"]["Wmax"]:.2f}',f'{r["geometry"]["parallel_stretch"]:.3f}',f'{r["diagnostics"]["core_width"]["minimum_chord_cells"]:.2f}',f'{min(r["geometry"]["outer_mean_outward_speeds"]):.3f}'])
    final_note=('All planned runs have completed. Numerical agreement must be read from the comparison tables.' if state['status']=='complete' else 'The matrix is still running; queued cases and refinements have no results yet.')
    if failures:final_note='The matrix stopped for a numerical failure. Saved evidence is retained.'
    viewer=viewer_html(records,proto)
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Vortex stress tests</title>
    <style>body{{font:16px/1.6 system-ui,sans-serif;background:#f4f6f4;color:#203132;margin:0}}main{{max-width:1200px;margin:auto;padding:30px}}section{{background:white;padding:24px;border-radius:12px;margin:22px 0}}h1,h2{{line-height:1.25}}.label{{color:#37695e;font-weight:600}}.table{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #d9e4df}}th{{background:#eef4f1}}a{{color:#146152}}code{{overflow-wrap:anywhere}}img{{max-width:100%;height:auto}}.warning{{border-left:5px solid #b66534}}@media(max-width:600px){{main{{padding:14px}}section{{padding:16px}}}}</style></head><body><main>
    <p><a href="change-influence-result.html">Imported flow report</a></p><p class="label">New controlled comparison · design 2</p>
    <h1>One core. Three surroundings.</h1><p><strong>{count} of {nplanned} runs complete.</strong> {esc(final_note)}</p>
    {viewer}
    <section><details><summary style="font-weight:600;font-size:20px;cursor:pointer">Starting fields and settings</summary><p>These are new starting fields for the requested comparison. A central localized vortex is paired with an aligned, compressive, or initially departing background. All start at maximum vorticity 80. The cube has side 6, viscosity 0.001, zero mean velocity, and zero external force. The endpoint is model time 0.40.</p>
    {table(['Case at 48³','Initial W','Origin background strain','Core width in cells','Minimum mean outward speed'],initial_table)}
    <p>The initial strain sign and outward motion are verified. Sustained strain, later dispersal, and decay are outcomes to measure. The reverse case is the compressive variant; bending is not imposed.</p>
    <p>Version 1 was rejected because reversing its surroundings only translated the whole field. Its short pilot is archived and excluded. Version 2 localizes the central tube axially. Its initial origin spin rates are {pre['distinct_initial_interactions']['aligned_origin_domega_z_dt']:.3f} and {pre['distinct_initial_interactions']['compressive_origin_domega_z_dt']:.3f}, respectively.</p>
    <p><a href="../work/adversarial-vortex-study/protocol.json">All settings and scope</a> · <a href="../work/adversarial-vortex-study/preflight.json">Initial and solver checks</a> · <a href="../work/adversarial-vortex-study/task-state.json">Task state</a></p></details></section>
    <section><h2>Recorded curves</h2><img src="adversarial-vortex-study.png" alt="Saved maximum-vorticity, integral and strain curves for three cases; unfinished cases are marked queued."><p>Solid lines: Fourier evolution. Dashed lines: fourth-order finite-difference evolution. Both use SSP RK3. The finite-difference pressure solve and shared filter use FFTs. Thus spatial transport is replicated, while the pressure backend and time integrator are shared.</p></section>
    <section><details><summary style="font-weight:600;font-size:20px;cursor:pointer">Run matrix · {count} / {nplanned} complete</summary>{table(['Case','Method','Grid','dt setting','Status','Saved t','Largest saved W','Last W','Last I'],summaries)}<p>Half-step repeats are included at 80³, 112³ and 160³ for both methods. The cases execute in the requested order. Peak timing here uses samples separated by 0.02.</p></details></section>
    <section><details><summary style="font-weight:600;font-size:20px;cursor:pointer">Completed numerical comparisons</summary>{table(['Control','Run A','Reference B','Max W curve difference','Max I curve difference'],comparison_rows)}<p>Differences are maximum absolute curve differences divided by the reference curve maximum. Full-field comparisons use the coarser retained Fourier modes and are recorded separately.</p><p><a href="../work/adversarial-vortex-study/comparisons.json">Comparison data and full-field differences</a></p></details></section>
    <section class="warning"><h2>Resolution and validity</h2>{table(['Completed run','First resolution flag at','Smallest global chord, cells','Largest high-band enstrophy'],warnings)}<p>Fewer than six cells across a measured chord or more than 1% of enstrophy in the highest retained band flags a resolution concern. These screening levels do not prove convergence or failure of the continuous equation.</p><p>{esc('; '.join(failures) if failures else 'No hard numerical failure has been recorded.')}</p></section>
    <section><h2>Measurement details</h2><ul><li>W: maximum magnitude of Fourier curl on grid points; native FD curl is recorded separately.</li><li>I: maximum spin integrated at every RK stage. Its finite numerical value is not a continuum bound.</li><li>Core width: twelve transverse half-peak chords through a sampled peak, in domain units and cells. Interpolation adds no resolution.</li><li>Central measurements use a fixed radius-one region and material markers. Neither is an automatic vortex identity tracker.</li><li>Exodus distances compare advected centreline markers using periodic distance. Viscosity can move vorticity relative to those material points. Separation and weakened central strain must both be observed before claiming exodus.</li><li>Energy and native-enstrophy balance residuals, fixed threshold volumes, mean momentum, divergence, eigenalignment and spin ratios are in each saved result.</li></ul><p><a href="../work/adversarial-vortex-study/runs/">Saved run folders</a></p></section></main></body></html>'''
    temp=OUTPUT/'adversarial-vortex-study.tmp.html';temp.write_text(page,encoding='utf-8');temp.replace(OUTPUT/'adversarial-vortex-study.html')
    atomic_json(ROOT/'report-status.json',{'completed':count,'planned':nplanned,'study_status':state['status'],
       'built_at_unix':time.time(),'comparisons':len(cache),'hard_failures':failures,
       'report':str(OUTPUT/'adversarial-vortex-study.html')})
    return state['status'],count

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');a=p.parse_args()
    while True:
        status,count=build();print(f'report: {count}/48 complete; {status}',flush=True)
        if not a.watch or status in ('complete','failed','stopped_for_numerical_failure'):break
        time.sleep(30)
