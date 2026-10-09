"""Build the report from measurements.json; no simulation or private documents."""
import os, json, html
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get(row,key):
    for part in key.split('.'):row=row[part]
    return row
def difference(a,b,key,end=.4):
    aa={round(r['t'],8):r for r in a if r['t']<=end+1e-9};bb={round(r['t'],8):r for r in b if r['t']<=end+1e-9}
    tt=sorted(aa.keys()&bb.keys()); av=np.array([get(aa[t],key) for t in tt]);bv=np.array([get(bb[t],key) for t in tt])
    return float(np.max(abs(av-bv))/max(np.max(abs(bv)),1e-30)*100)
def chart_style():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold',
        'figure.facecolor':'#f7f9fc','axes.facecolor':'white','grid.alpha':.2,'axes.grid':True,'savefig.facecolor':'#f7f9fc'})

def main():
    d=json.loads((ROOT/'measurements.json').read_text())
    r=d['original']['base']['rows'];sur=d['surroundings'];local=d['local_terms']
    comparisons=[]
    for case in ('aligned','compressive'):
        for category,a,b in [('grid',f'{case}-fourier-n112-base',f'{case}-fourier-n160-base'),
                             ('time step',f'{case}-fourier-n112-base',f'{case}-fourier-n112-half'),
                             ('method',f'{case}-fd4-n160-base',f'{case}-fourier-n160-base')]:
            if a not in sur or b not in sur:continue
            for end in (.1,.2,.4):
                row={'case':case,'comparison':category,'a':a,'b':b,'end':end}
                for key in ('Wmax','central_marker_Wmax','central_marker_parallel_stretch_spin_weighted','core_width.minimum_half_peak_chord'):
                    row[key]=difference(sur[a]['rows'],sur[b]['rows'],key,end)
                comparisons.append(row)
    summary={'new_simulations':0,'verified_fields':len(d['verification']),'complete_result_records':len(sur),
        'original_global_peak_in_central_count':sum(x['peak_in_central_region'] for x in r),
        'original_global_peak_in_join_band_count':sum(x['peak_in_join_band'] for x in r),
        'original_outputs':len(r),'original_largest_saved_peak':max(r,key=lambda x:x['Wmax']),
        'original_final':r[-1], 'original_half_step_curve_difference_percent':difference(r,d['original']['half']['rows'],'Wmax'),
        'comparisons':comparisons,'surroundings':{}}
    for case in ('aligned','compressive'):
        rows=sur[f'{case}-fourier-n160-base']['rows'];terms=local[f'{case}-fourier-n160-base']
        take=lambda r:{k:r[k] for k in ('t','Wmax','central_marker_Wmax','central_marker_parallel_stretch_spin_weighted','W_central_roi','core_width','global_width','high_band_enstrophy_fraction')}
        summary['surroundings'][case]={'initial':take(rows[0]),'final':take(rows[-1]),
            'lowest_saved_marker_peak':take(min(rows,key=lambda x:x['central_marker_Wmax'])),
            'final_winning_marker':{'index':terms[-1]['marker_max_index'],
                'initial_position':terms[0]['marker_locations'][terms[-1]['marker_max_index']],
                'initial_spin':terms[0]['marker_spin'][terms[-1]['marker_max_index']],
                'final_spin':terms[-1]['marker_max_spin'],
                'original_center_final_spin':terms[-1]['initial_origin_marker_spin']},
            'terms':[{k:v for k,v in x.items() if not isinstance(v,list)} for x in terms]}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    chart_style()
    fig,ax=plt.subplots(2,1,figsize=(10.8,7),sharex=True,layout='constrained')
    t=[x['t'] for x in r]
    for key,label,c,ls in [('Wmax','Whole box','#a44330','-'),('central_W','Fixed interior cylinder','#176b83','-'),('join_band_W','Band near the wrapping faces','#9a761c',':')]:
        ax[0].plot(t,[x[key] for x in r],label=label,color=c,ls=ls,lw=2)
    ax[0].plot(t,[x['Wmax'] for x in d['original']['half']['rows']],color='#a44330',ls='--',lw=1,label='Whole box, half time step')
    ax[0].set(ylabel='Maximum spin (1 / model time)',title='The whole-box peak never lies in the fixed interior cylinder')
    ax[0].legend(fontsize=9,ncol=2)
    peak=summary['original_largest_saved_peak'];ax[0].annotate(f"{peak['Wmax']:.1f} at t={peak['t']:.2f}\ninterior: {peak['central_W']:.1f}",(peak['t'],peak['Wmax']),xytext=(.09,530),arrowprops={'arrowstyle':'->','color':'#495566'})
    ax[1].plot(t,[x['width_at_global_peak']['minimum_chord_cells'] for x in r],color='#7255a0',lw=2,label='Width at whole-box peak')
    ax[1].axhline(6,color='#596779',ls='--',label='Six-cell screen')
    ax[1].set(xlabel='Model time',ylabel='Width (grid cells)',ylim=(0,6.7),title='Half-step agreement does not repair the narrow spatial peak')
    ax[1].legend(fontsize=9)
    fig.suptitle('Original matched-stretch start | 64³ | unforced | viscosity 0.001',fontsize=14)
    fig.savefig(ROOT/'original-timeline.png',dpi=155);plt.close(fig)

    fig,axes=plt.subplots(3,2,figsize=(13,10.8),sharex=True,layout='constrained')
    for col,case in enumerate(('aligned','compressive')):
        rows=sur[f'{case}-fourier-n160-base']['rows'];terms=local[f'{case}-fourier-n160-base'];ts=[x['t'] for x in rows];tt=[x['t'] for x in terms]
        for key,label,c,ls in [('Wmax','Whole box','#a44330','-'),('central_marker_Wmax','Moving central points','#176b83','-'),('W_central_roi','Fixed central cylinder','#85929c',':')]:
            axes[0,col].plot(ts,[x[key] for x in rows],color=c,lw=2,ls=ls,label=label)
        axes[0,col].plot(tt,[x['initial_origin_marker_spin'] for x in terms],color='#78559b',marker='s',ls='--',label='Point initially at center')
        axes[0,col].set(title='Initially '+case,ylabel='Maximum spin (1 / time)')
        axes[0,col].legend(fontsize=8)
        axes[1,col].plot(ts,[x['central_marker_parallel_stretch_spin_weighted'] for x in rows],color='#176b83',label='Stretching',lw=2)
        axes[1,col].plot(tt,[x['weighted_viscous'] for x in terms],color='#c78517',marker='o',label='Viscous contribution',lw=1.5)
        axes[1,col].plot(tt,[x['weighted_discrete_interpolation_residual'] for x in terms],color='#85828e',marker='x',ls=':',label='Discrete/interpolation residual')
        axes[1,col].axhline(0,color='#667',lw=.8)
        axes[1,col].set(ylabel='Weighted local rate (1 / time)',title='The stretching response changes sign')
        axes[1,col].legend(fontsize=8)
        axes[2,col].plot(tt,[x['half_local_spin_width']['minimum'] for x in terms],color='#176b83',marker='o',lw=2,label='At strongest moving point')
        axes[2,col].plot(ts,[x['core_width']['minimum_half_peak_chord'] for x in rows],color='#85929c',ls=':',label='At fixed-cylinder maximum')
        axes[2,col].axhline(6*6/160,color='#a44330',ls='--',label='Six cells at this grid')
        axes[2,col].set(xlabel='Model time',ylabel='Half-local-spin width (length)',title='Widths use a changing relative threshold')
        axes[2,col].legend(fontsize=8)
    fig.suptitle('Separate periodic surroundings study | 160³ | same central start, opposite backgrounds',fontsize=14)
    fig.savefig(ROOT/'central-timeline.png',dpi=155);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(12,7.2),sharex=True,layout='constrained')
    for col,case in enumerate(('aligned','compressive')):
        for n,color in ((80,'#a07a27'),(112,'#398196'),(160,'#53398a')):
            rows=sur[f'{case}-fourier-n{n}-base']['rows'];ts=[x['t'] for x in rows]
            axes[0,col].plot(ts,[x['central_marker_Wmax'] for x in rows],color=color,label=f'{n}³',lw=2)
            axes[1,col].plot(ts,[x['central_marker_parallel_stretch_spin_weighted'] for x in rows],color=color,lw=2)
        half=sur[f'{case}-fourier-n112-half']['rows']
        axes[0,col].plot([x['t'] for x in half],[x['central_marker_Wmax'] for x in half],'--',color='#398196',label='112³ half step')
        axes[1,col].plot([x['t'] for x in half],[x['central_marker_parallel_stretch_spin_weighted'] for x in half],'--',color='#398196')
        axes[0,col].set(title='Initially '+case,ylabel='Maximum spin at moving points')
        axes[0,col].legend(fontsize=9,ncol=2)
        axes[1,col].set(xlabel='Model time',ylabel='Weighted stretching (1 / time)')
        axes[1,col].axhline(0,color='#667',lw=.8)
    fig.suptitle('Early response agrees more closely; later stretching depends on the grid',fontsize=14)
    fig.savefig(ROOT/'refinement.png',dpi=155);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),sharex=True,sharey=True,layout='constrained')
    for ax,case in zip(axes,('aligned','compressive')):
        terms=local[f'{case}-fourier-n160-base'];tt=[x['t'] for x in terms];winner=terms[-1]['marker_max_index']
        for idx,label,color in ((32,'Point initially at center','#78559b'),(winner,'Point strongest at the end','#176b83')):
            ax.plot(tt,[x['marker_spin'][idx] for x in terms],marker='o',color=color,lw=2,label=label)
        rows=sur[f'{case}-fourier-n160-base']['rows']
        ax.plot([x['t'] for x in rows],[x['central_marker_Wmax'] for x in rows],color='#87919c',ls=':',label='Maximum over all 64 points')
        ax.set(title='Initially '+case,xlabel='Model time',ylabel='Spin at a fixed material label')
        ax.legend(fontsize=8)
    fig.suptitle('The rising maximum changes which point it describes | 160³ saved run',fontsize=14)
    fig.savefig(ROOT/'same-point.png',dpi=155);plt.close(fig)
    build_report(d,summary)

def build_report(d,s):
    table=[]
    for case in ('aligned','compressive'):
        rows=d['surroundings'][f'{case}-fourier-n160-base']['rows']
        for t in (0,.1,.2,.3,.4):
            r=min(rows,key=lambda x:abs(x['t']-t))
            table.append(f"| {case} | {t:.2f} | {r['central_marker_Wmax']:.2f} | {r['Wmax']:.2f} | {r['central_marker_parallel_stretch_spin_weighted']:.3f} |")
    compare=[]
    for r in s['comparisons']:
        if r['end'] not in (.1,.4):continue
        compare.append(f"| {r['case']} | {r['comparison']} | {r['end']:.2f} | {r['Wmax']:.3f}% | {r['central_marker_Wmax']:.3f}% | {r['central_marker_parallel_stretch_spin_weighted']:.3f}% |")
    text=f'''# What happens to the central tube?

**The initial surroundings do not keep their initial stretching sign. Both tested configurations first lose spin along the moving central points, then regain it. The later response is sensitive to spatial resolution.**

This is an analysis of saved unforced runs. It adds no force, boundary rule, new initial field or time evolution. The original matched-stretch start and the separate periodic surroundings study are shown separately. New analysis is authorized for publication; no private source documents or their extracted contents are included.

## 1. Original matched-stretch start

![Original timeline](original-timeline.png)

At all **41 saved times**, the sampled whole-box maximum lies outside the fixed interior cylinder. At the largest saved peak, time **0.24**, whole-box spin is **701.79**, while the interior maximum is **213.18**. The global maximum is in the defined band near the wrapping faces at **34 of 41** times. The largest peak itself is outside that band, at `(0, -2.53125, -0.65625)`; therefore the band count does not establish its cause.

The interior cylinder is `x²+y² <= 0.4², |z| <= 2.5`. The wrapping-face band is `max(|x|,|y|,|z|) >= 2.625` in a side-6 box. These are fixed spatial regions. They do not follow material or establish the identity of the original tube after motion.

The base and half-step whole-box curves differ by at most **{s['original_half_step_curve_difference_percent']:.3f}% of the reference curve's maximum**. Spatial widths remain below six cells. These findings locate the measured peak; they do not establish resolved concentration or attribute every later peak to the join. The [earlier tube-removal control](../join-isolation/README.md) covers only time 0 to 0.04.

## 2. Separate surroundings comparison

![Moving central points, rates and widths](central-timeline.png)

These runs use a different, smooth periodic starting construction, the same central tube in each case, viscosity 0.001 and zero external force. Initially aligned and initially compressive backgrounds differ by sign. The curves use the completed 160³ Fourier runs through time 0.40.

| Initial surroundings | Time | Maximum at moving central points | Whole-box maximum | Weighted stretching |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(table)}

The two labels describe initialization, not an enforced behavior. In the aligned case, weighted stretching is positive at the start, negative by 0.10, and positive again by 0.30. In the compressive case it is negative at the start and positive by 0.20. Both sets of central points lose spin early and later regain it. The aligned run has a small initial increase before its decline. The snapshots do not establish a recurring cycle or a permanent restoring response.

### Does the same point strengthen?

![Fixed material labels](same-point.png)

In the aligned 160³ run, the point initially at the center falls from **80 to 1.90** by time 0.40. The strongest tracked point at the end is a different point: it started at `(0, 0, -2.71875)` with spin **13.43** and ends at **100.16**. In the compressive run, the initially central point ends at **18.11**, while a point starting at `(0, 0, -3)` rises from **12.91 to 103.45**.

The 64 labels span the full initial periodic axis, including its initially weaker portions. The displayed maximum therefore changes which material point it describes. This is a recorded identity change, not evidence that one original high-spin center continuously sharpened. The late numbers are specific to this grid and do not survive all refinement checks.

The width around the strongest sampled moving point initially broadens and later narrows. This width uses half that point's current spin as its threshold; the selected point can change. It is not a tracked material boundary or a fixed-threshold volume.

### Stretching and viscosity

The rates are evaluated at the same 64 moving points. For each point, stretching is `omega · ((omega · grad)u) / |omega|²`, and the viscous contribution is `omega · (nu Laplacian omega) / |omega|²`. Their displayed averages use weights proportional to `|omega|²`. A negative stretching rate means the local velocity gradient reduces the vorticity magnitude there. Viscosity is not assumed to have a negative contribution at every point.

Saved fields permit these additional term measurements every 0.10; existing spin and strain records are spaced by 0.02. The gray residual reports the difference between the sampled discrete fluid RHS plus material advection and stretching plus viscosity. It includes finite spectral truncation and interpolation effects. It becomes material later and must not be interpreted as an extra physical force. Weighted averages of rates are not derivatives of the maximum-spin curve.

## 3. What survives refinement?

![Refinement comparison](refinement.png)

Each number below is the largest absolute curve difference divided by the maximum absolute value of the reference curve over that interval. For sign-changing stretching this is a curve-scale error, not a pointwise relative error. Grid compares 112³ with 160³; time step compares 112³ base with half step; method compares FD4 with Fourier at 160³.

| Initial surroundings | Comparison | Through time | Whole-box spin | Moving-point spin | Weighted stretching |
| --- | --- | ---: | ---: | ---: | ---: |
{chr(10).join(compare)}

Early agreement is closer than late agreement. The late increase and its stretching balance cannot yet be called spatially converged. Half-step agreement alone is insufficient. The finite-difference comparison shares the pressure backend, filter and time integrator with the Fourier implementation.

## What this establishes

- The original whole-box maximum is a different measurement from spin in a fixed interior region.
- In the separate surroundings study, local stretching changes sign during unforced evolution. An initially compressive background does not guarantee continued suppression.
- The observed early broadening and later narrowing provide a response to follow, but these records do not establish repeated breathing, persistent containment or smoothness for all time.
- The departing-surroundings case was still queued at the analysis snapshot; no result is inferred for it.

## Definitions and verification

Spin means the magnitude of vorticity, `|curl u|`, in inverse model-time units. Length and time are model units, not an SI calibration. Moving points were tracked passively during the original simulations. Viscosity can separate vorticity peaks from these points; the centreline samples do not cover the whole core. The fixed cylinder used by the surroundings study has radius 1 and includes the full periodic z direction, unlike the interior cylinder used above.

This analysis checked **{s['verified_fields']} saved fields**: 82 original fields and 30 surroundings fields. Stored peaks and energy were independently reproduced for the original fields. Surroundings fields reproduced stored global spin, moving-point spin and weighted stretching. Field and result hashes, source hashes, per-point rates, locations and measurements are included in [measurements.json](measurements.json). All values passed finite and consistency checks. These checks verify this analysis, not continuum spatial resolution.

Files: [summary](summary.json), [analysis code](analyze.py), [chart and report code](report.py), [source measurements](measurements.json).

To rebuild charts from the published measurements, install NumPy and Matplotlib, then run `python report.py` in this folder. To recompute from the raw fields, `analyze.py` requires the original workspace layout, saved simulation fields, NumPy and SciPy. It reads solver definitions but never calls a time-stepping routine. Raw binary fields are retained locally.
'''
    (ROOT/'README.md').write_text(text,encoding='utf-8')
    body=''.join(f'<section><h2>{html.escape(title)}</h2><img src="{name}" alt="{html.escape(title)}"></section>' for title,name in [('Does the same point strengthen?','same-point.png'),('Original matched-stretch start','original-timeline.png'),('Moving central points and their response','central-timeline.png'),('Resolution and time-step checks','refinement.png')])
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>What happens to the central tube?</title><style>body{font:17px/1.6 system-ui;background:#f1f5f9;color:#172536;max-width:1200px;margin:auto;padding:32px}h1{font-size:36px;line-height:1.2}section{background:white;border:1px solid #dce4ed;border-radius:14px;padding:20px;margin:24px 0}img{width:100%;height:auto}p{max-width:900px}a{color:#176b83}.tag{color:#526274;font-size:14px}</style><p class="tag">SAVED DATA · ZERO EXTERNAL FORCE · TWO DISTINCT STARTS</p><h1>What happens to the central tube?</h1><p><strong>The central points lose spin, then regain it. Their stretching response changes sign.</strong> This occurs in the separate surroundings study. Its later response still depends on the grid.</p><p>In the original matched-stretch run, the whole-box maximum lies outside the fixed interior cylinder at every saved time. A fixed region does not follow the moving tube.</p>'''+body+'''<p>Thirty surroundings fields and 82 original fields checked. No new simulation. <a href="README.md">Full definitions and comparison tables</a> · <a href="summary.json">Summary data</a></p></html>'''
    (ROOT/'index.html').write_text(page,encoding='utf-8')
    print(json.dumps({'report':str(ROOT/'README.md'),'charts':4,'comparisons':len(s['comparisons'])}))

if __name__=='__main__':main()
