"""Growth audit of one unchanged unforced start, using saved base and half-step fields."""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(ROOT/'source'))
from numerics import Flow

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def growth(t,w):
    q=np.diff(w)/np.diff(t)/w[:-1]**2
    a=np.diff(np.log(w))/np.diff(t)
    j=int(q.argmax())
    return dict(max_interval_square_rate=float(q[j]),max_interval=[float(t[j]),float(t[j+1])],
      square_rates=q.tolist(),log_rates=a.tolist(),declining_intervals=int(sum(np.diff(w)<0)),
      positive_log_rate_range=[float(a[a>0].min()),float(a[a>0].max())])

def fit(t,w):
    # Same predeclared rule for both records: first uninterrupted rising part;
    # train on first 60% of samples, evaluate remaining samples in W units.
    stop=int(np.flatnonzero(np.diff(w)<=0)[0])+1
    count=int(np.ceil(.6*stop));tt=t[:count];ww=w[:count]
    a=float(least_squares(lambda z:(w[0]*np.exp(z[0]*tt)-ww)/w[0],[5.],bounds=(0,np.inf)).x[0])
    upper=(1-1e-8)/(w[0]*tt[-1])
    b=float(least_squares(lambda z:(w[0]/(1-z[0]*w[0]*tt)-ww)/w[0],[.02],bounds=(0,upper)).x[0])
    ht=t[count:stop];hw=w[count:stop];pole=1/(b*w[0])
    qvalid=bool(np.all(ht<pole))
    return dict(training_end=float(tt[-1]),segment_end=float(t[stop-1]),held_out_samples=len(ht),a=a,b=b,
      linear_error=float(np.linalg.norm(w[0]*np.exp(a*ht)-hw)/np.linalg.norm(hw)),quadratic_pole=pole,
      quadratic_error=float(np.linalg.norm(w[0]/(1-b*w[0]*ht)-hw)/np.linalg.norm(hw)) if qvalid else None,
      quadratic_valid_at_all_held_out_times=qvalid)

def main():
    f=Flow(64,.001,workers=1);cases={};manifest=[]
    for label in ('base','half'):
        p=ROOT.parent/'runs'/f'baseline-n64-{label}'
        full=json.loads((p/'result.json').read_text());rows=[]
        assert full['status']=='complete' and full['job']['nu']==.001
        for i,row in enumerate(full['series']):
            if row['t']>.35+1e-12:break
            path=p/f'field-{i:03d}.npy';digest=sha(path);h=np.load(path)
            assert np.isfinite(h).all()
            # NumPy inverse transform independently checks the stored W and energy.
            u=np.fft.irfftn(h,s=(64,)*3,axes=(-3,-2,-1))
            wh=f.curl(h);w=np.fft.irfftn(wh,s=(64,)*3,axes=(-3,-2,-1));mag=np.sqrt(np.sum(w*w,axis=0))
            ix=np.unravel_index(mag.argmax(),mag.shape);W=float(mag[ix]);energy=float(.5*np.sum(u*u)*f.dx**3)
            assert abs(W-row['Wmax'])/W<1e-12 and abs(energy-row['energy'])/energy<1e-12
            rhs,_,_=f.rhs(h);domega=f.real(f.curl(rhs));rate=np.sum(w*domega,axis=0)/np.maximum(mag,1e-300)
            tied=mag>=W*(1-1e-12)
            # At a maximum of finitely many smooth sampled magnitudes, the right
            # derivative is the greatest derivative among exactly tied maximizers.
            # Tie tolerance is explicit; report its derivative spread.
            wrate=float(rate[tied].max());single=float(rate[ix])
            out=dict(row,semidiscrete_W_rate=wrate,rate_at_argmax=single,
              instantaneous_square_ratio=wrate/W**2,tied_maxima=int(tied.sum()),
              tied_rate_range=[float(rate[tied].min()),wrate],field_sha256=digest)
            rows.append(out);manifest.append(dict(case=label,file=path.name,sha256=digest,finite=True,W_relative_error=abs(W-row['Wmax'])/W,energy_relative_error=abs(energy-row['energy'])/energy))
        t=np.array([r['t'] for r in rows]);w=np.array([r['Wmax'] for r in rows]);instant=np.array([r['instantaneous_square_ratio'] for r in rows]);j=int(instant.argmax())
        cases[label]=dict(rows=rows,dt=full['dt'],growth=growth(t,w),fit=fit(t,w),
          max_saved_instantaneous_square_ratio=float(instant[j]),max_saved_instantaneous_time=float(t[j]),
          minimum_width_cells=min(r['width_at_global_peak']['minimum_chord_cells'] for r in rows),
          maximum_high_band_enstrophy=max(r['high_band_enstrophy_fraction'] for r in rows))
        save(ROOT/f'{label}-measurements.json',cases[label]);print(label,'checked',len(rows),'fields',flush=True)
    source=json.loads((ROOT/'inputs/matched-nostop.json').read_text());t=np.array([r['t'] for r in source['series']]);w=np.array([r['biggest'] for r in source['series']])
    quote=dict(growth=growth(t,w),fit=fit(t,w),W0=float(w[0]),Wend=float(w[-1]),peak=float(w.max()),reciprocal_initial=1/float(w[0]))
    bw=np.array([r['Wmax'] for r in cases['base']['rows']]);hw=np.array([r['Wmax'] for r in cases['half']['rows']]);
    bi=np.array([r['I'] for r in cases['base']['rows']]);hi=np.array([r['I'] for r in cases['half']['rows']]);
    result=dict(status='complete',equation='Unforced incompressible Navier-Stokes',force=0,nu=.001,L=6,n=64,
        new_fluid_runs=0,one_initial_field=True,time_interval=[0,.35],saved_fields_verified=len(manifest),quoted_table=quote,
        timestep_comparison=dict(W_curve_relative_l2=float(np.linalg.norm(bw-hw)/np.linalg.norm(hw)),max_pointwise_relative_W=float(np.max(abs(bw-hw)/hw)),I_curve_relative_l2=float(np.linalg.norm(bi-hi)/np.linalg.norm(hi))),
        case_summaries={k:{a:b for a,b in c.items() if a!='rows'} for k,c in cases.items()},
        limitation='A grid maximum and its semidiscrete RHS derivative are numerical quantities. Original join and width failures remain. Saved instantaneous rates are not a bound between records. No continuum or future-time growth law is established.')
    save(ROOT/'results.json',result);save(ROOT/'verification.json',dict(status='passed',fields=manifest,
        numerics_sha256=sha(ROOT/'source/numerics.py'),imported_table_sha256=sha(ROOT/'inputs/matched-nostop.json')))
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    colors={'base':'#28657b','half':'#c97835'}
    for label,c in cases.items():
        tt=np.array([r['t'] for r in c['rows']]);ww=np.array([r['Wmax'] for r in c['rows']]);col=colors[label]
        axs[0,0].plot(tt,ww,label=label+' timestep',color=col)
        axs[0,1].plot(tt[1:],100*np.array(c['growth']['square_rates']),color=col,label='interval average: '+label)
    c=cases['base'];tt=np.array([r['t'] for r in c['rows']]);ww=np.array([r['Wmax'] for r in c['rows']]);ft=c['fit']
    axs[0,1].scatter(tt,100*np.array([r['instantaneous_square_ratio'] for r in c['rows']]),s=13,c='#6b3e79',label='RHS at saved fields: base')
    axs[0,1].axhline(16.3941960195,color='#a43d40',ls=':',label='16.39% from imported table')
    axs[1,0].plot(tt,ww,'o-',ms=3,color=colors['base'],label='Saved unforced run')
    predt=np.linspace(0,ft['segment_end'],400)
    axs[1,0].plot(predt,ww[0]*np.exp(ft['a']*predt),label='Fitted linear rate',color='#187551')
    valid=predt<ft['quadratic_pole'];axs[1,0].plot(predt[valid],ww[0]/(1-ft['b']*ww[0]*predt[valid]),label='Fitted squared rate',color='#b95330')
    axs[1,0].axvline(ft['training_end'],color='gray',ls=':',label='End of training');axs[1,0].set_ylim(0,800)
    axs[1,1].plot(tt,[r['width_at_global_peak']['minimum_chord_cells'] for r in c['rows']],color='#683c70',label='Measured width')
    axs[1,1].axhline(6,color='#a43d40',ls='--',label='Six-cell screen')
    for ax,title,y in zip(axs.flat,['One unforced start, existing timestep control','16% is not a transferable rate ceiling','Limited growth-law comparison','Peak remains under-resolved'],['Maximum vorticity W','Rate / W² (%)','Maximum vorticity W','Half-peak width (cells)']):
        ax.set(xlabel='Model time',ylabel=y,title=title);ax.legend(fontsize=8)
    fig.suptitle('Zero external force · viscosity 0.001 · 64³ · saved interval 0–0.35',fontsize=14)
    fig.savefig(ROOT/'growth-check.png',dpi=150);plt.close(fig)
    b=cases['base'];last=b['rows'][-1];inst=100*b['max_saved_instantaneous_square_ratio'];avg=100*b['growth']['max_interval_square_rate']
    report=f'''# Unforced growth-rate check

**The quoted 16.39% interval calculation reproduces. It is not a ceiling on the instantaneous growth rate or a prediction of blowup time.**

This check uses the agreed **unforced incompressible Navier–Stokes equation**, viscosity 0.001, periodic side length 6. It reuses one 64³ starting field and its already completed timestep-halving control through time 0.35. No new fluid evolution, breathing force or imposed cycle was introduced.

![Growth and resolution checks](growth-check.png)

## Measured results

| Measurement | Original timestep | Half timestep |
|---|---:|---:|
| Starting maximum spin | {bw[0]:.6f} | {hw[0]:.6f} |
| Spin at time 0.35 | {bw[-1]:.6f} | {hw[-1]:.6f} |
| Largest saved interval rate / starting W² | {avg:.4f}% | {100*cases['half']['growth']['max_interval_square_rate']:.4f}% |
| Largest RHS rate / W² at saved fields | {inst:.4f}% | {100*cases['half']['max_saved_instantaneous_square_ratio']:.4f}% |
| Accumulated maximum spin at 0.35 | {bi[-1]:.6f} | {hi[-1]:.6f} |

The W curves differ by {100*result['timestep_comparison']['W_curve_relative_l2']:.4f}% in relative L2 norm over the saved times. The largest pointwise difference is {100*result['timestep_comparison']['max_pointwise_relative_W']:.4f}%.

The RHS measurement evaluates the unforced solver's velocity derivative, takes its curl, and projects it along vorticity at the grid maximum. It includes all terms of the numerical RHS. This is a derivative of the finite grid calculation. It is checked at saved times, not bounded between them. At maxima tied within relative 1e-12, the greatest derivative is used and the tied derivative range is recorded.

## What the earlier 16% means

The separate imported 16-sample table starts at 60.198351; this reproduced 64³ start is 59.279954. Its complete original generating settings were not supplied, so these are not claimed to be the same trajectory.

For the imported table, `ΔW / (Δt × W_start²)` reaches **{100*quote['growth']['max_interval_square_rate']:.6f}%** on time 0.0747–0.0997. This is an interval average. The statement “fastest instantaneous growth was only 16%” is unsupported by those samples.

`1/60.198351 = {quote['reciprocal_initial']:.8f}` is the pole of the hypothetical equation `W'=W²`. A proven inequality `W'≤W²` would provide a comparison bound only before that time. It would not predict a failure there. That inequality has not been established for this fluid. A smaller positive coefficient in `W'=cW²` also still permits a finite pole; a small percentage alone does not establish linear growth.

## Linear versus squared growth

Apply the same rule to both records: use the first uninterrupted sampled rise; fit the first 60% of its samples with initial W fixed, then evaluate the remaining samples. The models are `W'=aW` and `W'=bW²`, fitted in W units with one nonnegative coefficient each.

For the reproduced unforced run, training ends at {ft['training_end']:.2f}; the initial rising segment ends at {ft['segment_end']:.2f}. The linear model's error on {ft['held_out_samples']} held-out samples is **{100*ft['linear_error']:.2f}%**. The squared model has a pole at **{ft['quadratic_pole']:.6f}**, {'inside the held-out interval, so it fails to predict the remaining finite values' if not ft['quadratic_valid_at_all_held_out_times'] else 'after the held-out interval'}.

The imported table gives a linear-model held-out error of {100*quote['fit']['linear_error']:.2f}% and a squared-model pole at {quote['fit']['quadratic_pole']:.6f}. This is a limited comparison of two fixed growth laws. The rising segment was selected from the observed data; this is not cross-run validation. Neither law describes the later rises and falls as a constant positive growth rule.

## What remains unresolved

All {len(manifest)} saved fields were finite. W and energy were remeasured using NumPy inverse transforms and agree with saved records within 1e-12 relative error. The initial peak is associated with the background join; the minimum measured peak width is {b['minimum_width_cells']:.3f} cells and fails the six-cell screen. High-band enstrophy reaches {100*b['maximum_high_band_enstrophy']:.2f}%. The timestep comparison cannot repair that spatial limitation.

**Result: the rate statistic can be checked, but it does not establish a resolved linear-growth law or rule out a later singularity.** A separate join-isolation experiment is testing whether growth persists after removing the central tube. The [initial-field comparison](../periodic-start-check/README.md) documents the original join and the periodic candidate.

## Files

[Results](results.json) · [Base measurements](base-measurements.json) · [Half-step measurements](half-measurements.json) · [Field verification](verification.json) · [Analysis code](analyze.py)

The local audit reads the original immutable fields in `../runs/baseline-n64-base` and `../runs/baseline-n64-half`. Those full arrays remain local; their hashes are listed in verification.json. Published measurements and the imported table contain the numerical evidence used in this report.
'''
    (ROOT/'README.md').write_text(report,encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('case_summaries','quoted_table')},indent=2))

if __name__=='__main__':main()
