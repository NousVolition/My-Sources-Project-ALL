"""Analyze only completed join counterfactuals; preserve all original measurements."""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(ROOT/'source'))
from numerics import Flow
from run import extra,initial

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main():
    f=Flow(64,.001,workers=1);cases={};manifest=[];fields={}
    for kind in ('original','without-tube','periodic-join'):
        for half in ('base','half'):
            name=kind+'-'+half
            folder=ROOT.parent/'runs'/('baseline-n64-'+half) if kind=='original' else ROOT/'runs'/name
            r=json.loads((folder/'result.json').read_text());assert r['status']=='complete'
            rows=[]
            for i,row in enumerate(r['series'][:5]):
                path=folder/f'field-{i:03d}.npy';digest=sha(path)
                if 'field_sha256' in row:assert digest==row['field_sha256']
                h=np.load(path);assert np.isfinite(h).all()
                # Independent inverse FFT backend for saved peak and energy.
                wh=f.curl(h);w=np.fft.irfftn(wh,s=(64,)*3,axes=(-3,-2,-1));mag=np.sqrt(np.sum(w*w,axis=0))
                u=np.fft.irfftn(h,s=(64,)*3,axes=(-3,-2,-1));E=float(.5*np.sum(u*u)*f.dx**3)
                assert abs(mag.max()-row['Wmax'])/row['Wmax']<1e-12
                assert abs(E-row['energy'])/E<1e-12
                rec=dict(row);rec.update(extra(f,h));rows.append(rec)
                manifest.append(dict(case=name,t=row['t'],sha256=digest,finite=True))
                fields[name,i]=h
            cases[name]=dict(dt=r['dt'],rows=rows)
    orig=fields['original-base',0];bg=fields['without-tube-base',0];tube=orig-bg
    assert np.sqrt(f.inner(initial(f,orig,'periodic-join')-fields['periodic-join-base',0],initial(f,orig,'periodic-join')-fields['periodic-join-base',0])/f.inner(orig,orig))<1e-13
    paired=[]
    for i in range(5):
        h=fields['original-base',i];b=fields['without-tube-base',i]
        w=f.real(f.curl(h));wb=f.real(f.curl(b));m=np.sqrt(np.sum(w*w,axis=0));idx=np.unravel_index(m.argmax(),m.shape)
        wb_at=wb[(slice(None),)+idx]
        paired.append(dict(t=.01*i,
            W_original=float(m[idx]),W_without_tube=float(np.sqrt(np.sum(wb*wb,axis=0)).max()),
            without_tube_vorticity_at_original_max=wb_at.tolist(),
            without_tube_magnitude_at_original_max=float(np.linalg.norm(wb_at)),
            velocity_relative_l2=float(np.sqrt(f.inner(h-b,h-b)/f.inner(h,h))),
            vorticity_relative_l2=float(np.linalg.norm(w-wb)/np.linalg.norm(w)),
            vorticity_cosine=float(np.sum(w*wb)/(np.linalg.norm(w)*np.linalg.norm(wb)))))
    controls={}
    for kind in ('original','without-tube','periodic-join'):
        a=np.array([r['Wmax'] for r in cases[kind+'-base']['rows']]);b=np.array([r['Wmax'] for r in cases[kind+'-half']['rows']])
        rows=cases[kind+'-base']['rows']
        controls[kind]=dict(W_relative_l2=float(np.linalg.norm(a-b)/np.linalg.norm(b)),max_relative_W=float(np.max(abs(a-b)/b)),
            minimum_width_cells=min(r['width_at_global_peak']['minimum_chord_cells'] for r in rows),
            maximum_high_band_enstrophy=max(r['high_band_enstrophy_fraction'] for r in rows),
            maximum_abs_energy_budget_residual=max(abs(r['energy_budget_relative_residual']) for r in rows),
            maximum_divergence=max(r['divergence_max'] for r in rows))
    o=cases['original-base']['rows'];b=cases['without-tube-base']['rows'];p=cases['periodic-join-base']['rows']
    idx=np.unravel_index(np.sqrt(np.sum(f.real(f.curl(orig))**2,axis=0)).argmax(),(64,)*3)
    tw=f.real(f.curl(tube))[(slice(None),)+idx]
    summary=dict(status='complete',new_runs=4,reused_original_runs=2,fields_verified=len(manifest),end_time=.04,n=64,nu=.001,force=0,
        original_initial_W=o[0]['Wmax'],without_tube_initial_W=b[0]['Wmax'],original_final_W=o[-1]['Wmax'],without_tube_final_W=b[-1]['Wmax'],
        without_tube_growth_percent=100*(b[-1]['Wmax']/b[0]['Wmax']-1),
        final_peak_change_on_tube_removal_percent=100*(b[-1]['Wmax']/o[-1]['Wmax']-1),
        initial_tube_magnitude_at_original_max=float(np.linalg.norm(tw)),
        repair_initial_energy_change_percent=100*(p[0]['energy']/o[0]['energy']-1),
        repair_initial_W=p[0]['Wmax'],repair_final_W=p[-1]['Wmax'],
        timestep_controls=controls,paired_fields=paired,
        limitation='Short 64-grid attribution experiment. Repair changes the whole background energy and projected strain, so it does not isolate a seam-only effect. No attribution of the later t~0.3 peak, spatial convergence, or physical concentration claim.')
    save(ROOT/'analysis.json',summary);save(ROOT/'measurements.json',cases)
    save(ROOT/'verification.json',dict(status='passed',saved_fields=manifest,numerics_sha256=sha(ROOT/'source/numerics.py'),checks='Finite fields, stored field hashes, independent inverse-FFT W and energy within1e-12 relative; exact original tube and mean retained in repaired start.'))
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    colors=['#245d76','#ca762e','#18745b']
    for kind,label,col in zip(('original','without-tube','periodic-join'),('Original','Central tube removed','Periodic background + same tube'),colors):
        for half in ('base','half'):
            r=cases[kind+'-'+half]['rows'];t=[v['t'] for v in r]
            ax[0].plot(t,[v['Wmax'] for v in r],color=col,ls='-' if half=='base' else '--',label=label if half=='base' else None)
        ax[1].plot(t,[v['interior_central_cylinder_W'] for v in cases[kind+'-base']['rows']],color=col,label=label)
    ax[0].set(title='Early peak persists without the central tube',ylabel='Global maximum vorticity W',xlabel='Model time')
    ax[1].set(title='Fixed interior central region',ylabel='Maximum vorticity in central cylinder',xlabel='Model time')
    for a in ax:a.legend(fontsize=8)
    fig.suptitle('Same unforced equation · 64³ · ν=0.001 · base and half timesteps')
    fig.savefig(ROOT/'comparison.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(11,7),layout='constrained',sharex=True,sharey=True)
    vmax=max(cases[k+'-base']['rows'][-1]['Wmax'] for k in ('original','without-tube','periodic-join'))
    for j,(kind,label) in enumerate(zip(('original','without-tube','periodic-join'),('Original','Tube removed','Periodic background'))):
        for row,i in enumerate((0,4)):
            w=f.real(f.curl(fields[kind+'-base',i]));m=np.sqrt(np.sum(w*w,axis=0));projection=m.max(axis=1)
            im=axs[row,j].imshow(projection.T,origin='lower',extent=(-3,3,-3,3),vmin=0,vmax=vmax,cmap='magma',aspect='equal')
            axs[row,j].set(title=label+f', t={i*.01:.2f}',xlabel='x',ylabel='z')
    fig.colorbar(im,ax=axs.ravel().tolist(),label='Vorticity magnitude, maximum along y',shrink=.8)
    fig.savefig(ROOT/'locations.png',dpi=150);plt.close(fig)
    report=f'''# Does the join-associated peak need the central tube?

**No, over this short numerical interval. Removing the central tube leaves the initial maximum unchanged and the peak still rises {summary['without_tube_growth_percent']:.2f}% through time 0.04.**

![Peak comparison](comparison.png)

| Start | W at0 | W at0.04 | Starting energy | Center strain ∂u_z/∂z |
|---|---:|---:|---:|---:|
| Original | {o[0]['Wmax']:.6f} | {o[-1]['Wmax']:.6f} | {o[0]['energy']:.6f} | {o[0]['center_gradient'][2][2]:.6f} |
| Tube removed | {b[0]['Wmax']:.6f} | {b[-1]['Wmax']:.6f} | {b[0]['energy']:.6f} | {b[0]['center_gradient'][2][2]:.6f} |
| Periodic background, same tube | {p[0]['Wmax']:.6f} | {p[-1]['Wmax']:.6f} | {p[0]['energy']:.6f} | {p[0]['center_gradient'][2][2]:.6f} |

At0.04, removing the tube changes the global peak by {summary['final_peak_change_on_tube_removal_percent']:.4f}%. Its vorticity magnitude at the original run's maximum is {paired[-1]['without_tube_magnitude_at_original_max']:.6f}, versus {o[-1]['Wmax']:.6f} in the original run. The short-time peak therefore does not require the central tube.

![Peak locations](locations.png)

## What changed

The original base and half-step trajectories were reused. Four new controls remove the tube or replace the raw background by a smooth periodic sine-coordinate construction, each at the same original timestep and half timestep. All use the same unforced incompressible Navier–Stokes solver, positive viscosity0.001, side6 and spectral projection. No force, energy rescaling or change to existing simulations was made.

The periodic variant keeps the original sampled tube and mean exactly. It changes starting energy by {summary['repair_initial_energy_change_percent']:.2f}% and changes the projected central strain. It is a separate starting-field experiment, not an adopted replacement and not a seam-only modification. A lower peak there cannot be attributed solely to repairing the join.

The interior central measurement uses x²+y²≤0.4² and |z|≤2.5. It is a fixed spatial region, not an advected label identifying tube material. After evolution, the difference between two runs is a nonlinear response to removing the tube, not an additive decomposition of the original evolved flow.

## Checks and limits

All {len(manifest)} saved fields were checked for finite values, hashes, peak spin and energy. The tube-removal W curves differ by {100*controls['without-tube']['W_relative_l2']:.6f}% under timestep halving; the periodic-background curves differ by {100*controls['periodic-join']['W_relative_l2']:.6f}%. Original and control peaks remain subject to the existing width and spectral screens, recorded in the data.

This establishes that the early peak persists in the background without the central tube. It does not establish a resolved physical spike, isolate all effects of the raw join, or explain the much later peak near time0.3. Spatial refinement remains necessary for those conclusions.

[Measurements](measurements.json) · [Comparisons](analysis.json) · [Verification](verification.json) · [Protocol](protocol.json) · [Runner](run.py) · [Analysis](analyze.py)

Raw field arrays remain local with published hashes. The original baseline input fields are in the existing matched-stretch study. Running this audit requires those saved arrays; the portable reconstruction helper generates equivalent reference runs in a fresh directory.
'''
    report=report.replace('W at0','W at 0').replace('At0.04','At 0.04').replace('viscosity0.001','viscosity 0.001').replace('side6','side 6').replace('time0.3','time 0.3')
    report=report.replace('Original and control peaks remain subject to the existing width and spectral screens, recorded in the data.',f"Minimum peak widths are {controls['original']['minimum_width_cells']:.3f} cells in the original, {controls['without-tube']['minimum_width_cells']:.3f} with the tube removed, and {controls['periodic-join']['minimum_width_cells']:.3f} with the periodic background. All fail the six-cell screen. Spectral and budget checks are recorded in the data.")
    report += '\nTo reproduce in a new directory, install requirements.txt and run `python reproduce.py --output NEW_DIRECTORY`. This generates equivalent original reference trajectories and then the four controls; it does not overwrite the saved study.\n'
    (ROOT/'README.md').write_text(report,encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
