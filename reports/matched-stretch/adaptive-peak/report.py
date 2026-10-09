import hashlib,json,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    d=json.loads((ROOT/'results.json').read_text());initial=d['initial'];runs=d['runs'];common=d['common_available_end']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
    nn=[r['n'] for r in initial];width=[r['width_at_global_peak']['minimum_chord_cells'] for r in initial];core=[r['width_at_central_roi_peak']['minimum_chord_cells'] for r in initial]
    axes[0].plot(nn,width,'o-',label='At largest spin');axes[0].plot(nn,core,'o-',label='In central-tube region');axes[0].axhline(6,color='#a62729',linestyle='--',label='Minimum: 6 cells');axes[0].set_ylabel('Minimum transverse half-peak width / dx');axes[0].set_title('Initial width');axes[0].legend(fontsize=8)
    axes[1].plot(nn,[100*r['high_band_enstrophy_fraction'] for r in initial],'o-',color='#c87016');axes[1].axhline(1,color='#a62729',linestyle='--');axes[1].set_ylabel('Enstrophy near cutoff (%)');axes[1].set_title('Initial spectral check: limit 1%')
    axes[2].plot(nn,[min(3-abs(x) for x in r['peak_location']) for r in initial],'o-',color='#7147a5');axes[2].set_ylabel('Distance to nearest wrapping face');axes[2].set_title('Location of largest initial spin')
    for ax in axes:ax.set_xlabel('Grid size N');ax.set_xticks(nn);ax.grid(alpha=.2)
    fig.suptitle('Adaptive test stops at the starting field\nOriginal matched-stretch formula | L=6 | viscosity=0.001 | no new time evolution')
    fig.savefig(ROOT/'initial-resolution.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7.6),layout='constrained')
    fields=[('Wmax','Peak spin W'),('logW','log W'),('inverseW','1 / W'),('I','Accumulated peak spin I'),('log_growth_rate','Saved-time d(log W)/dt'),('max_to_mean_spin','Maximum / mean spin')]
    for n in ('64','80','128'):
        r=runs[n];rows=[x for x in r['rows'] if x['t']<=common+1e-12];tt=[x['t'] for x in rows]
        for ax,(key,title) in zip(axes.flat,fields):ax.plot(tt,[x[key] for x in rows],'.-',label=f'{n} cubed');ax.set_title(title);ax.set_xlabel('Model time');ax.grid(alpha=.2)
    axes[0,0].legend();fig.suptitle(f'Existing runs over their common available interval: 0 to {common:.2f}\nResolution gates fail; no singular-time extrapolation is fitted')
    fig.savefig(ROOT/'growth-checks.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,4.9),layout='constrained')
    for n in ('64','80','128'):
        rows=runs[n]['rows'];r=min(rows,key=lambda a:abs(a['t']-.12));p=r['profile']
        axes[0].plot(p['scaled_distance'],p['scaled_magnitude'],label=f'{n}: t={r["t"]:.4f}, width={p["width_cells"]:.2f} cells')
    for target in (0,.12,.24,.4):
        r=min(runs['64']['rows'],key=lambda a:abs(a['t']-target));p=r['profile'];axes[1].plot(p['scaled_distance'],p['scaled_magnitude'],label=f't={r["t"]:.2f}; {p["width_cells"]:.2f} cells')
    for ax in axes[:2]:ax.axhline(.5,color='grey',linestyle=':');ax.set_xlabel('Distance / half-width');ax.set_ylabel('Local spin / maximum spin');ax.set_ylim(0,1.3);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0].set_title('Across grids near t=0.12');axes[1].set_title('Across times on 64 cubed')
    for n in ('64','80','128'):
        rows=[r for r in runs[n]['rows'] if r['t']<=common+1e-12];axes[2].plot([r['t'] for r in rows],[r['profile']['width_cells'] for r in rows],'.-',label=n)
    axes[2].axhline(6,color='#a62729',linestyle='--');axes[2].set_ylabel('Minimum core width / dx');axes[2].set_xlabel('Model time');axes[2].set_title('Actual samples across the core');axes[2].legend();axes[2].grid(alpha=.2)
    fig.suptitle('Profiles rescaled by measured peak and transverse half-width\nTrilinear interpolation draws the lines; it does not increase grid resolution')
    fig.savefig(ROOT/'peak-profiles.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,7.2),layout='constrained')
    rows=[r for r in runs['64']['rows'] if r['t']<=common+1e-12];tt=[r['t'] for r in rows]
    for name,sign,label in [('stretching_RHS',1,'Stretching'),('transport_LHS',-1,'Minus transport'),('viscosity_RHS',1,'Viscosity')]:axes[0,0].plot(tt,[sign*r['peak_budget'][name]['along_vorticity'] for r in rows],'.-',label=label)
    axes[0,0].set_title('Signed terms at the same current peak');axes[0,0].set_ylabel('Rate of local spin magnitude');axes[0,0].legend()
    for n in ('64','80','128'):
        r=[v for v in runs[n]['rows'] if v['t']<=common+1e-12];t=[v['t'] for v in r]
        axes[0,1].plot(t,[v['alignment_degrees'] if v['alignment_degrees'] is not None else np.nan for v in r],'.-',label=n)
        axes[1,1].plot(t,[v['width_at_global_peak']['minimum_chord'] for v in r],'.-',label=n)
    axes[0,1].set_title('Vorticity angle to most stretching strain axis');axes[0,1].set_ylabel('Degrees; 0 = aligned');axes[0,1].legend()
    for axis,label in enumerate(('x','y','z')):axes[1,0].plot(tt,[r['peak_location'][axis] for r in rows],'.-',label=label)
    axes[1,0].set_title('64-grid global-maximum coordinates');axes[1,0].set_ylabel('Position in cube');axes[1,0].legend()
    axes[1,1].set_title('Core width in physical domain units');axes[1,1].set_ylabel('Minimum half-peak chord');axes[1,1].legend()
    for ax in axes.flat:ax.set_xlabel('Model time');ax.grid(alpha=.2)
    fig.suptitle('Colocated vorticity diagnostics | L=6, viscosity=0.001\nA fresh global maximum is selected at each output; this is not a particle trajectory')
    fig.savefig(ROOT/'local-mechanism.png',dpi=150);plt.close(fig)
    slices=json.loads((ROOT/'slices.json').read_text());fig,axes=plt.subplots(1,len(slices),figsize=(14,4.4),layout='constrained')
    for ax,s in zip(axes,slices):
        mag=np.array(s['mag']);im=ax.imshow(mag.T/mag.max(),origin='lower',extent=(-3,3,-3,3),cmap='magma',vmin=0,vmax=1);ax.set_title(f't={s["t"]:.2f}; z={s["z"]:.3f}');ax.set_xlabel('x');ax.set_ylabel('y')
    fig.colorbar(im,ax=list(axes),label='Spin / peak in this slice',shrink=.8);fig.suptitle('Actual 64-grid slices through the maximum | Original saved fields\nLater panels are unresolved coarse-run observations, not a cross-grid comparison')
    fig.savefig(ROOT/'peak-slices.png',dpi=150);plt.close(fig)
    table='\n'.join(f"| {r['n']} cubed | {r['width_at_global_peak']['minimum_chord_cells']:.3f} | {r['width_at_central_roi_peak']['minimum_chord_cells']:.3f} | {100*r['high_band_enstrophy_fraction']:.3f}% | Stop |" for r in initial)
    comp='\n'.join(f"| {k} | {100*v['W_curve_relative_error']:.3f}% | {'Pass' if v['passes_5_percent'] else 'Fail'} |" for k,v in d['comparison'].items())
    samplecounts=sum(len(r['rows']) for r in runs.values());peakcount=sum(len(r['sampled_local_maximum_indices']) for r in runs.values())
    text=f'''# Adaptive peak-spin test: original start

**The original start fails the requested resolution gate before time evolution.** The strongest-spin region is about two cells wide on every checked grid, including 256 cubed. The central-tube region gets better resolved, but it is not where the largest initial spin occurs.

![Initial resolution gate](initial-resolution.png)

| Grid | Width at largest spin, in cells | Width in central-tube region, in cells | Enstrophy near cutoff | Adaptive decision |
| --- | --- | --- | --- | --- |
{table}

The minimum width is a transverse half-peak chord, measured in 12 directions perpendicular to local vorticity. It is not the cube root of the total hot-region volume. The required minimum is 6 cells, with 8 preferred. The enstrophy-band limit is 1%; either failure stops extension. [Fixed protocol](protocol.json).

## What was tested

We retained the original formula, pressure projection, filter, mean, viscosity and unforced equation. Existing 64/128/256 initial audits were reused, 80 was remeasured from its saved initial field, and a new 112 initialization was generated. No new time evolution was started in this adaptive test.

The raw background strain has different curl traces at opposite wrapping faces. Fourier filtering makes each individual grid field smooth, but has not established one smooth periodic limiting start. The largest initial peaks lie near those faces. Their physical width decreases with grid spacing while staying near two cells. This is evidence of a grid-scale feature in the supplied construction.

The user chose to keep this original start for diagnosis. No smooth replacement or energy renormalization was introduced. Initial energies and means differ across grids and are recorded in [results.json](results.json).

## Peak growth on the common available interval

![Growth checks](growth-checks.png)

The saved 64, 80 and 128 curves overlap through time {common:.2f}. The 112 check has only time zero because its initial gate fails. Thus the complete requested ladder has no accepted evolution interval. The three existing curves are shown as **unresolved observations**, not as a resolved comparison.

| Existing grids | Relative W-curve L2 difference | 5% curve screen |
| --- | --- | --- |
{comp}

The available 64-grid half-step study already changes W by only 0.10865% through 0.4. Temporal agreement does not repair the spatial width and spectral failures. `log W`, `1/W`, `I` and the saved-time growth rate are plotted without fitting a singular time to unresolved data.

## Core profiles and location

![Rescaled profiles](peak-profiles.png)

Rescaling can make even a feature spanning only two cells look similar across grids. Such a collapse cannot override the width failure. The lines use periodic trilinear interpolation; no finer samples were measured between the original grid points.

![Local mechanism](local-mechanism.png)

The three vorticity terms are recorded at the same selected maximum. Transport is shown with its minus sign on the right-hand side. The angle uses the most stretching eigenvector of the local symmetric velocity gradient; an almost degenerate eigenvalue pair makes the angle undefined. The full discrete RHS is compared with the sum of these filtered terms.

The maximum is selected afresh at every output. It can switch locations, especially when several cells nearly tie. This is not a continuously tracked fluid parcel. The records include the number of near-tied maxima and exact peak coordinates.

![Saved flow slices](peak-slices.png)

## Coverage and remaining tests

- {samplecounts} saved fields remeasured for W, I, maximum/mean ratio, energy, enstrophy, peak position, physical width, fixed and half-peak volumes, spectra, divergence, local terms and alignment.
- {peakcount} local maxima detected at the saved output cadence. Peaks between saved times cannot be recovered; the report does not claim every within-step event was measured.
- Fixed and relative threshold volumes, shell totals and shell averages are saved together in [results.json](results.json).
- The original fixed 48-run matrix already covers three viscosities requested here, half timesteps, random divergence-free perturbations, amplitude and angle changes. Its status remains in the matched-stretch report. Its longer unresolved records are not an adaptive pass.
- The proposed longer normalized parameter matrix and geometric perturbation ensemble have not been run by this audit. They are held at the failed starting-resolution gate. Normalizing the old field would change it and was not done.

**Decision: unresolved concentration. Do not extend this original peak under the requested adaptive rule.** This diagnoses the numerical start and the sampled evolution; it does not establish a fluid singularity.

## Reproduce and inspect

The copied `source/numerics.py` has the same SHA-256 as the running study; it is not edited. `audit.py` reads immutable checkpoints and adds diagnostics. `report.py` generates these charts. `field-manifest.json` records hashes of the fields read. Full arrays stay in the local calculation workspace.

Run `python audit.py --study-root PATH_TO_MATCHED_STRETCH --legacy-root PATH_TO_FIXED_CUTOFF --max-128-time {runs['128']['rows'][-1]['t']:.2f}`, then `python report.py`. Python dependencies: NumPy, SciPy and Matplotlib. The main saved-time cadence is 0.01; the 80-grid replay uses its recorded exact times approximately 0.02 apart. Actual dt values and all selected times are stored per run.

This publication is a completed audit of saved observations. The fixed 128-grid simulation itself is still running toward 0.4; this audit freezes its available data through {runs['128']['rows'][-1]['t']:.2f}. The 64-grid record ends at 0.4 and the completed 80-grid replay ends at 0.12. Only 0 through 0.12 is used for their cross-grid curves.
'''
    (ROOT/'README.md').write_bytes(text.encode());print(json.dumps({'saved_fields':samplecounts,'sampled_local_maxima':peakcount,'common_end':common,'grid_comparisons':d['comparison']}))
if __name__=='__main__':main()
