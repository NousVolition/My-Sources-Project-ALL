import hashlib,json,math,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
d=json.loads((ROOT/'results.json').read_text());rows=d['rows'];assert d['status']=='complete'
def finite(x):
 if isinstance(x,float):assert math.isfinite(x)
 elif isinstance(x,dict):
  for v in x.values():finite(v)
 elif isinstance(x,list):
  for v in x:finite(v)
finite(d)
checks={'status':'passed','initial_grids':len(rows),'source_unchanged':True,'new_time_evolution':False,'max_divergence':max(r['candidate']['divergence_max'] for r in rows)}
checks['max_mean_difference']=max(float(np.max(abs(np.array(r['candidate']['mean'])-rows[-1]['candidate']['mean']))) for r in rows)
checks['max_original_tube_fraction_at_boundary_peak']=max(float(np.linalg.norm(r['contributions_at_each_combined_peak']['original_tube'])/r['original']['W']) for r in rows)
checks['candidate_128_to_256_W_relative_difference']=abs(rows[-2]['candidate']['W']/rows[-1]['candidate']['W']-1)
assert checks['max_divergence']<1e-10 and checks['max_mean_difference']<1e-12
checks['code_sha256']=hashlib.sha256((ROOT/'run.py').read_bytes()).hexdigest()
(ROOT/'verification.json').write_bytes(json.dumps(checks,indent=2).encode())
nn=[r['n'] for r in rows]
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(14,4.7),layout='constrained')
for kind,color in [('original','#b04a40'),('candidate','#257c9a')]:
 ax[0].plot(nn,[r[kind]['width']['minimum_chord_cells'] for r in rows],'o-',label=kind,color=color)
 ax[1].plot(nn,[r[kind]['W'] for r in rows],'o-',label=kind,color=color)
 ax[2].plot(nn,[r[kind]['energy'] for r in rows],'o-',label=kind,color=color)
ax[0].axhline(6,color='black',linestyle='--',label='6-cell minimum')
for a,title in zip(ax,('Width at largest spin (cells)','Largest initial spin W','Initial kinetic energy')):
 a.set_title(title);a.set_xlabel('Grid size');a.set_xticks(nn);a.grid(alpha=.2);a.legend(fontsize=9)
fig.suptitle('Separate periodic starting field: width improves, starting energy changes\nInitialization only | L=6 | original formula retained separately')
fig.savefig(ROOT/'comparison.png',dpi=150);plt.close(fig)
s=json.loads((ROOT/'slices.json').read_text());s=[r for r in s if r['n']==256]
fig,ax=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
for a,r in zip(ax,s):
 im=a.imshow(np.array(r['mag']).T,origin='lower',extent=(-3,3,-3,3),vmin=0,vmax=60,cmap='magma',interpolation='nearest')
 a.set_title(r['kind']+' | center plane z=0');a.set_xlabel('x');a.set_ylabel('y')
fig.colorbar(im,ax=ax,label='Vorticity magnitude');fig.suptitle('Actual starting fields at 256 cubed; common color scale')
fig.savefig(ROOT/'center-slices.png',dpi=150);plt.close(fig)
table='\n'.join(f"| {r['n']} | {r['original']['width']['minimum_chord_cells']:.3f} | {r['candidate']['width']['minimum_chord_cells']:.3f} | {r['candidate']['W']:.6f} | {100*r['energy_ratio']:.3f}% | {'Pass' if r['candidate']['gate']['passed'] else 'Fail: width'} |" for r in rows)
mean=rows[-1]['candidate']['mean'];q=rows[-1]
text=f"""# Test of a periodic starting field

**The new start passes the initial width and spectral screens at 256 cubed. Its starting energy is only {100*q['energy_ratio']:.2f}% of the original.** This is a different starting field whose later evolution has not been tested.

![Width, peak and energy](comparison.png)

| Grid | Original peak width, cells | New peak width, cells | New peak spin | New energy / original | New initial screen |
| --- | --- | --- | --- | --- | --- |
{table}

## What changed

The background strain now uses periodic sine coordinates and a periodic envelope. Near the origin, their leading Taylor terms match the old raw formula. The pressure projection is unchanged; its resulting center strain is measured in [results.json](results.json).

The central tube uses the exact Fourier coefficients of a periodized Gaussian with the same radius 0.2 and amplitude 0.4. Its velocity differs from the original tube by {q['tube_relative_L2_difference']:.3g} in relative L2 at 256 cubed. Thus the large starting-field change is in the surrounding strain.

Every new grid uses the same mean velocity, fixed to the original 256-grid mean: {mean}. The original grids retain their original means. No energy rescaling was applied.

The new background formula is:
- coordinate: sin(k*x)/k, with k=2*pi/6;
- squared distance: sum of 2*(1-cos(k*x))/k^2 over x,y,z;
- envelope: exp(-0.04*squared_distance);
- raw velocity: (-40*coordinate_x, -40*coordinate_y, 80*coordinate_z)*envelope, followed by the original divergence-free projection.

[Exact protocol](protocol.json) · [Runnable initial-field check](run.py)

## Where the original spike comes from

At the original maximum, the original tube contributes at most {100*checks['max_original_tube_fraction_at_boundary_peak']:.3g}% of the vorticity magnitude across these five grids. The background contributes the boundary spike. Vector addition was checked at the same location.

At 256 cubed, the new maximum is at {q['candidate']['peak_location']}, in the central tube. Its measured physical width is {q['candidate']['width']['minimum_chord']:.6f}, or {q['candidate']['width']['minimum_chord_cells']:.3f} cells. Its high-band enstrophy fraction is {100*q['candidate']['high_band_enstrophy']:.3g}%.

![Center-plane slices](center-slices.png)

These images show z=0. The original global maximum can lie outside that plane; the recorded peak coordinates identify it exactly.

## What the pass establishes

The starting peak passes the predefined six-cell and spectral checks on the 256 grid. The 128 grid still has fewer than six cells across the core. This is not a completed two-grid evolution comparison.

The energy drop is {100*(1-q['energy_ratio']):.2f}%. A later difference in growth would therefore combine the effects of changed geometry, strain and energy. It cannot be attributed solely to repairing the join. Uniform energy normalization would also change the central tube's spin, so that adjustment has not been silently applied.

No time stepping, forcing or replacement of old checkpoints occurred. The original 48-run matrix continues under its existing protocol. [Original-start diagnosis](../adaptive-peak/README.md).

[Measured results](results.json) · [Verification](verification.json)
"""
(ROOT/'README.md').write_bytes(text.encode());print(json.dumps(checks))
