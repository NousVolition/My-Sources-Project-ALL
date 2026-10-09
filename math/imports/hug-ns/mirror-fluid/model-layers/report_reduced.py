import hashlib,json,math,os
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parent
d=json.loads((ROOT/'reduced-results.json').read_text())
assert d['status']=='complete'
def finite(x):
 if isinstance(x,float):assert math.isfinite(x)
 elif isinstance(x,dict):
  for v in x.values():finite(v)
 elif isinstance(x,list):
  for v in x:finite(v)
finite(d);tight={};perm=[1,0,3,2]
for name,r in d['four_node'].items():
 K=np.array(r['K']);pr=(np.sum(abs(K),axis=1)>0).astype(float);y0=r['states'][0];t=np.array(r['t'])
 sol=solve_ivp(lambda tt,z:pr*(.5*z-z**3+K@z-K.sum(axis=1)*z),(0,40),y0,t_eval=t,method='DOP853',rtol=1e-12,atol=1e-14)
 assert sol.success;tight[name]=sol.y.T
res=max(float(np.max(abs(y[:,perm]-tight[name.replace('-s1','-s-1')]))) for name,y in tight.items() if name.endswith('-s1'))
checks=dict(d['checks'],four_node_tighter_mirror_max=res,number_of_reduced_runs=sum(len(d[k]) for k in ('ball','four_node','listening','spatial')),finite=True)
assert res<1e-8 and checks['ball_tolerance_refinement_max']<1e-7
assert checks['spatial_dt_relative_error']<1e-5 and max(checks['spatial_grid_relative_errors'].values())<1e-6
checks['status']='passed';checks['source_sha256']=hashlib.sha256((ROOT/'run_reduced.py').read_bytes()).hexdigest()
(ROOT/'reduced-verification.json').write_bytes(json.dumps(checks,indent=2).encode())
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
for A in (0,.15,.5):
 r=d['ball'][f'r1-A{A}-s1'];t=np.array(r['t']);x=np.array(r['x']);v=np.array(r['v']);late=t>=60
 axes[0].plot(t,x,label=f'A={A}');axes[1].plot(x[late],v[late],label=f'A={A}');s=np.array(r['strobe']);axes[2].plot(np.arange(len(s)),s[:,0],'.-',label=f'A={A}')
for ax,title in zip(axes,('Signed lean','Late phase portrait (time 60–80)','One sample per breathing cycle')):
 ax.set_title(title);ax.grid(alpha=.2);ax.legend(fontsize=9)
axes[0].set_xlabel('Model time');axes[0].set_ylabel('x');axes[1].set_xlabel('x');axes[1].set_ylabel('v');axes[2].set_xlabel('Cycle');axes[2].set_ylabel('x at phase zero')
fig.suptitle('Driven double-well ball | r=1, damping=0.4, period=4\nReduced model with an explicitly chosen signed force')
fig.savefig(ROOT/'ball.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for delta,absent in ((0,False),(.05,False),(.2,False),(0,True)):
 r=d['four_node'][f'd{delta}-absent{absent}-s1'];label='D absent' if absent else f'delta={delta}'
 axes[0].plot(r['t'],r['q_AB'],label=label);axes[1].plot(r['t'],r['E_model'],label=label)
for ax,title in zip(axes,('A–B signed difference / 2','Unsigned two-pair difference')):
 ax.set_title(title);ax.set_xlabel('Model time');ax.grid(alpha=.2);ax.legend()
fig.suptitle('Four coupled scalar states | same initial A and B\nChanging one reciprocal coupling pair or removing D')
fig.savefig(ROOT/'four-node.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for h in (0,.1,.4):
 r=d['listening'][f'h{h}-b0.0001'];y=np.array(r['states'])
 axes[0].plot(r['t'],r['q_AB'],label=f'h={h}');axes[1].plot(r['t'],y[:,2],label=f'h={h}')
for ax,title in zip(axes,('A–B signed difference / 2','Listening variable q')):
 ax.axhline(0,color='grey',linestyle=':');ax.set_title(title);ax.set_xlabel('Model time');ax.grid(alpha=.2);ax.legend()
fig.suptitle('Two-way feedback in the reduced ODE\nInitial A–B half-difference = 0.0001; exact zero-bias controls remain zero')
fig.savefig(ROOT/'feedback.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
for ax,(A,local) in zip(axes.flat,((0,False),(.15,False),(.5,False),(.15,True))):
 r=d['spatial'][f'n256-A{A}-loc{local}-s1']
 im=ax.imshow(r['snapshots'],origin='lower',extent=(0,2*np.pi,-2,34),aspect='auto',vmin=-1,vmax=1,cmap='coolwarm',interpolation='nearest')
 ax.set_title(f'A={A}'+(' + localized initial change' if local else ''));ax.set_xlabel('Position');ax.set_ylabel('Model time (samples every 4)')
fig.colorbar(im,ax=axes,label='Signed local amplitude')
fig.suptitle('One-dimensional reaction–diffusion model | 256 points\nReflection: a(x) -> -a(-x); signed patches are scalar states, not fluid vortices')
fig.savefig(ROOT/'spatial.png',dpi=150);plt.close(fig)
ballrows='\n'.join(f"| {r} | {a} | {d['ball'][f'r{r}-A{a}-s1']['crossings']} | {d['ball'][f'r{r}-A{a}-s1']['late_crossings']} | {d['ball'][f'r{r}-A{a}-s1']['last_x']:.5f} |" for r in (-1,1) for a in (0,.15,.5))
feedbackrows='\n'.join(f"| {h} | {d['listening'][f'h{h}-b0.0001']['q_AB'][-1]:.6f} | {d['listening'][f'h{h}-b0.0001']['states'][-1][2]:.6f} |" for h in (0,.1,.4))
text=f"""# Enlarging the lean model in layers

**Four reduced-model layers are complete: {checks['number_of_reduced_runs']} runs. Seventeen separate driven fluid controls are underway.** The reduced models show their own measured behavior; a match to fluid behavior must be checked separately.

## 1. Motion and breathing

dx/dt = v; dv/dt = r*x - x^3 - 0.4*v + A*sin(2*pi*t/4).

The paired reflected run negates x, v and the applied force. Keeping a signed additive force unchanged while negating only the start would not preserve this symmetry.

![Ball time series, phase portrait and cycle samples](ball.png)

| r | A | Crossings through 80 | Crossings during 60–80 | Final x |
| --- | --- | --- | --- | --- |
{ballrows}

There is no chaos classification from these finite records. The once-per-cycle values and crossing counts are measurements. The chosen double-well force supplies its side states.

## 2. A–B–C–D coupling

The four-state model uses f(z)=0.5*z-z^3 and the reciprocal coupling matrix in each record. Reflection swaps A/B and C/D and transforms the entire coupling matrix. Both outgoing and incoming cross couplings obey that transformation.

q_AB=(A-B)/2; q_CD=(C-D)/2; E_model=sqrt(q_AB^2+q_CD^2).

![Four-node comparison](four-node.png)

The coupling changes are 0, 0.05 and 0.2. The unmatched-source case removes D and all its edges; its reflected partner removes C. An imposed coupling difference is a controlled input, not spontaneous direction selection.

## 3. Listening and feedback

dA/dt=-A-A^3+0.3*(B-A)+h*q;
dB/dt=-B-B^3+0.3*(A-B)-h*q;
dq/dt=-0.2*q+0.8*(A-B)-q^3.

![Feedback](feedback.png)

| Feedback h | Final (A-B)/2 from +0.0001 | Final q |
| --- | --- | --- |
{feedbackrows}

All exact zero-bias controls remain zero. Opposite biases give opposite responses. This feedback is in the reduced model only. No q feedback has been added to the existing Navier–Stokes studies.

## 4. Spatial lean

da/dt=0.5*a-a^3+0.05*Laplacian(a)+A*sin(2*pi*t/4)*sin(x), on a periodic line of length 2*pi.

![Spatial states](spatial.png)

Here reflection maps a(x) to -a(-x). The odd spatial forcing sin(x) is invariant under that full transformation. Negating a alone, without reflecting position or changing the forcing, would be a different test.

The 128/256 final-field relative differences are at most {max(checks['spatial_grid_relative_errors'].values()):.3g}. Halving dt from 0.005 to 0.0025 changes the tested strong-drive final field by {checks['spatial_dt_relative_error']:.3g}. These controls concern this scalar equation.

## 5. Fluid comparison: underway

The new fluid matrix uses the existing mirrored starting field, L=6, viscosity=0.01 and Heun integration through 0.4. It tests two forcing patterns:

- **Even gap:** opening and closing preserves reflection and does not itself supply a signed lean.
- **Odd lean:** the force has the antisymmetric template's direction; the reflected test reverses that force.

The amplitudes are 0.15 and 0.5 and periods are 0.2 and 0.1. The force has units U0/time and is included at both Heun stages with dt. This differs from the older uploaded per-step kick, whose strength changed with timestep.

Tests include mirrored pairs at 49 cubed, stronger-drive 65-grid and half-step controls, and a zero-bias even-gap control. Existing unforced controls are reused. [Fluid protocol](fluid-protocol.json).

The comparison records E, D, peak spin, its integral, spectra, divergence, energy input and the center-line half-peak gap. A gap diagnostic is not an imposed wall. Signed D and normalized unsigned E are not numerically identified with scalar x or E_model.

A scalar-to-fluid calibration will be fitted on the low-amplitude odd-drive run and tested on the other drive settings without refitting. Four-node, listening and spatial coefficients have no established numerical mapping to the fluid. Their symmetry predictions are compared with fluid controls; their coefficients are not presented as fluid measurements.

## Verification and code

- Ball tolerance-refinement discrepancy: {checks['ball_tolerance_refinement_max']:.3g}.
- Four-node mirror discrepancy after tighter integration: {res:.3g}.
- Spatial mirror discrepancy: {checks['spatial_mirror_max']:.3g}.
- All stored numbers are finite.

[Fixed reduced protocol](protocol.json) · [Reduced results](reduced-results.json) · [Verification](reduced-verification.json) · [Reduced runner](run_reduced.py) · [Fluid runner](run_fluid.py).

These are defined mechanism tests. A reduced model earns predictive use only when its held-out fluid comparison succeeds; structural symmetry alone does not establish that.
"""
(ROOT/'README.md').write_bytes(text.encode());print(json.dumps(checks))
