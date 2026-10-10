"""Make original scientific figures and a self-contained interactive report."""
from pathlib import Path
import base64
import html
import json
import hashlib
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from hug_model import Parameters, PRESETS, geometry, equilibria

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
COLORS=['#15758d','#b96938','#676f81','#77862f']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'figure.facecolor':'white','axes.facecolor':'white','axes.prop_cycle':plt.cycler(color=COLORS)})
def read(name):return json.loads((DATA/name).read_text(encoding='utf-8'))
def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=150,bbox_inches='tight')
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
def trajectory(key):
    z=np.load(DATA/f'preset-{key}.npz');return z['time'],z['state']

fig,axs=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for key in ['damped','overdamped']:
    t,y=trajectory(key);axs[0,0].plot(t,y[:,0],label='Light damping' if key=='damped' else 'Strong damping')
axs[0,0].set(title='Damping changes the return',xlabel='Model time',ylabel='Lean q');axs[0,0].legend()
for key in ['right','left','balanced']:
    t,y=trajectory(key);axs[0,1].plot(t,y[:,0],label={'right':'+0.04 start','left':'−0.04 start','balanced':'Exactly zero'}[key])
axs[0,1].set(title='Two possible resting shapes',xlabel='Model time',ylabel='Lean q');axs[0,1].legend()
t,y=trajectory('memory');p=PRESETS['memory'][1]
axs[1,0].plot(t,y[:,0],label='Lean');axs[1,0].plot(t,y[:,2],label='Imprint')
axs[1,0].axhline(.25/1.5,color=COLORS[2],ls=':',label='Imprint at stiffness change')
axs[1,0].set(title='Memory changes the restoring tendency',xlabel='Model time',ylabel='Model value');axs[1,0].legend()
t,y=trajectory('opening');p=PRESETS['opening'][1]
gap=[geometry(ti,yi,p,3)['gap'] for ti,yi in zip(t,y)]
axs[1,1].plot(t,gap,label='Tip separation');axs[1,1].plot(t,y[:,2],label='Imprint')
axs[1,1].set(title='Original opening rule remains latched',xlabel='Model time',ylabel='Model value');axs[1,1].legend()
save(fig,'hug-response')

fig,axs=plt.subplots(1,3,figsize=(11,3.8),constrained_layout=True)
for ax,key in zip(axs,['original','right','opening']):
    t,y=trajectory(key);p=PRESETS[key][1]
    ax.add_patch(plt.Rectangle((-3,-3),6,6,fill=False,color='#999999',lw=.8))
    ax.add_patch(plt.Rectangle((-2.55,-2.55),5.1,5.1,fill=False,color='#aaaaaa',lw=.6,ls=':'))
    for j,ti in enumerate([0,4,24]):
        i=int(round(ti/.01));g=geometry(ti,y[i],p)
        for k,half in enumerate(['left','right']):
            pts=np.array(g[half]);ax.plot(pts[:,0],pts[:,1],color=COLORS[j],label=f't = {ti}' if k==0 else None)
    ax.set(xlim=(-3.2,3.2),ylim=(-3.2,3.2),aspect='equal',xlabel='x · model length',ylabel='y · model length',title={'original':'Pressure and imprint','right':'Added lean','opening':'Pressure-triggered opening'}[key]);ax.legend(fontsize=8)
fig.suptitle('Existing arm geometry, with an optional area-preserving lean\nSquare: periodic-domain guide; dotted square: gate collar edge',fontsize=12)
save(fig,'hug-shapes')

fig,axs=plt.subplots(1,3,figsize=(11,3.5),constrained_layout=True)
for ax,b in zip(axs,[0.,.2,2.]):
    x=np.linspace(-np.pi,np.pi,17);v=np.linspace(-2.5,2.5,13);X,V=np.meshgrid(x,v)
    U=V;W=-np.sin(X)-b*V;norm=np.sqrt(U*U+W*W);norm[norm==0]=1
    ax.quiver(X,V,U/norm,W/norm,color='#b9bdc2',scale=30)
    for v0 in [0.,2.2]:
        z=np.load(DATA/f'pendulum-{b}-{v0}.npz');Y=z['state'];a=(Y[:,0]+np.pi)%(2*np.pi)-np.pi
        a[1:][abs(np.diff(a))>np.pi]=np.nan
        ax.plot(a,Y[:,1],label=f'Initial v = {v0}')
    ax.set(xlim=(-np.pi,np.pi),ylim=(-2.7,2.7),title=f'Pendulum damping = {b:g}',xlabel='Wrapped angle θ · radians',ylabel='Angular speed · model units')
    ax.set_xticks([-np.pi,0,np.pi],['−π','0','π']);ax.legend(fontsize=8)
save(fig,'pendulum')

fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for row in read('numerical-orders.json'):
    axs[0].loglog(row['steps'],row['errors'],'o-',label=f"{row['method'].upper()} · order {row['global_order']:.2f}")
for method in ['euler','heun','rk4']:
    rows=[r for r in read('hug-integrators.json') if r['method']==method]
    axs[1].loglog([r['dt'] for r in rows],[r['max_error'] for r in rows],'o-',label=method.upper())
axs[0].set(title='Known exponential solution',xlabel='Time step',ylabel='Endpoint absolute error')
axs[1].set(title='Hug compared with independent solver',xlabel='Time step',ylabel='Maximum state error')
for ax,ticks in zip(axs,[[.025,.05,.1,.2],[.01,.02,.04,.08]]):
    ax.legend();ax.grid(alpha=.18)
    ax.set_xticks(ticks,[f'{x:g}' for x in ticks]);ax.xaxis.set_minor_formatter(NullFormatter())
save(fig,'convergence')

fig,axs=plt.subplots(1,3,figsize=(11,3.8),constrained_layout=True)
for ax,r in zip(axs,[-1,0,1]):
    q=np.linspace(-1.6,1.6,150);v=np.linspace(-1.5,1.5,120);Q,V=np.meshgrid(q,v)
    ax.streamplot(q,v,V,r*Q-Q**3-V,color='#9199a0',density=.65,linewidth=.7,arrowsize=.8)
    for eq,eig in equilibria(r,1):
        index=-1 if r>0 and eq==0 else 1
        ax.scatter(eq,0,c=COLORS[0] if index==1 else COLORS[1],s=35,zorder=3)
        ax.annotate(f'index {index:+d}',(eq,0),xytext=(0,13),textcoords='offset points',ha='center',fontsize=9)
    ax.set(title=f'After imprint fades · r = {r}',xlabel='Lean q',ylabel='Lean speed v',xlim=(-1.6,1.6),ylim=(-1.5,1.5))
save(fig,'stability')

fig,axs=plt.subplots(1,3,figsize=(11,3.8),constrained_layout=True)
for ax,law in zip(axs,['logistic','gompertz','allee']):
    for row in read('growth.json'):
        if row['law']==law:
            a=np.array(row['samples']);ax.plot(a[:,0],a[:,1],label=f"Start {row['x0']}")
    ax.set(title=law.capitalize(),xlabel='Model time',ylabel='Scalar state',ylim=(0,1.6));ax.legend(fontsize=8)
save(fig,'growth')

fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
t=np.linspace(0,4,400)
for j,delay in enumerate([0,1,2]):
    for sign in [-1,1]:axs[0].plot(t,sign*(2*np.maximum(t-delay,0)/3)**1.5,color=COLORS[j],label=f'Departure time {delay}' if sign==1 else None)
axs[0].plot(t,t*0,color='#333333',ls='--',label='Remain at zero')
axs[0].set(title='Same start, many valid solutions',xlabel='Time',ylabel='x');axs[0].legend(fontsize=8)
r=np.linspace(-1,1.5,300);positive=r>0
axs[1].plot(r[~positive],r[~positive]*0,color=COLORS[0],label='Stable branch')
axs[1].plot(r[positive],r[positive]*0,ls='--',color=COLORS[1],label='Unstable branch')
axs[1].plot(r[positive],np.sqrt(r[positive]),color=COLORS[0]);axs[1].plot(r[positive],-np.sqrt(r[positive]),color=COLORS[0])
axs[1].set(title='Buckling-style change of equilibrium',xlabel='Restoring parameter r',ylabel='Equilibrium lean q');axs[1].legend(fontsize=8)
save(fig,'uniqueness-buckling')

v=read('verification.json');presets=read('presets.json');summary=read('preset-checks.json')
fragment=(ROOT/'visual-template.html').read_text(encoding='utf-8').replace('__HUG_DATA__',json.dumps(presets,separators=(',',':'),allow_nan=False))
(ROOT/'hug-dynamics-fragment.html').write_text(fragment,encoding='utf-8')

def figure(name,caption):
    image=base64.b64encode((FIG/(name+'.png')).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{image}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'

mapping=[
 ('1–2','Pendulum cylinder and damping','Independent periodic-angle benchmark; six completed trajectories. It is not substituted for the hug’s nonperiodic lean.'),
 ('3','Vector-field index','Measured winding around equilibria of the lean equation; sink +1, saddle −1. This measures vector rotation, not passage through a membrane.'),
 ('4–6','Linear stability and zero derivative','Eigenvalues classify the ordinary cases. A zero eigenvalue is explicitly inconclusive; the cubic and quartic restoring terms must be examined.'),
 ('7–8','Existence and uniqueness (duplicate page)','Six analytic delayed-departure solutions plus the stationary solution verify the cube-root counterexample. That nonunique law is excluded from the hug.'),
 ('9','Overdamped motion','Small-mass spring controls approach the first-order relaxation solution. Added velocity is what permits overshoot and oscillation in the hug extension.'),
 ('10–11, 13','Euler, improved Euler and RK4','Executable methods, known-solution order checks, independent-solver comparisons and time-step refinement. Improved Euler is Heun, already used in the imported fluid solver.'),
 ('12','Autocatalytic, Gompertz and Allee laws','Twelve scalar trajectories; logistic/autocatalytic and Gompertz compared with analytic solutions. Candidate feedback modules only; no chemistry or tumor data fitted.'),
 ('14','Buckling','Cubic lean equation has one stable central rest for r < 0 and two stable side rests for r > 0. Opposite seeds and exact-zero controls tested.')]
table='<table><tr><th>Your images</th><th>Concept</th><th>What was built and checked</th></tr>'+''.join('<tr>'+''.join(f'<td>{html.escape(x)}</td>' for x in row)+'</tr>' for row in mapping)+'</table>'
rows=''.join(f'<tr><td>{html.escape(r["label"])}</td><td>{r["final_q"]:.6g}</td><td>{r["peak_abs_q"]:.6g}</td><td>{r["reference_error"]:.3g}</td></tr>' for r in summary)
body=f'''<h1>The hug with pressure memory and dynamics</h1>
<p>This extension starts with your repository’s mirrored arms, pressure threshold, fading imprint and reciprocal axis stretch. It adds a signed lean and its velocity, so the hug can settle, oscillate or approach a buckled shape. The diagrams below come from completed calculations.</p>
<p><strong>{v['passed']} of {v['checks']} checks passed.</strong> The stress sweep contains 81 configurations, each at two time steps. These results validate the implementation against mathematical controls. They do not establish that this is a measured material law or a validated fluid-boundary model.</p>
{fragment}
<h2>What changed</h2>
<p>The original pressure-memory model remains available as the first experiment. The new lean is a proposed extension with explicitly chosen coefficients. A square guide makes the existing starting gate visible. Turning its display off does not change the calculated motion. The gate prepares a fluid start in the repository; it is not an evolving wall or source of force in this reduced experiment.</p>
{figure('hug-shapes','The same arm construction is reused. The added shear has determinant one, preserving the area of a closed shape. Once the arms open, no enclosed-area claim is made.')}
<h2>Equations and choices</h2>
<pre>p(t) = p_max sin²(πt/4), 0 &lt; t &lt; 4; otherwise zero
m′ = (p − m)/τ, τ = 0.8 while loading, 6 while unloading
q′ = v
v′ = r q − q³ − damping v + coupling m q
s = 1 + m [0.04 + 0.02 sin(2πt/3)]
(X,Y) = (s x + 0.35 q y/s, y/s)</pre>
<p>The base arm coordinates (x,y), pressure-normalized opening threshold 1, and 0.8-duration opening transition come from the original code. Opening stays latched after the threshold is reached. The smooth pressure pulse, lean equation, memory-to-stiffness coupling and shear coefficient are new choices. No pressure, length, temperature or time calibration is supplied; all extension results use model units.</p>
<p>The pressure and memory enter the lean equation symmetrically: a scalar pressure does not choose left or right. With q=v=0 exactly, the solution stays balanced even when the central state becomes unstable. A small signed disturbance is needed to choose a side. The first-order memory law alone is not an autonomous oscillator; the added velocity supplies a second state.</p>
{figure('hug-response','Response curves for the original model and the added motion. Memory persists after the pressure pulse. Pressure-triggered opening remains released.')}
<table><tr><th>Experiment</th><th>Final lean at 24</th><th>Peak absolute lean</th><th>Maximum error against DOP853</th></tr>{rows}</table>
<h2>How your pages were incorporated</h2>{table}
{figure('stability','Phase portraits after imprint has faded. These are autonomous slices of the lean equation, not portraits of the full time-dependent memory system. At r=0 linearization has a zero eigenvalue; the nonlinear restoring term matters.')}
{figure('pendulum','Independent pendulum controls. The angle wraps from π to −π; opposite plot edges are identified, giving the cylinder interpretation. Conservative motion preserves E = v²/2 − cos θ; damping removes energy.')}
{figure('uniqueness-buckling','Left: analytic solutions of x′ = real cube root of x all start at zero, so uniqueness fails. Right: equilibria of the proposed cubic lean law; this is the mathematical buckling pattern, not a calibrated beam model.')}
{figure('growth','Separate scalar response-law demonstrations. With the chosen constants, the autocatalytic equation equals the logistic equation. An Allee threshold separates return toward zero from approach to the upper rest state.')}
<h2>Numerical error and energy accounting</h2>
<p>The largest state difference between time steps 0.02 and 0.01 across the 81 configurations was <strong>{v['maxima']['sweep_step_difference']:.8g}</strong> in absolute model units. The largest error across eight displayed presets against the independent DOP853 calculation was <strong>{v['maxima']['preset_reference']:.8g}</strong>. DOP853 tolerances were tightened independently. These are numerical discrepancies, not experimental errors or statistical confidence intervals.</p>
<pre>E = v²/2 − r q²/2 + q⁴/4
E(t) − E(0) = ∫ coupling m q v dt − ∫ damping v² dt</pre>
<p>The largest sweep energy-budget residual was {v['maxima']['sweep_energy_budget']:.8g}. This budget covers the lean oscillator. It explicitly counts work transferred through the prescribed memory coupling; it is not the energy of a fluid or a proof of a passive membrane.</p>
{figure('convergence','The known exponential solution shows first-, second- and fourth-order global convergence. RK4 one-step local error is fifth order. In the hug, the loading/unloading switch can reduce observed formal order; measured agreement is reported directly.')}
<h2>What can be added to your model</h2>
<p>Damping, an inertial lean, a buckling parameter, local stability checks, vector-field index checks and numerical convergence controls now run together in the reduced hug extension. The pendulum and scalar growth laws are executable comparison modules. Their variables have not been assigned a physical meaning inside the hug, so they are not silently substituted into it.</p>
<p>The supplied cubic and memory laws are locally Lipschitz for bounded states and prescribed continuous pressure, giving local uniqueness between opening events. The opening event has a deterministic latched rule. This contrasts with the cube-root example. No all-time Navier–Stokes regularity statement follows.</p>
<h2>Connection still to be tested</h2>
<p>The existing Navier–Stokes gate acts at initialization and the domain is periodic. This package does not add a moving wall or a permeable interface to that solver. Such an extension needs a tracked boundary, a fluid-to-boundary force rule, a defined transport law, equal-and-opposite force/energy accounting and matched tests with feedback switched off. No claim that the square boundary itself creates memory, that matter crosses these arms, or that the model reproduces experimental fluid behavior is supported by these reduced tests.</p>
<p>The parameter sweeps describe sensitivity within chosen equations. They are not samples from measured uncertainty. Model accuracy, spatial convergence of a coupled fluid-boundary calculation, and temperature dependence remain untested.</p>
<h2>Reproduce</h2><pre>python -m pip install -r requirements.txt
python run_study.py
python build_report.py</pre>
<p>The data folder contains every sweep trajectory at both resolutions, displayed trajectories and independent references, scalar and pendulum controls, protocol, source hashes, numerical checks and environment versions. Original source copies are preserved in source/. No external service is needed to open this report.</p>
<h2>Sources</h2>
<ul><li>Your 14 supplied screenshots, with images 7 and 8 duplicates: pendulum, index, stability, uniqueness, damping, integration, scalar growth and buckling. Figure and exercise numbers are mapped above; the edition was not supplied. The screenshots are not redistributed.</li>
<li><a href="https://github.com/NousVolition/My-Sources-Project-ALL/blob/f5a2106f69bd1c4a0f8eb7a64f5291952c4c5e25/math/notes/soft-envelope.md">Your pressure, imprint and soft-envelope model</a>.</li>
<li><a href="https://github.com/NousVolition/My-Sources-Project-ALL/blob/f5a2106f69bd1c4a0f8eb7a64f5291952c4c5e25/math/imports/hug-ns/corrected_v1/solver.py">Your corrected initial gate and fluid solver</a>.</li>
<li><a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html">SciPy solve_ivp documentation</a>: independent high-accuracy ODE reference, successful completion checked.</li>
<li><a href="https://math.mit.edu/~dyatlov/18.03/">MIT differential equations course</a>: pendulum companion system, critical points and phase portraits.</li></ul>'''
css='''body{font:16px/1.6 system-ui,sans-serif;margin:0 auto;padding:28px;max-width:1100px;color:#21343b;background:#fff}h1,h2{line-height:1.2;font-weight:600}h2{margin-top:40px}p{max-width:92ch}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;border-bottom:1px solid #d9e0e3;padding:9px;vertical-align:top}pre{white-space:pre-wrap;background:#f3f5f6;padding:16px;overflow-wrap:anywhere}img{max-width:100%;height:auto}figure{margin:28px 0}figcaption{font-size:14px;color:#4b616a}a{color:#126782}button,input,select{font:inherit}button{padding:9px 14px;cursor:pointer}select{max-width:100%;padding:8px}input[type=range]{width:100%}.viz-controls,.viz-row{display:flex;flex-wrap:wrap;gap:14px;align-items:center}.form-label{display:block;margin:12px 0}.form-check{display:flex;gap:8px;align-items:center}.form-select{display:block}.text-small{font-size:14px}.text-muted{color:#54666f}.tabular-nums{font-variant-numeric:tabular-nums}:root{--foreground:#21343b;--muted-foreground:#60727a;--border:#cbd4d8;--viz-series-1:#15758d}#hug-dynamics-lab{max-width:736px;margin:32px auto} @media(max-width:600px){body{padding:16px}table{font-size:12px}th,td{padding:6px}}'''
document='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>The hug dynamics lab</title><style>'+css+'</style><main>'+body+'</main></html>'
(ROOT/'report.html').write_text(document,encoding='utf-8')

readme=f'''# The hug with pressure memory and dynamics

A tested extension of the existing hug geometry and pressure-memory prototype, using the supplied pages on damping, stability, pendulums, numerical methods and buckling.

**{v['passed']}/{v['checks']} checks passed.** 81 parameter configurations were run at two time steps. The largest time-step discrepancy was {v['maxima']['sweep_step_difference']:.8g}; the eight displayed presets differed from an independent solver by at most {v['maxima']['preset_reference']:.8g} in model units.

![Hug response](figures/hug-response.png)

Download and open [report.html](report.html) to play the original hug, damped motion, buckling, memory feedback and pressure opening. GitHub's file viewer displays its source; the downloaded report works offline.

The unchanged original pressure model is a control. The extension adds a signed lean and velocity, cubic restoring force, damping and an optional memory-dependent stiffness. The original pressure limit still opens the arms irreversibly. The square guide shows the existing initial gate; it is not an active wall.

## What the supplied pages add

- Damping, oscillation, buckling and stability calculations now run in the reduced hug extension.
- Euler, Heun and RK4 are compared with exact and independent reference solutions.
- Pendulum phase and energy, vector-field winding, degenerate equilibria, nonuniqueness, logistic/autocatalytic, Gompertz and Allee examples have executable checks or analytic controls. The report maps all 14 images, including the duplicate.
- Every sweep keeps both time resolutions in `data/`, alongside reference trajectories, protocol, provenance and verification.

## Limits

This is a proposed reduced dynamics model, with chosen coefficients and no material calibration. It does not yet couple an evolving boundary to Navier–Stokes, simulate permeability, or predict measured water behavior. The mathematical consistency checks do not establish those physical claims. The original fluid calculations remain unchanged.

## Reproduce

```sh
python -m pip install -r requirements.txt
python run_study.py
python build_report.py
```

[Model](hug_model.py) · [Completed verification](data/verification.json) · [Fixed protocol](data/protocol.json) · [Original source hashes](data/sources.json) · [Full report](report.html)

The source baseline is My-Sources-Project-ALL commit `f5a2106f69bd1c4a0f8eb7a64f5291952c4c5e25`. Original pressure/geometry tests are executed without changing their assertions. The report cites the supplied pages and primary numerical documentation.
'''
(ROOT/'README.md').write_text(readme,encoding='utf-8')
print(json.dumps(dict(report_bytes=(ROOT/'report.html').stat().st_size,fragment_bytes=(ROOT/'hug-dynamics-fragment.html').stat().st_size,figures=7)))
