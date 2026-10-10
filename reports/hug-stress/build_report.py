"""Build static scientific figures and a self-contained offline stress report."""
from pathlib import Path
import json,base64,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from stress_solver import model

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';FIG=ROOT/'figures'
def read(name):return json.loads((DATA/name).read_text())
def arrays(name):return np.load(DATA/(name+'.npz'))
def main():
    FIG.mkdir(exist_ok=True)
    rows=read('screen.json');summary=read('summary.json');tight=read('long-refinement.json')
    byname={r['name']:r for r in rows}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.titleweight':'bold','savefig.facecolor':'white'})
    blue='#166985';orange='#bd582c';green='#248369';red='#b62e42';gray='#677780'
    def save(fig,name):
        fig.savefig(FIG/(name+'.png'),dpi=170,bbox_inches='tight');fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    passed=[summary['coarse_pass'],summary['fine_pass']];stopped=[summary['coarse_stopped'],summary['fine_stopped']]
    missed=[224-p-s for p,s in zip(passed,stopped)]
    axs[0].bar(['Step 0.02','Step 0.01'],passed,color=green,label='Met accuracy target')
    axs[0].bar(['Step 0.02','Step 0.01'],missed,bottom=passed,color=orange,label='Finished, missed target')
    axs[0].bar(['Step 0.02','Step 0.01'],stopped,bottom=np.array(passed)+missed,color=red,label='Numerical guard stop')
    for i in range(2):
        for y,n in [(passed[i]/2,passed[i]),(passed[i]+missed[i]/2,missed[i]),(224-stopped[i]/2,stopped[i])]:axs[0].text(i,y,str(n),ha='center',va='center',color='white',weight='bold')
    axs[0].set(ylabel='Configurations',title='Smaller steps help, but do not ensure accuracy',ylim=(0,255));axs[0].legend(loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False)
    valid=[r for r in rows if r['fixed'][1]['failure'] is None]
    for label,selection,color in [('Other completed runs',[r for r in valid if not r['false_convergence']],blue),('Misleading step agreement',[r for r in valid if r['false_convergence']],red)]:
        axs[1].scatter([r['parameters']['pulse_duration'] for r in selection],[max(1e-15,r['fixed'][1]['scaled_reference_error']) for r in selection],s=24,alpha=.75,color=color,label=label)
    axs[1].axhline(1e-3,color='black',ls='--',lw=1,label='Accuracy target')
    axs[1].set(xscale='log',yscale='log',xlabel='Pressure pulse duration (model time)',ylabel='Scaled state error',title='52 cases agree with the wrong answer')
    axs[1].legend(loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False)
    save(fig,'stress-overview')

    a=arrays('missed-pulse');fix=arrays('pulse-resolved-0.01');p=model.Parameters(**byname['missed-pulse']['parameters'])
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    ts=np.linspace(0,.001,401)
    axs[0].plot(ts,[model.pressure(t,p) for t in ts],color=orange,label='Applied pressure')
    axs[0].scatter([0,.005],[0,0],color=red,label='First RK4 sample locations, h = 0.01')
    axs[0].set(xlim=(-.00015,.0052),xlabel='Model time',ylabel='Pressure (model units)',title='The entire pulse fits between samples');axs[0].legend(fontsize=8)
    axs[1].plot(a['reference_time'],a['reference_state'][:,2],color=blue,lw=2.5,label='Pulse-resolving reference')
    axs[1].plot(fix['time'],fix['state'][:,2],color=green,ls='--',label='RK4 with pulse subdivision')
    axs[1].plot(a['fine_time'],a['fine_state'][:,2],color=red,label='Both original step sizes: zero')
    axs[1].set(xlim=(0,24),xlabel='Model time',ylabel='Fading imprint m',title='True peak imprint: 0.0624614');axs[1].legend(fontsize=8)
    save(fig,'missed-pulse')

    a=arrays('stiff-damping');fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for tag,color in [('coarse',orange),('fine',red)]:axs[0].plot(a[tag+'_time'],a[tag+'_state'][:,0],'.-',color=color,label='RK4 '+str(byname['stiff-damping']['fixed'][0 if tag=='coarse' else 1]['dt']))
    axs[0].axhline(2,color=blue,ls='--',label='Reference stays close to 2 at early times')
    axs[0].set(yscale='symlog',xlabel='Model time',ylabel='Lean q (symmetric log scale)',title='Numerical runaway at strong damping');axs[0].legend(fontsize=8)
    axs[1].plot(a['reference_time'],a['reference_state'][:,0],color=blue,label='DOP853')
    b=arrays('stiff-damping-radau');axs[1].plot(b['time'],b['state'][:,0],color=green,ls='--',label='Independent Radau')
    axs[1].set(xlabel='Model time',ylabel='Lean q',title='The equations give bounded, slow relaxation');axs[1].legend()
    save(fig,'stiff-damping')

    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for name,label,color in [('boundary-no-feedback','Pressure stretch, no lean feedback',blue),('boundary-escape','Pressure stretch with lean feedback',orange)]:
        a=arrays(name);axs[0].plot(a['reference_time'],a['extent'],color=color,label=label)
    axs[0].axhline(3,color=red,ls='--',label='Square guide at ±3');axs[0].set(xlabel='Model time',ylabel='Maximum coordinate magnitude',title='29 configurations extend outside the guide');axs[0].legend(fontsize=8)
    name='boundary-no-feedback';a=arrays(name);i=int(np.argmax(a['extent']));p=model.Parameters(**byname[name]['parameters']);g=model.geometry(a['reference_time'][i],a['reference_state'][i],p,points=257)
    for key in ('left','right'):
        arm=np.array(g[key]);axs[1].plot(arm[:,0],arm[:,1],color=blue,lw=3)
    axs[1].add_patch(Rectangle((-3,-3),6,6,fill=False,edgecolor=red,linestyle='--'))
    axs[1].set(xlim=(-7.6,7.6),ylim=(-3.5,3.5),xlabel='x (model units)',ylabel='y (model units)',title=f'No-feedback shape at t = {a["reference_time"][i]:.2f}');axs[1].set_aspect('equal')
    save(fig,'square-guide')

    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    threshold=sorted([r for r in rows if r['group']=='threshold'],key=lambda r:r['parameters']['pressure'])
    for row,color in zip(threshold,[blue,orange,green]):
        p=model.Parameters(**row['parameters']);t=np.linspace(1.7,3,501);op=model.opening_time(p)
        b=np.ones_like(t) if op is None else np.array([1-model.smoothstep((x-op)/.8) for x in t])
        axs[0].plot(t,b,color=color,ls='--' if p.pressure>1 else '-',label=f'Peak pressure {p.pressure:.8f}')
    axs[0].set(xlabel='Model time',ylabel='Arm bend (1 = closed, 0 = open)',title='The prescribed opening decision has a threshold');axs[0].legend(fontsize=8)
    for name,label,color in [('bias-0','Exactly balanced',gray),('bias-1','Starting lean −10⁻¹²',blue),('bias-2','Starting lean +10⁻¹²',orange)]:
        a=arrays(name);axs[1].plot(a['reference_time'],a['reference_state'][:,0],color=color,label=label)
    axs[1].set(xlabel='Model time',ylabel='Lean q',title='Tiny signed starts choose opposite resting shapes');axs[1].legend(fontsize=8)
    save(fig,'threshold-and-bias')

    fig,axs=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    for ax,name in zip(axs,['long-conservative','long-weak-damping']):
        a=arrays(name);b=arrays(name+'-tight-Radau');p=model.Parameters(**byname[name]['parameters']);scale=np.maximum(1,np.max(abs(b['state'][:,:3]),axis=0))
        for tag,col in [('coarse',orange),('fine',blue)]:
            # Plot the same common observation grid, avoiding interpolation artifacts at off-grid pulse sample times.
            t=a[tag+'_time'];idx=np.rint(b['time']/(t[1]-t[0])).astype(int);mask=np.abs(t[idx]-b['time'])<1e-8
            err=np.max(abs(a[tag+'_state'][idx[mask],:3]-b['state'][mask,:3])/scale,axis=1)
            ax.plot(b['time'][mask],np.maximum(1e-16,err),color=col,lw=.8,label='RK4 '+str(.02 if tag=='coarse' else .01))
        ax.axhline(1e-3,color='black',ls='--',label='Accuracy target');ax.set(yscale='log',ylim=(1e-10,.05),xlabel='Model time',ylabel='Scaled state error',title=name.replace('long-','').replace('-',' ').capitalize());ax.legend(fontsize=8)
    save(fig,'long-runs')

    captions={
      'stress-overview':'The pass target is a maximum scaled q, v, m error of 0.001. Counts describe this deliberately harsh finite test set; they are not probabilities of failure in use. Stopped attempts are counted separately from inaccurate completed trajectories.',
      'missed-pulse':'A pulse lasting 0.001 model time with peak pressure 100 is missed by RK4 steps 0.02 and 0.01. Their comparison is exactly zero, although both miss an imprint peak of 0.0624614. Subdividing the pulse into at least 64 steps reduces the checked error to 1.58×10⁻⁷ without changing the equations.',
      'stiff-damping':'With damping 500 and no pressure, both original step sizes trigger the numerical guard at t=0.08. Independent adaptive calculations instead take q from 2 to 1.63094 by t=24. The large damping makes explicit time integration stiff; this is not physical failure.',
      'square-guide':'The greatest sampled coordinate extent is 8.39098, against the guide at 3. Even with lean feedback disabled, pressure stretch reaches 6.97067. The existing square has no confining force, so leaving it is expected from these equations. It does not establish rupture or escaped fluid.',
      'threshold-and-bias':'Changing peak pressure from 0.99999999 to 1 changes the latched opening outcome. This is an imposed rule, not an experimentally identified material transition. Starts of ±10⁻¹² eventually reach opposite equilibria; exact zero remains balanced. Sensitivity near an unstable equilibrium does not by itself imply chaos.',
      'long-runs':'Long runs reach t=1000. Curves compare the original RK4 states with tightened Radau calculations. Initial adaptive cross-checks failed the stricter 2×10⁻⁶ reference target; the original failures remain in the data. The addendum tightens both methods by 100× and adds an analytic elliptic-function solution for the conservative case.'}
    img=lambda n:'data:image/png;base64,'+base64.b64encode((FIG/(n+'.png')).read_bytes()).decode()
    sections=''.join(f'<section><h2>{title}</h2><img src="{img(n)}" alt="{html.escape(title)}"><p>{captions[n]}</p></section>' for n,title in [
        ('stress-overview','How far the existing solver can be trusted'),('missed-pulse','Two steps can agree and both be wrong'),('stiff-damping','More damping can require a different solver'),('square-guide','The invisible boundary is not yet a responding wall'),('threshold-and-bias','Thresholds and small disturbances'),('long-runs','Long-term checks and reference uncertainty')])
    reference_rows=''.join(f'<tr><td>{x["name"]}</td><td>{x["scaled_discrepancy"]:.3e}</td><td>{x["original_reference_change"]:.3e}</td><td>{x["rk4_errors"]["fine"]:.3e}</td></tr>' for x in tight)
    htmltext=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>The hug under stress</title>
<style>body{{overflow-wrap:anywhere;font:17px/1.65 system-ui,sans-serif;color:#24353e;background:#f5f7f8;margin:0}}main{{max-width:1050px;margin:auto;padding:40px 24px}}h1{{font-size:42px;line-height:1.1}}h2{{line-height:1.3}}.lead{{font-size:21px}}section{{background:white;border:1px solid #dce3e6;border-radius:14px;padding:25px;margin:28px 0}}img{{width:100%;height:auto}}.numbers{{display:flex;flex-wrap:wrap;gap:18px}}.number{{flex:1;min-width:130px;background:#e6eef1;padding:18px;border-radius:10px}}.number strong{{font-size:32px;display:block}}table{{border-collapse:collapse;font-size:14px;width:100%}}td,th{{text-align:left;padding:9px;border-bottom:1px solid #ccd7dc}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#eef2f4;padding:16px}}a{{color:#106782}}.scroll{{overflow:auto}}</style>
<main><p>HUG DYNAMICS · EXTREME NUMERICAL TESTS · 2026-10-10 UTC</p><h1>The hug under stress</h1>
<p class="lead">The reduced hug equations survive the checked conditions, but the earlier fixed-step calculation does not. The square still needs an actual response law if it is meant to contain the shape.</p>
<div class="numbers"><div class="number"><strong>224</strong>configurations</div><div class="number"><strong>448</strong>original RK4 attempts</div><div class="number"><strong>52</strong>misleading agreements</div><div class="number"><strong>29</strong>outside the guide</div></div>
<p>These are completed numerical experiments on the proposed reduced model. All values use chosen model units. No material, permeability, water measurement or fluid-boundary behavior has been validated here.</p><p><strong>Added after the Lorenz request:</strong> <a href="lorenz-report.html">Play the new Lorenz-driven hug</a>, with a separate 84-case screen and sensitivity tests. The 224 cases below test the original prescribed pulse.</p>
{sections}
<section><h2>Independent checks and uncertainty</h2><p>All 224 configurations have pulse-resolving DOP853 reference calculations. A selected 43 were also run with the independent implicit Radau method. At the initial tolerances, 41 of 43 met the independent comparison target; the two long runs did not. Eleven targeted controls passed. None of these counts removes the fixed-step failures.</p>
<div class="scroll"><table><tr><th>Tightened long run</th><th>Radau vs DOP853</th><th>Original reference change</th><th>Fine RK4 vs tight Radau</th></tr>{reference_rows}</table></div>
<p>The conservative run additionally agrees with its analytic Jacobi-cn solution: tight DOP853 error {tight[0]['exact_errors']['DOP853']:.3e}; tight Radau error {tight[0]['exact_errors']['Radau']:.3e}. Both refined cross-checks {'met' if all(x['passed'] for x in tight) else 'did not meet'} the unchanged 2×10⁻⁶ target. The coarse/fine pass classifications for these long runs remain unchanged.</p>
<p>Error means maximum state difference over saved times, with each q, v, m component divided by max(1, its reference peak magnitude). It is not a relative percentage near zero. The energy residual is separately scaled by max(1, peak energy, work, loss). Peaks and first exits are sampled in time; no exact event-time error bound or confidence interval is claimed. Independent solver agreement estimates numerical sensitivity, not measurement uncertainty.</p></section>
<section><h2>What was stressed</h2><p>The fixed matrix covers 144 combinations: stiffness r = −4, 0, 4; damping = 0, 0.05, 4, 500; peak pressure = 0.5, 1.0001, 100; pulse duration = 0.004 or 4; memory coupling = 0 or 8. Another 64 configurations use a seeded Latin-hypercube design (seed 20261010), followed by 16 adversarial, threshold, bias and long-run controls. Every original trial uses steps 0.02 and 0.01.</p>
<p>The sampled design extends stiffness from −10 to 10, damping approximately 0.001 to 501, pressure 0.001 to 100, pulse durations 0.001 to 50.1, and coupling 0.001 to 31.6. Starting lean spans very small signed disturbances to order-one changes; initial velocity also varies. Long controls reach t=1000. This is a reproducible exploration, not exhaustive coverage.</p></section>
<section><h2>The equations tested</h2><pre>q′ = v
v′ = r q − q³ − γ v + κ m q
m′ = (p − m) / τ, with τ = 0.8 when p &gt; m, otherwise 6
p(t) = P sin²(πt/T) for 0 &lt; t &lt; T; zero otherwise
E = v²/2 − r q²/2 + q⁴/4
E(t) − E(0) = ∫ κ m q v dt − ∫ γ v² dt</pre>
<p>The original mirrored arms, breathing stretch and pressure-opening threshold are retained from commit 991704f2c4d708d2ab3fd29532ce176426ef78c9. Opening is irreversible in this prototype, and pressure is still prescribed after opening. The lean-energy identity is not a budget for an enclosing material or fluid.</p>
<p>There is also an analytic finite-time bound for these reduced equations. For bounded nonnegative pressure, m(0)=0 gives 0≤m≤P. With γ≥0 and κ≥0, define V=v²/2+q⁴/4+q²/2. Then V′=(r+1+κm)qv−γv²≤(|r+1|+κP)V, so V(t)≤V(0)exp((|r+1|+κP)t). The locally Lipschitz vector field therefore continues for every finite time under these assumptions. Numerical runaway here is not a finite-time singularity. This result is specific to this small system and says nothing about global regularity of three-dimensional Navier–Stokes.</p></section>
<section><h2>How this differs from the Lorenz image</h2><p>The supplied Lorenz screenshot discusses sensitivity to initial conditions in a particular convection model. Numerical runaway from a time step that is too large is a different phenomenon. The hug's opposite resting states from tiny signed starts show sensitivity near an unstable equilibrium, but these tests do not establish a chaotic attractor.</p><p>The Lorenz variables and control parameter are not the hug's variables and stiffness parameter, even when both are written with the letter r. A Lorenz bifurcation threshold cannot be transferred to these hug equations. A chaos claim would need additional diagnostics, such as a positive long-time Lyapunov exponent that persists under solver refinement, together with bounded nonperiodic dynamics. No such claim is made here. See <a href="https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml">Lorenz's original 1963 study</a>.</p></section>
<section><h2>What this means for the hug idea</h2><p>The existing model represents yielding lean, damping, imposed opening and a fading imprint. Its numerical implementation needs pulse resolution and a stiffness-aware solver. Those changes preserve its equations.</p>
<p>A breathable enclosing shield additionally needs a specified boundary response, a pressure difference that changes with deformation and venting, and a transport law for what passes through. Contact, irreversible yielding and fracture would each need their own definitions if desired. Clipping the drawing to the square would hide the missing mechanics rather than test them.</p>
<p>No active boundary was added in this study. There is no temperature variable, calibrated length/time scale, grid of fluid velocities, density, viscosity, or nuclear-quantum effect in these equations. These results do not validate the earlier quantum-water or Navier–Stokes claims, and the separate fluid experiments have not been altered.</p></section>
<section><h2>Reproduce and inspect</h2><pre>python -m pip install -r requirements.txt
python run_stress.py
python refine_long_runs.py
python build_report.py
python verify_package.py</pre><p>The initial study intentionally exits with status 1 because the two original long reference comparisons fail their target. The refinement is a recorded follow-up, not a replacement of those results. For a full recomputation, work in a copy and remove its data directory first; otherwise the screen reuses saved case results. Independent comparisons and controls rerun.</p>
<p>Raw trajectories, failure records, parameter design, environment versions, source hashes and the editorial correction to a factorial-count comment are saved in data/. NPZ trajectories have time and state arrays; state columns are q, v, imprint, accumulated work and dissipation. No failed trials were removed.</p>
<p><a href="https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/hug-stress">GitHub package</a> · <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html">SciPy solver documentation</a> · <a href="https://dlmf.nist.gov/22.13">NIST elliptic-function differential equations</a></p></section></main></html>'''
    (ROOT/'report.html').write_text(htmltext,encoding='utf-8')
    print('Built 6 scientific figures in PNG/SVG and report.html')

if __name__=='__main__':main()
