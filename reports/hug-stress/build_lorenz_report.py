"""Build the added mode's recorded playback and scientific evidence."""
from pathlib import Path
import json,base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'lorenz-data';FIG=ROOT/'figures'
def read(n):return json.loads((DATA/n).read_text())
def main():
    presets=read('presets-matched.json');summary=read('summary.json');extra=read('sweep-independent.json');event=read('event-check.json')
    for record in presets:
        a=np.load(DATA/record.get('trajectory_file',record['key']+'.npz'))
        record['trace']=np.round(a['state'][:,[0,2]],6).tolist()
    fragment=(ROOT/'lorenz-template.html').read_text(encoding='utf-8').replace('__LORENZ_DATA__',json.dumps(presets,separators=(',',':')))
    (ROOT/'lorenz-hug-fragment.html').write_text(fragment,encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.titleweight':'bold'})
    def save(fig,name):
        for ext in ('png','svg'):fig.savefig(FIG/(name+'.'+ext),dpi=160,bbox_inches='tight')
        plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained');a=np.load(DATA/'chaotic.npz')
    axs[0].plot(a['state'][:,0],a['state'][:,2],lw=.5,color='#176d89');axs[0].set(xlabel='Lorenz X (model units)',ylabel='Lorenz Z (model units)',title='The added Lorenz drive at ρ = 28')
    axs[1].plot(a['time'],a['pressure'],color='#bb602f',lw=.8,label='Changing pressure');axs[1].plot(a['time'],a['state'][:,5],color='#176d89',lw=2,label='Fading imprint')
    axs[1].set(xlabel='Model time',ylabel='Model value',title='Memory filters the changing pressure');axs[1].legend();save(fig,'lorenz-pressure')
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for key,label,color in [('steady','Steady pressure','#697785'),('chaotic','Lorenz pressure + lean coupling','#176d89'),('no-feedback','Lorenz pressure; lean coupling off','#bb602f')]:
        a=np.load(DATA/('no-feedback-matched.npz' if key=='no-feedback' else key+'.npz'));axs[0].plot(a['time'],a['state'][:,3],label=label,color=color)
    axs[0].set(xlabel='Model time',ylabel='Lean q',title='The new drive changes the hug response');axs[0].legend(fontsize=8)
    for i in range(4):
        a=np.load(DATA/f'lyapunov-{i}.npz');axs[1].plot(a['time'],a['estimate'],lw=.8,label=f'Run {i+1}')
    axs[1].axhline(0,color='black',ls='--',lw=.8);axs[1].set(xlabel='Time after discarded transient',ylabel='Largest Lyapunov estimate (1/model time)',ylim=(-.2,1.3),title='Positive sensitivity persists under refinement');axs[1].legend(fontsize=8)
    save(fig,'lorenz-response')
    def img(n):return 'data:image/png;base64,'+base64.b64encode((FIG/(n+'.png')).read_bytes()).decode()
    style='''<style>:root{--foreground:#23363d;--muted-foreground:#5c6c74;--border:#c7d1d6;--viz-series-1:#16748b;--viz-series-2:#bc5d2e}body{font:17px/1.65 system-ui,sans-serif;margin:0;background:#f5f7f8;color:var(--foreground);overflow-wrap:anywhere}main{max-width:1050px;margin:auto;padding:32px 24px}section{background:white;padding:24px;border:1px solid #d9e0e4;border-radius:12px;margin:24px 0}h1{line-height:1.2}img{max-width:100%;height:auto}pre{white-space:pre-wrap;background:#edf2f4;padding:14px}.viz-controls{display:flex;gap:14px;flex-wrap:wrap;align-items:end}.form-label{display:block;min-width:0}.viz-controls>.form-label{flex:1 1 180px}.form-select{display:block;width:100%;max-width:100%;padding:8px;font:inherit}.btn{padding:9px 14px;font:inherit;background:#263c46;color:white;border:0;border-radius:7px}.form-range{display:block;width:100%}.text-small{font-size:14px}.tabular-nums{font-variant-numeric:tabular-nums}a{color:#176d89}</style>'''
    content=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Lorenz-driven hug</title>{style}<main>
<h1>Lorenz dynamics added to the hug</h1><p>The Lorenz system now drives changing pressure on the existing hug. The hug responds through its fading imprint, lean and original latched opening rule. This is a new, explicit one-way coupling.</p>
<section>{fragment}</section>
<p>Playback uses computed trajectories. The five modes include a steady input, chaotic Lorenz input, the same drive with lean coupling disabled, and stronger pressure that opens the arms. The stronger-pressure preset starts just above the opening limit, so opening begins at time zero. The square remains a reference guide.</p>
<section><h2>How the connection works</h2><pre>X′ = σ(Y−X)
Y′ = X(ρ−Z)−Y
Z′ = XY−βZ
pressure = mean + amplitude × tanh((Z−max(ρ−1,0))/10)
memory′ = (pressure−memory)/τ
lean″ = r × lean − lean³ − damping × lean′ + coupling × memory × lean</pre>
<p>Here σ=10 and β=8/3. Lorenz ρ is different from the hug's stiffness r. The bounded pressure mapping is an added modeling choice, not a measured relation. The hug does not feed back into the Lorenz system. Values use model units, and no permeability or active wall law is implied.</p>
<p>The two low-variation presets intentionally have the same mean pressure after their Lorenz transient, so they produce the same hug response despite different Lorenz coordinates. The chaotic input changes the time dependence. The coupling-off playback reuses the identical recorded Lorenz trajectory, pressure and memory, then solves the independent lean equation. This keeps the long-time control matched even when separate chaotic calculations drift apart. The original separately integrated control remains saved.</p></section>
<section><h2>Completed tests</h2><p><strong>84 new driven configurations</strong> span ρ=0.5, 10, 24, 24.74, 28, 100, 160; three damping values; coupling off/on; and two pressure amplitudes. All complete to t=20 with bounded pressure and imprint. The greatest scaled lean-energy residual is {summary['controls']['max_scaled_energy_residual']:.3e}.</p>
<p>All five playback presets pass independent Radau and smaller-step comparisons on [0,5]. Seven further short-horizon checks cover every tested ρ with low damping and active coupling; the largest scaled state discrepancy is {max(x['scaled_error'] for x in extra):.3e}, below their 1e-5 target. A separate opening event after startup occurs at t={event['opening_times'][0]:.8f}, agreeing across the three calculations within {event['maximum_difference']:.3e}.</p>
<p>These short-horizon comparisons check implementation accuracy. Long chaotic trajectories need not match point by point: no claim of accurate timing of every turn at t=40 is made. The 84-case sweep is a finite stress screen, not exhaustive parameter coverage or a measurement comparison.</p>
<img src="{img('lorenz-pressure')}" alt="Lorenz X-Z trajectory and the pressure-memory response">
<img src="{img('lorenz-response')}" alt="Lean response controls and four finite-time Lyapunov estimates">
<p>For the Lorenz drive at ρ=28, four finite-time largest Lyapunov estimates range from <strong>{summary['lyapunov_min']:.4f} to {summary['lyapunov_max']:.4f}</strong> per model time. Two starting states and two solver resolutions are used, discarding 50 time units and accumulating for 300. Positive values persist across these choices. This numerical evidence concerns the Lorenz drive; it does not prove an independently chaotic hug or a physical fluid attractor. The range is sensitivity across runs, not a statistical confidence interval.</p>
<p>Analytic fixed-point checks reproduce the Lorenz origin's loss of stability at ρ=1 and the nonzero equilibria's Hopf threshold ρ={summary['controls']['hopf_threshold']:.8f}. That threshold does not mean every larger parameter has identical chaotic behavior. The high-ρ stress cases are not labeled chaotic without separate diagnostics.</p></section>
<section><h2>Reproduce</h2><pre>python run_lorenz.py
python align_lorenz_control.py
python validate_lorenz_sweep.py
python build_lorenz_report.py</pre><p><a href="https://github.com/NousVolition/My-Sources-Project-ALL/tree/main/reports/hug-stress">Code and complete data on GitHub</a>. Lorenz-data state columns are X,Y,Z,q,v,imprint,work,dissipation. Each run saves the actual pressure too. Source hashes, tolerances, events and the finite-time sensitivity records are included.</p>
<p><a href="report.html">Earlier 224-case stress report</a> · <a href="https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml">Lorenz's original 1963 paper</a> · <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html">Numerical solver documentation</a></p></section></main></html>'''
    (ROOT/'lorenz-report.html').write_text(content,encoding='utf-8')
    print('Built Lorenz playback, report, and two scientific figures')
if __name__=='__main__':main()
