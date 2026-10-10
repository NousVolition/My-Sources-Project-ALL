"""Build a matched comparison from completed simulations; never invent trajectories."""
from pathlib import Path
import json,base64
from dataclasses import asdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import two_way_hug as tw
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'two-way-data';FIG=ROOT/'figures'
def read(name):return json.loads((DATA/name).read_text())
def main():
    summary=read('summary.json');paired=read('paired.json');opening=read('opening-tests.json');controls=read('controls.json')
    assert summary['controls_passed'] and opening['all_passed']
    records={};cases=[]
    def record(name,g,opened,stride):
        a=np.load(DATA/(name+'.npz'));t=a['time'][::stride];y=a['state'][::stride];p=a['pressure'][::stride]
        eta=[tw.resistance(tt,mm,g,opened) for tt,mm in zip(t,y[:,5])]
        bends=[tw.bend(tt,opened) for tt in t]
        records[name]=dict(samples=np.round(np.column_stack([t,y[:,0],y[:,2],y[:,3],y[:,5],p,eta,bends]),6).tolist())
    for r in paired:record(r['name'],r['strength'],r['metrics']['opening'],5)
    for seed in range(3):
        for g in [.5,2.,8.]:
            cases.append(dict(key=f'closed-{seed}-{g:g}',label=f'Closed · return {g:g} · start {seed+1}',one=f'paired-{seed}-g0',two=f'paired-{seed}-g{g:g}'))
    for r in opening['cases']:record(r['name'],r['strength'],r['opening'],10)
    for seed in range(2):
        for pressure in [.99,1.,1.01]:
            items=[r for r in opening['cases'] if r['seed']==seed and r['requested_initial_pressure']==pressure]
            cases.append(dict(key=f'opening-{seed}-{pressure}',label=f'Opening · initial pressure {pressure:.2f} · start {seed+1}',one=items[0]['name'],two=items[1]['name']))
    payload=dict(cases=cases,records=records)
    fragment=(ROOT/'two-way-template.html').read_text(encoding='utf-8').replace('__TWO_WAY_DATA__',json.dumps(payload,separators=(',',':')))
    assert len(fragment.encode())<1_000_000
    (ROOT/'two-way-fragment.html').write_text(fragment,encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for g,color,label in [(0,'#a66439','One-way'),(8,'#087d88','Two-way, return 8')]:
        a=np.load(DATA/f'paired-0-g{g}.npz');t=a['time'];y=a['state']
        axs[0,0].plot(y[:,0],y[:,2],color=color,lw=.6,label=label)
        axs[0,1].plot(t,a['pressure'],color=color,lw=.8,label=label)
        axs[1,0].plot(t,y[:,5],color=color,lw=1.1,label=label)
        axs[1,1].plot(t,y[:,3],color=color,lw=1.1,label=label)
    axs[0,0].set(xlabel='Lorenz X',ylabel='Lorenz Z',title='The return resistance changes the drive')
    for ax,title,ylabel in [(axs[0,1],'Changed drive returns different pressure','Pressure'),(axs[1,0],'The hug retains a different imprint','Imprint'),(axs[1,1],'The lean changes in response','Lean q')]:
        ax.set(xlabel='Model time',ylabel=ylabel,title=title)
    axs[0,0].legend();fig.suptitle('Matched starting state · all values in model units')
    fig.savefig(FIG/'two-way-response.png',dpi=160);plt.close(fig)
    labels=['Lorenz activity RMS','Z fluctuations (SD)','Mean pressure','Mean imprint','Lean RMS']
    keys=['activity_rms','z_std','mean_pressure','mean_imprint','lean_rms']
    fig,ax=plt.subplots(figsize=(11,5),layout='constrained')
    for j,s in enumerate(summary['effects']):
        vals=[s['effects'][k]['mean_percent'] for k in keys]
        ranges=[s['effects'][k]['percent_changes']+s['effects'][k]['refined_percent_changes'] for k in keys]
        err=np.array([[v-min(r) for v,r in zip(vals,ranges)],[max(r)-v for v,r in zip(vals,ranges)]])
        ax.errorbar(np.arange(5)+(j-1)*.17,vals,yerr=err,fmt='o',capsize=3,label=f'Return {s["strength"]:g}')
    ax.axhline(0,color='#65757a',lw=.7);ax.set(xticks=np.arange(5),xticklabels=labels,ylabel='Change from matched one-way run (%)',title='Three starting states and two solver resolutions · time 20–80')
    ax.legend();fig.savefig(FIG/'two-way-effects.png',dpi=160);plt.close(fig)
    rows=''.join('<tr><td>'+str(s['strength'])+'</td>'+''.join(f'<td>{s["effects"][k]["mean_percent"]:+.2f}%</td>' for k in keys)+'</tr>' for s in summary['effects'])
    img=lambda name:'data:image/png;base64,'+base64.b64encode((FIG/name).read_bytes()).decode()
    max_error=max(c['independent_error'] for c in controls)
    report=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Two-way hug and Lorenz tests</title>
<style>:root{{--foreground:#23363d;--background:#fff;--border:#c7d1d6;--viz-series-1:#087d88;--viz-series-2:#a66439}}body{{font:16px/1.6 system-ui,sans-serif;max-width:1080px;margin:auto;padding:24px;color:var(--foreground)}}img{{max-width:100%;height:auto}}.viz-controls{{display:flex;gap:16px;flex-wrap:wrap;align-items:end}}.form-label{{display:block;min-width:0}}.form-select{{display:block;max-width:100%;padding:8px;font:inherit}}.form-range{{display:block;width:100%}}.btn{{font:inherit;padding:10px}}.viz-row{{display:flex;gap:20px;flex-wrap:wrap}}.table-responsive{{overflow:auto}}table{{border-collapse:collapse;width:100%;margin:16px 0}}th,td{{text-align:left;padding:7px;border-bottom:1px solid #ddd}}pre{{white-space:pre-wrap;background:#f0f4f5;padding:16px}}.text-small{{font-size:14px}}a{{color:#087080}}</style>
<h1>The hug and Lorenz drive affect each other</h1>
<p>This experiment adds a proposed return resistance to the original one-way model. Lorenz motion changes pressure, pressure changes the hug's imprint and lean, and the imprint and opening state change the resistance applied back to Lorenz X and Y. No physical calibration is claimed.</p>
<p><strong>Meaning of the hug:</strong> the user's description of a breathable form means a clay-like form that yields and molds easily around what presses into it. It is not a request for a valve or a repeated opening-and-closing cycle. The inherited arm-opening threshold is a feature of this simplified implementation, not a defining property of that idea. These tests explore reciprocal coupling in the existing reduced equations; they do not yet implement clay-like contact or general reshaping.</p>
{fragment}
<h2>Completed comparisons</h2><p>12 matched configurations use three distinct starting states and four return strengths, each run at two solver resolutions to time 80. Statistics below discard time 0–20. Another 12 configurations test Lorenz parameters 10, 28 and 100, hug damping 0.05 and 50, and return strengths 2 and 20 to time 20. Seven short control configurations use DOP853, Radau and a refined DOP853 run, including return strength 100 and hug damping 500. A constant-pressure control uses two Lorenz starting states. Twelve opening configurations each use two resolutions.</p>
<p><strong>83 integrations completed:</strong> 24 paired, 12 stress, 21 short controls, 2 constant-pressure, and 24 opening runs; plus one original one-way integration for the zero-return comparison. All checked bounds and passivity conditions pass. These counts refer to computations, not independent physical measurements.</p>
<div class="table-responsive"><table><thead><tr><th>Return strength</th>{''.join('<th>'+l+'</th>' for l in labels)}</tr></thead><tbody>{rows}</tbody></table></div>
<p>Entries are mean paired percentage changes over three starts at the base solver resolution. Lorenz activity RMS is sqrt(mean(X²+Y²)); fluctuations are the time standard deviation of Z. RMS lean is sqrt(mean(q²)). The intervals in the figure span all three starts at both resolutions; they are sensitivity ranges, not confidence intervals.</p>
<img src="{img('two-way-effects.png')}" alt="Percentage changes and sensitivity ranges for activity, fluctuations, pressure, imprint and lean at three feedback strengths">
<p>At return strength 8, the Z fluctuations almost disappear over the measured window, but nonzero Lorenz motion remains. Imprint decreases about 36% and lean RMS about 31%. At strength 2, mean pressure rises while imprint falls: the imprint charges rapidly and fades slowly, so its response depends on the pressure history, not just average pressure. Weak-return differences are small and sensitive to starting state and numerical phase drift. These finite-time calculations do not establish a global stability or chaos theorem.</p>
<img src="{img('two-way-response.png')}" alt="Matched one-way and two-way trajectories, pressure, imprint and lean">
<h2>What was added</h2><pre>eta(t) = g × bend(t)² × imprint / (1 + imprint)
X′ = sigma(Y − X) − eta X
Y′ = X(rho − Z) − Y − eta Y
Z′ = XY − beta Z
pressure = mean + amplitude tanh((Z − max(rho − 1, 0))/10)
imprint′ = (pressure − imprint)/tau
tau = 0.8 when charging; 6 when fading
q″ = r q − q³ − damping q′ + memory_coupling imprint q</pre>
<p>The remaining hug parameters and opening geometry are inherited. Lean is an output of the coupling; the added return force depends on imprint and opening, not directly on lean. The resistance rule is a modeling proposal used to test two-way influence, not the user's definition of moldability. The square is still a reference guide. This closure does not model clay-like contact, permeability, a wall force, quantum effects or resolved Navier–Stokes flow.</p>
<h2>Opening is a one-time release in this version</h2><p>At the first upward crossing of pressure 1, the solver stops exactly at the event and resumes with the existing smooth opening law. Bend falls from 1 to 0 over 0.8 model time. It stays at 0 afterward, so return resistance vanishes and stays absent even if pressure drops. The arms can still lean and retain imprint. This imposed opening is not a model of clay-like molding.</p>
<p>The 12 opening tests start at pressure 0.99, 1.00 and 1.01, with two different Lorenz configurations and return off/on. The largest state discrepancy after reducing the maximum step fourfold and tightening tolerance is {opening['maximum_error']:.3e}; opening times differ by at most {opening['max_event_error']:.3e}. Every completed release has zero return resistance. This verifies the imposed rule; it does not validate permanent opening as a property of the user's idea.</p>
<h2>Numerical checks and limits</h2><p>The largest short independent-method scaled state discrepancy is {max_error:.3e}, below the preselected 1e-5 target. Turning return strength to zero reproduces the original one-way system within {controls[0]['original_one_way_error']:.3e}. With pressure modulation disabled, changing the Lorenz initial state changes hug q, velocity and imprint by at most {summary['constant_pressure_hug_independence_error']:.3e}.</p>
<p>The new term contributes −eta(X²+Y²) to the derivative of (X²+Y²+Z²)/2, so it only removes this quadratic activity. The maximum normalized balance residual is {summary['maximum_scaled_budget_residual']:.3e}. This is a diagnostic of the model equations, not a measured physical energy transfer: the removed activity is not deposited in a material-energy equation for the hug.</p>
<p>Long chaotic paths may separate under refinement. The reported changes compare time-window statistics, with the spread from starting states and refinement retained. There are only three long-run starts, no statistical confidence intervals, no temperature or material parameters, and no measured reciprocal-coupling benchmark. The new rule is a testable mathematical proposal, not a verified law of your hug.</p>
<h2>Reproduce</h2><pre>python run_two_way.py
python test_two_way_opening.py
python build_two_way_report.py
python verify_two_way.py</pre><p>The original one-way sources and results are retained. <a href="lorenz-report.html">Original Lorenz report</a> · <a href="README.md">Package overview</a> · <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html">Solver and event documentation</a></p></html>'''
    (ROOT/'two-way-report.html').write_text(report,encoding='utf-8')
    print(json.dumps(dict(comparisons=len(cases),records=len(records),fragment_bytes=len(fragment.encode()))))
if __name__=='__main__':main()
