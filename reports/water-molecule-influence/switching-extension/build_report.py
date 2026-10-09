"""Build self-contained HTML, repository documentation, and a file manifest."""
from pathlib import Path
import base64, hashlib, json
from html import escape

HERE=Path(__file__).resolve().parent
D=HERE/'data'


def main():
    results=json.loads((D/'results.json').read_text())
    phase=json.loads((D/'phase_results.json').read_text())
    audit=json.loads((D/'audit.json').read_text())
    assert results['checks_passed'] and phase['checks_passed'] and audit['passed']
    rev,per=results['reversible'],results['persistent']
    off,on=[z['s'] for z in rev['folds']]
    short,long=per['pulse_runs'][:2]
    bottleneck=phase['periods'][-1]
    max_solver_error=max(r['max_DOP853_Radau_difference'] for r in results['convergence'])
    max_period_error=max(r['ODE_relative_error'] for r in phase['periods'])
    max_slip_error=max(r.get('largest_full_cycle_relative_error',0) for r in phase['locking'])
    mapping=[
        ('Biochemical switch','Activity x = Aₚ/K','Two stable states at the same stimulus; switching depends on history.'),
        ('Phase oscillator','Angle or phase difference on a circle','Rotation, locking and repeated slips; the basic locking model has one stable phase per cycle.'),
        ('Earlier water experiment','Response score of each molecule','Transient influence rankings under a chosen force field; these ODE tests do not establish water memory or synchrony.')]
    switch_summary=[
        ('Slow-sweep switch on, a=1.8',f's = {on:.6f}'),
        ('Slow-sweep switch off, a=1.8',f's = {off:.6f}'),
        ('Same-input stable activity levels',f'{rev["basin_finals"][0]:.6f} and {rev["basin_finals"][1]:.6f} at s={rev["basin_stimulus"]:.6f}'),
        ('Critical pulse duration, a=3, s=0.3',f'τ = {per["critical_pulse_duration_at_s03"]:.6f}'),
        ('Short pulse',f'τ={short["pulse_duration"]:.6f}: returns below 10⁻⁹'),
        ('Long pulse',f'τ={long["pulse_duration"]:.6f}: settles at x={long["x_after_100_without_stimulus"]:.6f}')]
    phase_summary=[
        ('Uniform runners, T₁=60 s and T₂=75 s','First lapping time: 300 s'),
        ('Nonuniform oscillator, ω=1, a=0','T = 6.283185'),
        ('Near bottleneck, ω=1, a=0.9999',f'T = {bottleneck["ODE_period"]:.6f}'),
        ('Time spent within θ=π/2 ± 0.2 rad',f'{100*bottleneck["fraction_time_in_04_rad_bottleneck"]:.2f}% of a cycle at a=0.9999'),
        ('Period divergence exponent',f'{phase["near_threshold_log_log_exponent"]:.6f} (predicted −1/2)'),
        ('Stable firefly phase locking','|Ω−ω| < A; marginal threshold at equality'),
        ('Full/reduced pendulum, μ=1.2','Maximum angle error falls from 0.299112 to 0.003528 rad as ε falls from 0.1 to 0.001')]
    def md_table(headers,rows):
        clean=lambda value:value.replace('|','\\|')
        return '| '+' | '.join(map(clean,headers))+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(clean,row))+' |' for row in rows)
    def table(headers,rows):
        return '<table><tr>'+''.join('<th>'+escape(h)+'</th>' for h in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+escape(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</table>'
    model='''The supplied equation is dAₚ/dt = kₚ S A + β Aₚⁿ/(Kⁿ+Aₚⁿ) − k_d Aₚ. The excerpt does not specify an equation or conservation law for A, so this implementation treats A as a fixed reservoir. Set x=Aₚ/K, τ=k_d t, s=kₚSA/(k_dK), and a=β/(k_dK). Then x′=s+a xⁿ/(1+xⁿ)−x. All reported switch parameters are illustrative and dimensionless; no concentrations, reaction rates or physical times were fitted.'''
    roots=f'''For n=2, equilibria solve x³−(s+a)x²+x−s=0, and stability follows from f′(x)=2ax/(1+x²)²−1. A fold additionally satisfies (1+x²)²=2ax. The feedback needed for any fold pair is a>8/(3√3)≈{results['critical_feedback_a_for_n2']:.6f}; whether both thresholds lie at nonnegative stimulus depends on a. At a=1.8 the two stable branches coexist for {off:.6f}<s<{on:.6f}. Starting on different sides of the unstable middle root at the same s={rev['basin_stimulus']:.6f} yields stable activities {rev['basin_finals'][0]:.6f} and {rev['basin_finals'][1]:.6f}.'''
    persistence=f'''At a=3, n=2 and zero stimulus, x=0 and x=(3+√5)/2≈{per['high_state']:.6f} are stable, separated by the unstable x=(3−√5)/2≈{per['separator']:.6f}. A rectangular s=0.3 pulse from x=0 crosses this boundary only if it lasts longer than τ*={per['critical_pulse_duration_at_s03']:.6f}. We compute τ* independently by integrating dx/f(x;s=0.3) up to the separator, then check it with ODE integration. Pulses of 0.95τ* and 1.05τ* land in opposite basins and remain distinct after 100 further units at s=0. The mathematical lower fold is at negative s; merely returning a nonnegative stimulus to zero therefore does not force the high state off. Turning feedback off resets the state, so this is persistence under fixed model conditions, not thermodynamic irreversibility.'''
    controls=f'''The sweep uses 601 equally spaced stimulus levels from 0 to 0.3 and back. For a=1.8, increasing the dwell per level from 5 to 25 to 125 gives loop areas 0.065049, 0.058563 and 0.056232, approaching the independently calculated quasistatic area {rev['static_loop_area']:.6f}. Finite-rate switching brackets differ slightly from exact folds because relaxation slows near the threshold; the numerical slowest loop still differs in area by about 1.17%. For a=1.4, below the fold threshold, the apparent loop shrinks to 1.33×10⁻⁹ at dwell 125. With feedback removed, the numerical relaxation agrees with its exact exponential solution and the loop shrinks toward zero. These controls distinguish slow response from bistability.'''
    n1='''A second control is deliberately subtle: changing the Hill exponent to n=1 while keeping a=3 still permits persistent high activity after a pulse, but x=0 is now unstable (f′(0)=2). Even a 10⁻⁶ seed grows to x=2 without a stimulus. There are not two stable zero-input states. A persistent trace alone is therefore insufficient evidence of bistable memory. In contrast, the same tiny seed decays to zero for the n=2 bistable case.'''
    oscillator='''Uniform phases obey θ̇₁=2π/T₁ and θ̇₂=2π/T₂. Their difference advances by 2π after T_lap=T₁T₂/(T₂−T₁), assuming T₂>T₁ and equal initial phases. For a nonuniform oscillator θ̇=ω−a sinθ with 0≤a<ω, a full turn takes T=2π/√(ω²−a²). Expanding near a=ω from below gives T≈π√(2/ω)(ω−a)⁻¹ᐟ². The bottleneck is near θ=π/2, where speed is smallest. At the threshold the period diverges and a degenerate equilibrium appears; beyond it, fixed phases replace positive circulation.'''
    locking='''For the supplied firefly model Θ̇=Ω and θ̇=ω+A sin(Θ−θ), define φ=Θ−θ, τ=At and μ=(Ω−ω)/A, with A>0. Then φ′=μ−sinφ. When |μ|<1, the stable phase is arcsinμ modulo 2π; the other equilibrium is unstable. Phase locking means a constant phase difference, which need not be zero. Equality |μ|=1 is a degenerate threshold, not robust exponential locking. For |μ|>1, the phase slips repeatedly, with period 2π/√(μ²−1). We measure complete cycles; an average over a short window containing fractional cycles can differ noticeably from the asymptotic drift. Exact unstable initial equilibria are exceptions to attraction.'''
    pendulum='''The mechanical equation is Iθ̈+bθ̇+mgL sinθ=Γ with I=mL². Using τ=(mgL/b)t, μ=Γ/(mgL) and ε=I mgL/b² gives εθ″+θ′=μ−sinθ. The first-order limit is the same phase equation. We compare ε=0.1, 0.01 and 0.001 at μ=0.8 and 1.2 from θ=0 and zero physical angular velocity, on 0≤τ≤40. The full model has an initial velocity-relaxation layer. Decreasing ε improves the angle agreement in both tested cases; this does not justify dropping inertia for arbitrary damping, torque or initial velocity.'''
    topology='''The shared mechanism is sensitivity near a threshold: a saddle-node can destroy an equilibrium and leave a slow bottleneck. This can delay a switch or greatly lengthen a rotation period. The outcomes still differ. A first-order autonomous flow on a line cannot have a nonconstant periodic orbit; a phase variable lives on a circle and can circulate continuously. A scalar activity can have two stable equilibria, while the basic first-order phase-locking equation has one stable fixed phase per cycle. Neither observation identifies a permanent leading molecule.'''
    validation=f'''DOP853 integrations use rtol=10⁻¹⁰ and atol=10⁻¹² for switching; pulse discontinuities are integrated as separate segments. Radau with tighter tolerances agrees with the two near-threshold pulse traces to within {max_solver_error:.2e} in x. Independent fixed-step RK4 endpoint errors decrease on halving the step and remain below 3×10⁻¹⁰ in these tests. ODE rotation periods agree with analytic values to relative error at most {max_period_error:.2e}; complete phase-slip intervals agree to {max_slip_error:.2e}. These are observed numerical errors for the recorded parameter choices, not general solver guarantees. The pendulum uses Radau for its stiff velocity dynamics. All saved arrays were checked for finite values; loop areas and final pulse states were reconstructed from them.'''
    limits='''These are deterministic phenomenological ODE tests. A clamped A pool does not enforce conservation of A+Aₚ, and the Hill term is not a full phosphorylation reaction mechanism. Noise, finite molecule numbers, ATP supply, spatial transport and measured kinetics are absent. Persistent states may switch under noise or changed parameters. The magnetization excerpt concerns collective order, a different observable from individual influence rank. None of the present simulations demonstrates phosphorylation, hysteresis or synchronization in the earlier rigid-water system. Useful next steps are a specified conserved-pool reaction model, calibrated rates, stochastic switching-lifetime tests, and sensitivity to parameter and sweep-rate uncertainty.'''
    commands='''python -m pip install -r switching-extension/requirements.txt
python switching-extension/switching_test.py
python switching-extension/phase_tests.py
python switching-extension/make_figures.py
python switching-extension/build_report.py'''
    readme=f'''# Reproducible switching, hysteresis, and phase-locking tests

Tests the supplied biochemical positive-feedback equation and the later oscillator, pendulum and firefly examples. **These are illustrative ODE models, distinct from the earlier water molecular dynamics.**

{md_table(['Model','State measured','What was tested'],mapping)}

## Hysteresis and pulse switching

![Switching results](switching_results.png)

{md_table(['Test','Recorded result'],switch_summary)}

### Model and assumptions

{model}

{roots}

{persistence}

### Controls against false memory

{controls}

{n1}

## Rotation and synchronization

![Phase results](phase_results.png)

{md_table(['Test','Recorded result'],phase_summary)}

{oscillator}

{locking}

{pendulum}

### How the examples connect

{topology}

## Reproduce

Use Python 3.12. From the parent experiment directory:

```sh
{commands}
```

Both simulation scripts accept `--out PATH` to save a new run separately. Figure and report scripts read this folder’s `data/` by default. The scripts are deterministic; no random seed is required. The report labels use x=Aₚ/K and τ=k_d t, not physical concentration or seconds, except the explicitly labeled runner example.

{validation}

[Full self-contained report](report.html) · [Switch protocol](data/protocol.json) · [Switch measurements](data/results.json) · [Phase protocol](data/phase_protocol.json) · [Phase measurements](data/phase_results.json) · [Saved switching trajectories](data/trajectories.npz) · [Saved phase trajectories](data/phase_trajectories.npz) · [Checksums](manifest_sha256.json)

## Limitations and next validation

{limits}

## Sources and provenance

The model equations and examples come from the user-supplied textbook excerpts. The choices of fixed A, parameters, numerical protocols and new figures are explicit additions in this experiment; the cropped original exercises are not claimed to be fully solved. The scanned pages are not redistributed. For the broader distinction between positive feedback and verified bistability, see the primary research paper [Angeli, Ferrell and Sontag (2004)](https://pmc.ncbi.nlm.nih.gov/articles/PMC357011/).

[Original water experiment](../README.md) · [Damping and rotating-hoop tests](../damping-extension/README.md)
'''
    (HERE/'README.md').write_text(readme,encoding='utf-8')
    def picture(name):
        return 'data:image/png;base64,'+base64.b64encode((HERE/name).read_bytes()).decode()
    model,roots,persistence,controls,n1,oscillator,locking,pendulum,topology,validation,limits = map(escape,
        [model,roots,persistence,controls,n1,oscillator,locking,pendulum,topology,validation,limits])
    html=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Switching, hysteresis and synchronization</title>
<style>body{{font:17px/1.65 system-ui,sans-serif;max-width:1100px;margin:45px auto;padding:0 24px;color:#17343d;background:#fafcfd}}h1{{font-size:40px;line-height:1.15}}h2{{font-size:26px;margin-top:40px}}h3{{font-size:21px}}.lead{{font-size:21px}}.box{{padding:18px;background:#e7f2f3;border-left:4px solid #087e83}}img{{width:100%;margin:18px 0}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #ccdbe0}}th{{background:#e7eff2}}pre{{padding:18px;overflow:auto;background:#ecf2f4}}a{{color:#087e83}}</style>
<p>REPRODUCIBLE MODEL TESTS · 9 OCTOBER 2026</p><h1>When a signal leaves a memory.<br>When an oscillator locks its phase.</h1>
<p class="lead">A temporary pulse can place a feedback model in a different stable state. A driven oscillator can settle to a fixed phase difference or slip repeatedly. These experiments distinguish the mechanisms and test the equations in the supplied pages.</p>
<div class="box">Model scope: deterministic concentration and phase equations with illustrative parameters. These calculations do not establish memory, phosphorylation or synchronization in the earlier water simulations.</div>
{table(['Model','State measured','Meaning'],mapping)}
<h2>1. A biochemical switch with history dependence</h2><img src="{picture('switching_results.png')}" alt="Hysteresis branches, pulse-triggered persistent activation, sweep-rate controls and two stable outcomes at identical input">
{table(['Test','Result'],switch_summary)}<h3>The equation actually simulated</h3><p>{model}</p><div class="box">x′ = s + a xⁿ/(1+xⁿ) − x<br>x=Aₚ/K, τ=k_d t, s=kₚSA/(k_dK), a=β/(k_dK).</div><p>{roots}</p><p>{persistence}</p>
<h3>Controls: persistence alone is not enough</h3><p>{controls}</p><p>{n1}</p>
<h2>2. Uniform and nonuniform oscillators</h2><img src="{picture('phase_results.png')}" alt="Period divergence, bottleneck time, locked and slipping phases, and convergence of the inertial pendulum toward its first-order limit">
{table(['Test','Result'],phase_summary)}<p>{oscillator}</p><h3>Firefly entrainment</h3><p>{locking}</p><div class="box">φ′ = μ − sinφ<br>Stable locking: |μ|&lt;1. Stable offset: φ*=arcsinμ (mod 2π).<br>Slipping: |μ|&gt;1, with T_slip=2π/√(μ²−1).</div><h3>The overdamped pendulum</h3><p>{pendulum}</p>
<h2>3. What the examples do—and do not—share</h2><p>{topology}</p><p>{limits}</p>
<h2>Methods and numerical checks</h2><p>{validation}</p><p>Transient switch traces run for 100 scaled time units after stimulus removal. The same-input basin test runs for 300 units. Phase locking is followed for 300 units. The reported divergence slope uses the exact period over 10⁻⁷≤1−a≤10⁻³ and is supported by separately integrated periods; it is not a fit to physical observations. No uncertainty intervals based on independent experimental samples are claimed.</p>
<h2>Reproduce</h2><p>From the parent experiment directory, with Python 3.12:</p><pre>{commands}</pre><p>Saved trajectories, equilibrium checks, protocols and SHA-256 checksums accompany this self-contained report. Use --out with either simulation script for a separate run; the figure and report builders read the delivered data folder by default.</p>
<h2>Sources</h2><p>Equations follow the supplied textbook excerpts; the fixed-pool assumption and numerical parameter choices are specified above. For primary research on testing bistability and hysteresis in feedback systems, see <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC357011/">Angeli, Ferrell and Sontag (2004)</a>. All displayed figures were generated from these recorded calculations.</p></html>'''
    (HERE/'report.html').write_text(html,encoding='utf-8')
    (HERE/'requirements.txt').write_text('numpy==2.5.3\nscipy==1.18.1\nmatplotlib==3.11.2\n',encoding='utf-8')
    files={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.rglob('*')
           if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest_sha256.json'}
    (HERE/'manifest_sha256.json').write_text(json.dumps(files,indent=2),encoding='utf-8')
    print(json.dumps({'files':len(files)+1,'switch_on':on,'switch_off':off,
                      'critical_pulse_duration':per['critical_pulse_duration_at_s03'],'checks_passed':True},indent=2))


if __name__=='__main__':main()
