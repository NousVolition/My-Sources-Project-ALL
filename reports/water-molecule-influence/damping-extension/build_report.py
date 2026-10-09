"""Reanalyze saved damping arrays and build figures, HTML, README and hashes."""
from pathlib import Path
import base64, hashlib, json
from itertools import cycle
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE/'data'


def read(name):
    return json.loads((DATA/name).read_text())


def main():
    oscillator, validation, hoop, protocol, results = [read(n) for n in ['oscillator.json', 'integrator_validation.json', 'hoop.json', 'protocol.json', 'results.json']]
    runs = results['runs']
    assert len(runs) == len(protocol['gamma_ps_inverse'])
    times = np.array(protocol['times_ps'])*1000
    osc, hp = np.load(DATA/'oscillator.npz'), np.load(DATA/'hoop.npz')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    fig.suptitle('What changes when damping is added?', fontsize=18, fontweight='bold')
    for row in oscillator['rows']:
        key = f'g{row["gamma"]:g}'
        ax[0, 0].plot(osc[key+'_time'], osc[key+'_position'], label=f'{row["regime"]}: γ={row["gamma"]:g}')
    ax[0, 0].axhline(0, color='#999', lw=.7)
    ax[0, 0].set(title='A | Spring: weak damping still permits oscillation', xlabel='Dimensionless time', ylabel='Displacement')
    ax[0, 0].legend(fontsize=8, frameon=False)
    audit = []
    for color, run in zip(cycle(['#087e83', '#cc642c', '#4964a8', '#915da8']), runs):
        gamma = run['gamma_ps_inverse']
        saved = np.load(DATA/f'water_gamma_{gamma:g}.npz')
        scores, pairs = saved['scores_ps'], saved['pair_response_norms_ps']
        assert scores.shape == (216, 4)
        assert np.isfinite(pairs).all() and (pairs >= 0).all()
        assert np.all(pairs[np.arange(216), :, np.arange(216)] == 0)
        assert np.allclose(scores, pairs.sum(axis=2), rtol=1e-14, atol=1e-14)
        for k, row in enumerate(run['rows']):
            assert np.isclose(row['mean_score_ps'], scores[:, k].mean(), rtol=1e-12)
            assert np.isclose(row['spearman_vs_undamped'], spearmanr(scores[:, k], saved['undamped_scores_ps'][:, k]).statistic, rtol=1e-12)
        label=f'γ={gamma:g} ps⁻¹'
        ax[0, 1].plot(times, [r['mean_over_undamped'] for r in run['rows']], '-o', color=color, label=label)
        ax[1, 0].plot(times, [r['spearman_vs_undamped'] for r in run['rows']], '-o', color=color, label=label)
        ax[1, 1].plot(np.r_[0, times], run['energy_change_kJ_mol_per_water'], '-o', color=color, label=label)
        audit.append({'gamma': gamma, 'saved_array_checks_passed': True})
    ax[0, 1].axhline(1, color='#999', ls=':')
    ax[0, 1].set(title='B | Water: mean response relative to no damping', xlabel='Time after impulse (fs)', ylabel='Damped / undamped mean score')
    ax[1, 0].set(title='C | Water: does the undamped ranking survive?', xlabel='Time after impulse (fs)', ylabel='Spearman rank correlation', ylim=(-1.05, 1.05))
    ax[1, 1].set(title='D | Water: energy removed by friction', xlabel='Time (fs)', ylabel='Energy change per water (kJ/mol)')
    for a in [ax[0, 1], ax[1, 0], ax[1, 1]]:
        a.legend(frameon=False)
    for a in ax.flat:
        a.spines[['top', 'right']].set_visible(False)
    fig.savefig(HERE/'damping_results.png', dpi=170)
    fig.savefig(HERE/'damping_results.svg')
    plt.close(fig)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    fig.suptitle('The rotating hoop: equilibria and fast–slow motion', fontsize=18, fontweight='bold')
    below, above = np.linspace(0, 1, 100), np.linspace(1, 3, 200)
    ax[0, 0].plot(below, below*0, color='#087e83', lw=2, label='Stable')
    ax[0, 0].plot(above, above*0, '--', color='#cc642c', label='Unstable')
    ax[0, 0].plot(above, np.arccos(1/above), color='#087e83', lw=2)
    ax[0, 0].plot(above, -np.arccos(1/above), color='#087e83', lw=2)
    ax[0, 0].set(title='A | Rotation changes the stable equilibria', xlabel='q = rω²/g', ylabel='Equilibrium angle φ (rad)')
    ax[0, 0].legend(frameon=False)
    for beta, color in [(2., '#cc642c'), (20., '#087e83')]:
        key = f'q2_b{beta:g}'
        ax[0, 1].plot(hp[key+'_tau'], hp[key+'_phi'], color=color, label=f'Full inertia, β={beta:g}')
    key = 'q2_b20'
    ax[0, 1].plot(hp[key+'_tau'], hp[key+'_overdamped'], 'k--', lw=1.2, label='First-order limit')
    ax[0, 1].set(title='B | Stronger damping approaches the reduction', xlabel='Slow time τ = t/T', ylabel='Angle φ (rad)')
    ax[0, 1].legend(frameon=False, fontsize=8)
    xx = np.linspace(-2*np.pi, 2*np.pi, 1000)
    for a, beta, label in [(ax[1, 0], 2., 'C'), (ax[1, 1], 20., 'D')]:
        a.plot(xx, np.sin(xx)*(2*np.cos(xx)-1), 'k--', lw=1.2, label='C: Ω = f(φ)')
        for i in range(6):
            key=f'phase_b{beta:g}_i{i}'
            x, y = hp[key+'_phi'], hp[key+'_omega']
            a.plot(x, y, lw=1.1, alpha=.8)
            a.plot(x[0], y[0], 'o', ms=3, color='#555')
            index = min(20, len(x)-2)
            a.annotate('', xy=(x[index+1], y[index+1]), xytext=(x[index], y[index]), arrowprops={'arrowstyle':'->', 'color':'#555'})
        a.set(title=f'{label} | β={beta:g}, ε={1/beta**2:g}: relaxation toward C', xlabel='Angle φ (rad)', ylabel='Scaled angular velocity Ω', xlim=(-2*np.pi, 2*np.pi), ylim=(-3, 3))
        a.legend(frameon=False, fontsize=8, loc='upper right')
    for a in ax.flat:
        a.spines[['top', 'right']].set_visible(False)
    fig.savefig(HERE/'hoop_results.png', dpi=170)
    fig.savefig(HERE/'hoop_results.svg')
    plt.close(fig)

    headers = ['Drag γ (ps⁻¹)', 'Mean response / undamped at 200 fs', 'Rank correlation with undamped', 'Leading molecule', 'Original force leader rank']
    rows = [[f'{r["gamma_ps_inverse"]:g}', f'{r["rows"][-1]["mean_over_undamped"]:.4f}', f'{r["rows"][-1]["spearman_vs_undamped"]:.3f}', str(r['rows'][-1]['leader']), str(r['rows'][-1]['force_leader_rank'])] for r in runs]
    md_table = '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(r)+' |' for r in rows)
    html_table = '<table><tr>'+''.join(f'<th>{h}</th>' for h in headers)+'</tr>'+''.join('<tr>'+''.join(f'<td>{c}</td>' for c in r)+'</tr>' for r in rows)+'</table>'
    checks=[]
    for r in runs:
        checks.append(f'γ={r["gamma_ps_inverse"]:g}: half-step L2 difference at 200 fs {100*r["half_timestep"]["L2_by_time"][-1]:.4f}%; at 5 fs {100*r["half_timestep"]["L2_by_time"][0]:.2f}%. Half-kick L2 difference at 200 fs {100*r["half_kick"]["L2_by_time"][-1]:.4f}%; maximum individual half-kick change over all checked times {100*max(r["half_kick"]["max_individual_relative_by_time"]):.2f}%. No-interaction residual {r["no_interactions_max_score_ps"]:.2e} ps; largest per-time permutation L2 error {max(r["permutation"]["L2_by_time"]):.2e}.')
    hoop_pair=[r for r in hoop['rows'] if r['q']==2]
    command='''python damping-extension/damped_leapfrog.py
python damping-extension/rotating_hoop.py
python damping-extension/validate_integrator.py
python damping-extension/damped_water.py --gamma 5 50 --platform OpenCL
python damping-extension/build_report.py'''
    gammas = ', '.join(f'{g:g}' for g in protocol['gamma_ps_inverse'])
    controls_count = len(protocol['controls_source_ids'])
    methods=f'''The water experiment uses one saved 216-molecule TIP3P state, with exactly the same positions and velocities as run {protocol['replica']+1} of the preceding impulse study. It tests every source molecule with paired ±0.01 nm/ps kicks on all three axes, at γ={gammas} ps⁻¹: {216*6*len(runs):,} primary perturbed trajectories plus controls. Each trajectory lasts 200 fs, with 1 fs time steps and 10⁻⁸ rigid-constraint tolerance. The no-damping scores come from the preceding experiment and are rechecked on {controls_count} selected sources. Friction acts equally on every atom as −mγv. No random force or thermal bath is added: positive damping deliberately dissipates energy and is not equilibrium water at 300 K.'''
    scheme='''The literal leapfrog code stores velocity at half time steps. With h=Δt, a=(1−γh/2)/(1+γh/2), c=h/(1+γh/2), its update is v[n+1/2]=a v[n−1/2]+c F(x[n])/m; x[n+1]=x[n]+h v[n+1/2]. Initialize v[−1/2]=(1+γh/2)v[0]−hF(x[0])/(2m). This centers the drag force and is second order. A large time step can still introduce numerical oscillations: γh>2 makes a negative, so damping is not a replacement for time-step validation.'''
    water_scheme='''The water code retains the previously corrected, synchronized velocity-Verlet/RATTLE integrator, rather than reintroducing mismatched half-step velocities. It adds exact half drag steps v←exp(−γh/2)v before and after each constrained Verlet step. Both implementations approximate m ẍ=F−mγẋ. They are separate second-order schemes, not asserted to be identical at finite h. At γ=0 the water integrator reproduces the previous undamped calculation.'''
    correction='''The compensating kicks also decay, so subtracting the old ballistic displacement t would create a false response. We replace t by Gγ,h(t)=h exp(−γh/2)[1−exp(−γt)]/[1−exp(−γh)], the exact free-translation response of this split integrator. Its continuous-time limit is [1−exp(−γt)]/γ; its γ=0 limit is t. The reported score remains the sum of Frobenius norms of other molecules’ corrected COM derivatives, in ps; it is not an energy-transfer fraction.'''
    scaling='''Using the book’s tangential-force convention, mr φ̈+b φ̇=mg sinφ(q cosφ−1), with q=rω²/g. Choose T=b/(mg), τ=t/T and ε=m²gr/b²=1/β², where β=b/(m√(gr)). Then ε φ″+φ′=f(φ), f(φ)=sinφ(q cosφ−1). With Ω=φ′, the phase equations are φ′=Ω and εΩ′=f(φ)−Ω. The curve C: Ω=f(φ) is the zero-angular-acceleration curve (a nullcline); it is not exactly invariant for finite ε. Small ε produces rapid velocity relaxation followed by slower angular motion near C, away from problematic scalings and after the initial transient.'''
    equilibrium='''For positive damping, the bottom φ=0 is stable for q<1 and unstable for q>1; the top φ=π is unstable. Two stable equilibria ±arccos(1/q) appear for q>1 (angles modulo 2π). At q=1 the bottom is nonhyperbolic and nonlinearly stable. Damping does not move these equilibria: it changes the approach to them. The symmetric branches form a supercritical pitchfork. These are properties of the hoop model, not evidence of a phase transition in the water experiment.'''
    limits='''One conditional water microstate cannot establish ensemble behavior. The two drag rates are sensitivity tests; “stronger damping” does not prove every molecular mode is overdamped. They correspond to free-velocity decay times of 200 and 20 fs, but interacting motions have additional timescales. Four selected sources receive half-step, half-kick, no-interaction and permutation controls; full-ranking convergence remains untested. Saved norm matrices support score reanalysis, but signed response blocks and every perturbed trajectory frame are not retained. Next checks are more independent states, full-ranking convergence, a range of mode-dependent damping scales, and—if equilibrium liquid water is the goal—Langevin noise matched to friction, with coupled noise for paired trajectories.'''
    readme=f'''# Damping added to leapfrog and the water-response experiment

Implements velocity-proportional drag, validates it against an exact spring solution, and follows the supplied rotating-hoop pages through equilibrium, scaling and fast–slow phase-space tests.

![Damping results](damping_results.png)

## Water results

{methods}

{md_table}

The leading ID is zero based; ranks run from 1 to 216. These are model-dependent responses, not permanent leaders.

## Two clearly identified numerical methods

{scheme}

{water_scheme}

{correction}

## The rotating hoop from the supplied pages

![Rotating hoop](hoop_results.png)

{scaling}

{equilibrium}

At q=2 and the same initial angle 0.2 rad, zero physical angular velocity, the maximum angle disagreement with the first-order model over slow time 0–12 is {hoop_pair[0]['max_angle_error_vs_overdamped']:.6f} rad at β=2, versus {hoop_pair[1]['max_angle_error_vs_overdamped']:.6f} rad at β=20. Halving the latter simulation’s step changes the full second-order angle by at most {hoop['half_step_angle_difference_q2_beta20']:.2e} rad. The phase portraits use six additional initial states at q=2.

## Validation

- Spring tests cover undamped, underdamped, critical and overdamped cases. Halving the step reduces position error approximately fourfold for both implementations.
- The actual OpenMM integrator reproduces free exponential velocity decay and its discrete displacement to within {max(validation['free_drag_discrete_position_error_nm'], validation['free_drag_velocity_error_nm_ps']):.2e} in the tested numerical units.
- {' '.join(checks)}

## Reproduce

Install the parent experiment’s pinned `requirements.txt`; run from `reports/water-molecule-influence/`:

```sh
{command}
```

Use `--platform Reference` for slower portable water execution. `damped_water.py` accepts `--gamma`, `--replica` (0–2), and `--out`. A separate `--out` avoids replacing the delivered water data; report generation reads this folder’s `data/` by default. The literal spring integrator is the reusable `leapfrog(force, x0, v0, gamma, dt, steps, mass)` function in `damped_leapfrog.py`.

[Full self-contained report](report.html) · [Protocol](data/protocol.json) · [Water measurements](data/results.json) · [Hoop measurements](data/hoop.json) · [Checksums](manifest_sha256.json)

## Limits and next steps

{limits}

## Sources

The spring and rotating-hoop equations are derived from the user-supplied textbook excerpts; the original scanned pages are not redistributed. The computational constraint implementation follows [OpenMM CustomIntegrator documentation](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html). The distinction between friction-only dynamics and thermal Langevin dynamics is documented in the [OpenMM integrator guide](https://docs.openmm.org/latest/userguide/theory/04_integrators.html).

[Original experiment](../README.md) · [Undamped finite-time study](../response-extension/README.md)
'''
    (HERE/'README.md').write_text(readme, encoding='utf-8')
    def image(name):
        return 'data:image/png;base64,'+base64.b64encode((HERE/name).read_bytes()).decode()
    html=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Damping, leapfrog, and slow motion</title>
<style>body{{font:17px/1.65 system-ui,sans-serif;max-width:1080px;padding:0 24px;margin:45px auto;color:#16313b;background:#fafcfd}}h1{{font-size:40px;line-height:1.15}}h2{{margin-top:40px;font-size:25px}}.lead{{font-size:21px}}.box{{padding:18px;background:#e8f2f3;border-left:4px solid #087e83}}img{{width:100%;margin:18px 0}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #ccd9df}}th{{background:#e7eef1}}pre{{padding:18px;background:#eef3f5;overflow:auto}}a{{color:#087e83}}</style>
<p>REPRODUCIBLE DAMPING EXTENSION · 9 OCTOBER 2026</p><h1>Damping changes the motion.<br>The forces set the equilibria.</h1>
<p class="lead">Velocity damping is now implemented and tested. A spring demonstrates when oscillations stop; the rotating hoop shows when a first-order approximation becomes accurate; a matched water experiment measures the change in molecular response.</p>
<div class="box">Drag force: F<sub>drag</sub>=−mγv. Friction removes energy. Without a compensating random thermal force, it does not maintain the original temperature.</div>
<img src="{image('damping_results.png')}" alt="Damped spring curves and measured water response and energy changes">
<h2>Water: measured effects of two drag rates</h2><p>{methods}</p>{html_table}<p>IDs are arbitrary, zero based. The original force leader was molecule 116; the undamped 200-fs response leader was 173. Rankings cover all 216 molecules within one state.</p>
<h2>The actual integration steps</h2><p>{scheme}</p><p>{water_scheme}</p><p>{correction}</p>
<h2>Following the rotating-hoop pages</h2><p>{scaling}</p><p>{equilibrium}</p>
<img src="{image('hoop_results.png')}" alt="Pitchfork equilibria, inertial and overdamped motion, and phase portraits showing relaxation toward the velocity nullcline">
<p>At q=2, the maximum full/overdamped angle discrepancy over slow time 0–12 falls from <strong>{hoop_pair[0]['max_angle_error_vs_overdamped']:.6f} rad at β=2</strong> to <strong>{hoop_pair[1]['max_angle_error_vs_overdamped']:.6f} rad at β=20</strong>. Both begin at φ=0.2 with zero physical angular velocity; the reduced first-order equation has no independently specified velocity, so an initial relaxation layer is expected. The half-step check changes the β=20 full solution by at most {hoop['half_step_angle_difference_q2_beta20']:.2e} rad.</p>
<p>The lower plots correspond to the latest supplied figure. Initial points are marked with dots. At smaller ε the initial motion is almost vertical, followed by slower motion close to C. The top equilibrium φ=π, and its periodic copies, are unstable and omitted from the branch diagram’s narrow vertical range.</p>
<h2>Checks and practical limits</h2><p>Analytic spring and free-drag tests passed for both implementations; the position error decreases by approximately fourfold when the time step is halved.</p><ul>{''.join('<li>'+s+'</li>' for s in checks)}</ul><p>{limits}</p>
<h2>Reproduce and inspect</h2><p>Install the parent experiment’s pinned dependencies, then run from its directory:</p><pre>{command}</pre><p>Use Reference in place of OpenCL for portable, slower water calculations. The scripts, saved norm matrices, input states, numerical checks and SHA-256 manifest accompany this report. The original experimental data are preserved.</p>
<h2>Sources</h2><p>Equations follow the supplied textbook excerpts. New figures are generated from the recorded calculations. Implementation references: <a href="https://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html">OpenMM CustomIntegrator</a> and <a href="https://docs.openmm.org/latest/userguide/theory/04_integrators.html">OpenMM integrator guide</a>.</p></html>'''
    (HERE/'report.html').write_text(html, encoding='utf-8')
    (DATA/'audit.json').write_text(json.dumps({'array_audits': audit, 'passed': True}, indent=2), encoding='utf-8')
    manifest={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest_sha256.json'}
    (HERE/'manifest_sha256.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'files':len(manifest)+1, 'water_200fs':rows, 'audits':'passed'}, indent=2))


if __name__ == '__main__':
    main()
