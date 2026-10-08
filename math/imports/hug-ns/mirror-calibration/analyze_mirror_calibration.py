"""Verify complete signed-coordinate and full-field reflection controls."""
from pathlib import Path
import argparse
import json
import os
import time
import numpy as np
from mirror_calibration import OUT, HERE, build, inner, m, mirror, save


def read(name):
    path = OUT/name/'result.json'
    return json.loads(path.read_text()) if path.exists() else None


def compare_curves(a, b):
    count = min(len(a['series']), len(b['series']))
    out = {'through': a['series'][count-1]['t'], 'fields': {}}
    for name in ('s_fraction_of_initial_norm', 'E_mirror_mismatch'):
        va = np.array([r[name] for r in a['series'][:count]])
        vb = np.array([r[name] for r in b['series'][:count]])
        out['fields'][name] = {'relative_max_difference': float(np.max(abs(va-vb))/max(np.max(abs(vb)), 1e-30)),
                              'relative_curve_l2': float(np.linalg.norm(va-vb)/max(np.linalg.norm(vb), 1e-30))}
    return out


def analyze():
    state = json.loads((OUT/'task-state.json').read_text())
    pairs = []; zero = []; time_checks = []; grid_checks = []; mode_checks = []
    for n in (33, 49):
        for half in (False, True):
            label = 'half' if half else 'base'
            control = read(f'n{n}-{label}-a0.0')
            if control:
                zero.append({'n': n, 'half': half, 'through': control['series'][-1]['t'],
                             'max_E': max(r['E_mirror_mismatch'] for r in control['series']),
                             'max_abs_signed': max(abs(r['s_fraction_of_initial_norm']) for r in control['series'])})
            for amp in (.001, .005):
                a = read(f'n{n}-{label}-a{amp}'); b = read(f'n{n}-{label}-a{-amp}')
                if not a or not b: continue
                count = min(len(a['series']), len(b['series']))
                full_errors = []
                for i in range(count):
                    positive = np.load(OUT/f'n{n}-{label}-a{amp}'/f'field-{i:02d}.npy')
                    negative = np.load(OUT/f'n{n}-{label}-a{-amp}'/f'field-{i:02d}.npy')
                    diff = negative-mirror(positive)
                    full_errors.append(float(np.sqrt(np.sum(diff*diff)/np.sum(positive*positive))))
                pairs.append({'n': n, 'half': half, 'amplitude': amp, 'through': a['series'][count-1]['t'],
                              'max_abs_signed_pair_sum': max(abs(p['s_fraction_of_initial_norm']+q['s_fraction_of_initial_norm']) for p, q in zip(a['series'], b['series'])),
                              'max_abs_E_pair_difference': max(abs(p['E_mirror_mismatch']-q['E_mirror_mismatch']) for p, q in zip(a['series'], b['series'])),
                              'max_relative_full_field_reflection_error': max(full_errors)})
                last = a['series'][-1]
                mode_checks.append({'n': n, 'half': half, 'amplitude': amp, 'through': last['t'],
                                    'residual_asymmetry_norm_fraction': last['residual_asymmetry_fraction'],
                                    'relative_E_error_if_only_one_mode_used': abs(last['E_single_mode_prediction']-last['E_mirror_mismatch'])/last['E_mirror_mismatch'],
                                    'exact_decomposition_error': last['exact_decomposition_error']})
        for amp in (.001, -.001, .005, -.005):
            base = read(f'n{n}-base-a{amp}'); half = read(f'n{n}-half-a{amp}')
            if base and half: time_checks.append({'n': n, 'amplitude': amp, **compare_curves(base, half)})
    for label in ('base', 'half'):
        for amp in (.001, -.001, .005, -.005):
            a, b = read(f'n33-{label}-a{amp}'), read(f'n49-{label}-a{amp}')
            if a and b: grid_checks.append({'timestep': label, 'amplitude': amp, **compare_curves(a, b)})
    scaling = []
    for n in (33, 49):
        for label in ('base', 'half'):
            a, b = read(f'n{n}-{label}-a0.001'), read(f'n{n}-{label}-a0.005')
            if a and b:
                count = min(len(a['series']), len(b['series']))
                scaling.append({'n': n, 'timestep': label, 'through': a['series'][count-1]['t'],
                    'max_fivefold_signed_amplitude_scaling_error': max(abs(q['s_fraction_of_initial_norm']/(5*p['s_fraction_of_initial_norm'])-1)
                         for p, q in zip(a['series'], b['series'])),
                    'max_fivefold_E_scaling_error': max(abs(q['E_mirror_mismatch']/(5*p['E_mirror_mismatch'])-1)
                         for p, q in zip(a['series'], b['series']))})
    report = {'status': state['status'], 'completed': len(state['completed']), 'planned': state['planned'],
              'zero_bias_controls': zero, 'opposite_pairs': pairs, 'time_step_comparisons': time_checks,
              'grid_comparisons': grid_checks, 'single_mode_capture': mode_checks, 'small_amplitude_scaling': scaling,
              'normalization': 'a=<u-Mu,phi>/(2<phi,phi>), s_calibrated=a/norm(u_sym(0)); E retains all antisymmetric modes.',
              'scope': 'Complete fluid calibration runs with an explicit Gaussian-derived template. These are not an exact reproduction of uploaded trajectories whose template and time loop were omitted.',
              'reduced_law': 'No cubic coefficients fitted or validated by these symmetry checks.',
              'updated_epoch': time.time()}
    if state['status'] == 'complete':
        assert len(pairs) == 8 and len(time_checks) == 8 and len(grid_checks) == 8
        assert max(z['max_E'] for z in zero) < 1e-10
        assert max(p['max_relative_full_field_reflection_error'] for p in pairs) < 1e-10
        report['reflection_checks'] = 'passed'
        report['time_step_screen'] = 'passed (<0.1%)' if max(v['relative_max_difference'] for r in time_checks for v in r['fields'].values()) < .001 else 'needs investigation'
        report['two_grid_agreement_screen'] = 'passed (<1%)' if max(v['relative_max_difference'] for r in grid_checks for v in r['fields'].values()) < .01 else 'needs investigation'
        report['grid_scope'] = 'Agreement between two grids for these observables; does not establish full-field convergence or convergence order.'
    save(OUT/'analysis.json', report)
    print(json.dumps({'status': report['status'], 'completed': report['completed'],
                      'max_pair_field_error': max((p['max_relative_full_field_reflection_error'] for p in pairs), default=None),
                      'timestep_screen': report.get('time_step_screen'), 'grid_screen': report.get('two_grid_agreement_screen')}), flush=True)
    return report


def rhs_checks():
    results = []
    for n in (33, 49):
        base, phi, norm, ops = build(n)
        u = base+.005*norm*phi
        r = np.array(m.rhs(u, *ops[:-1]))
        mirrored_rhs = np.array(m.rhs(mirror(u), *ops[:-1]))
        error = np.linalg.norm(mirrored_rhs-mirror(r))/np.linalg.norm(r)
        assert error < 1e-10
        results.append({'n': n, 'rhs_reflection_relative_error': float(error)})
    save(OUT/'rhs-reflection-check.json', {'status': 'passed', 'checks': results})


def plots(report):
    os.environ.setdefault('MPLCONFIGDIR', str(HERE/'mpl-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for n, color in ((33, '#267f88'), (49, '#ba6633')):
        for sign, style in ((1, '-'), (-1, '--')):
            run = read(f'n{n}-half-a{sign*.001}') or read(f'n{n}-base-a{sign*.001}')
            if not run: continue
            rows = run['series']; t = [r['t'] for r in rows]
            label = f'{n}³, {sign*.1:+.1f}%'
            axes[0, 0].plot(t, [r['s_fraction_of_initial_norm'] for r in rows], style, color=color, label=label)
            axes[0, 1].plot(t, [r['E_mirror_mismatch'] for r in rows], style, color=color, label=label)
            if sign == 1:
                axes[1, 0].plot(t, [r['E_mirror_mismatch'] for r in rows], color=color, label=f'{n}³ all asymmetry')
                axes[1, 0].plot(t, [r['E_single_mode_prediction'] for r in rows], ':', color=color, label=f'{n}³ chosen pattern only')
                axes[1, 1].plot(t, [r['residual_asymmetry_fraction'] for r in rows], color=color, label=f'{n}³')
    for ax, title, ylabel in zip(axes.flat, ('Opposite signed responses', 'Equal mirror-error magnitudes', 'What one signed coordinate captures', 'Asymmetry outside the chosen pattern'),
                                ('Signed amplitude / initial field norm', 'Mirror mismatch E', 'Mirror mismatch E', 'Fraction of asymmetry norm')):
        ax.set(title=title, xlabel='Model time', ylabel=ylabel); ax.grid(alpha=.18); ax.legend(fontsize=8)
    fig.suptitle('Completed mirror calibration · fixed viscosity 0.01 · no external force', fontsize=15)
    output = HERE/'outputs'
    output.mkdir(exist_ok=True)
    fig.savefig(output/'mirror-calibration.png', dpi=150); plt.close(fig)
    content = '''<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Mirror calibration results</title><style>body{font:17px/1.55 system-ui;max-width:1100px;margin:30px auto;padding:20px;color:#173a40;background:#f5f5f1}img{width:100%}a{color:#176875}section{background:white;padding:20px;margin:20px 0;border-radius:12px}</style><h1>Mirror calibration</h1>'''
    content += f'<section><b>{report["completed"]} / {report["planned"]} runs complete.</b><p>Zero bias, ±0.1% and ±0.5%; 33³ and 49³; base and half timesteps; model time 0.20.</p><p>Reflection checks: {report.get("reflection_checks", "in progress")}. Timestep agreement: {report.get("time_step_screen", "in progress")}. Two-grid observable agreement: {report.get("two_grid_agreement_screen", "in progress")}.</p></section>'
    content += '<img src="mirror-calibration.png" alt="Signed amplitudes, unsigned mirror mismatch, and asymmetry outside the selected pattern"><section><h2>The connection</h2><p>The signed coordinate measures the component along the chosen antisymmetric pattern. E measures all mirror mismatch. The starting signed coordinate equals the requested perturbation amplitude after correcting the factor of two.</p><p>During evolution, some asymmetry leaves the chosen pattern. The residual is measured explicitly, so one coordinate is not assumed to capture everything.</p><p>These runs verify the measurement and reflection relationships. They do not establish the proposed cubic evolution equation or resolve peak-vorticity behavior.</p></section><section><a href="../mirror-calibration/analysis.json">Detailed checks</a> · <a href="../mirror-calibration/protocol.json">Exact construction and definitions</a></section>'
    (output/'mirror-calibration.html').write_text(content, encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--wait', action='store_true'); args = parser.parse_args()
    rhs_checks()
    while True:
        report = analyze()
        if report['status'] == 'complete' or not args.wait: break
        time.sleep(20)
    plots(report)
