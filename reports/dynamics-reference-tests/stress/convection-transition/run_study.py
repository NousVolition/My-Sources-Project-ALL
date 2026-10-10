"""Execute the declared bounded onset study and save every comparison."""
import csv
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.stats import t as student_t
from convection import ROOT, PLAN, RAC, Q, modal_matrix, leading_rate, run, COLUMNS

DATA = ROOT/'data'


def dump(name, value):
    (DATA/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def main():
    started = time.perf_counter()
    DATA.mkdir(exist_ok=True)
    results = []
    arrays = {}
    def case(label, **kwargs):
        result, rows, maps, map_times = run(**kwargs)
        result['label'] = label
        results.append(result)
        np.savez_compressed(DATA/(label+'.npz'), records=rows, fields=maps, field_times=map_times)
        arrays[label] = (rows, maps, map_times)
        print(label, 'rate=', result['fitted_rate'], flush=True)
        return result

    for ratio in PLAN['rayleigh_over_critical']:
        for seed in PLAN['seeds']:
            case(f'primary_r{ratio:g}_s{seed}', ratio=ratio, seed=seed)
        for n in [16, 32]:
            case(f'grid{n}_r{ratio:g}', ratio=ratio, n=n)
        case(f'halfstep_r{ratio:g}', ratio=ratio, dt=.005)
    for ratio in PLAN['quarter_step_ratios']:
        case(f'quarterstep_r{ratio:g}', ratio=ratio, dt=.0025)
        case(f'halfamplitude_r{ratio:g}', ratio=ratio, amplitude=5e-9)
    for n in [16, 24, 32]:
        case(f'nonlinear_grid{n}', ratio=1.2, n=n, amplitude=.02, dt=.0025, duration=1.)
    case('zero_disturbance', ratio=1.2, amplitude=0.)
    case('buoyancy_off', ratio=1.2, buoyancy=False)
    case('coarse_negative_control', ratio=1.05, dt=.2)
    case('repeat_selected_grid32', ratio=1.05, n=32)

    # Added after the unchanged energy tolerance rejected dt=.01 near onset.
    for ratio in PLAN['rayleigh_over_critical']:
        for seed in PLAN['seeds']:
            case(f'refined_primary_r{ratio:g}_s{seed}', ratio=ratio, seed=seed, dt=.005)
        for n in [16,32]:
            case(f'refined_grid{n}_r{ratio:g}', ratio=ratio, n=n, dt=.005)
    for ratio in [.95,1.05]:
        case(f'fixed_absolute_temperature_r{ratio:g}', ratio=ratio, dt=.005, thermal_initial_multiplier=1/ratio)
    case('repeat_refined_grid32', ratio=1.05, n=32, dt=.005)

    for result in results:
        tolerance=PLAN['acceptance']
        result['rate_accuracy_passed']=(result['rate_error'] <= tolerance['growth_rate_absolute_error']) if result['rate_error'] is not None else None
        result['energy_budgets_passed']=max(result['relative_kinetic_budget'],result['relative_thermal_budget'])<=tolerance['energy_budget_relative']
        result['boundary_checks_passed']=max(result['relative_divergence'],result['relative_wall'])<=tolerance['wall_relative']

    by = {r['label']: r for r in results}
    checks = []
    def check(name, observed, limit, passed):
        checks.append(dict(name=name, observed=observed, limit=limit, passed=bool(passed)))
    a = PLAN['acceptance']
    accepted = [r for r in results if r['label'] not in ('coarse_negative_control',) and not r['label'].startswith('nonlinear')]
    rates = [r['rate_error'] for r in accepted if r['rate_error'] is not None and r['buoyancy']]
    check('small-signal rates match independent linear theory', max(rates), a['growth_rate_absolute_error'], max(rates) <= a['growth_rate_absolute_error'])
    for ratio in [.8, .95, 1.05, 1.2]:
        vals = [r['fitted_rate'] for r in results if r['label'].startswith('refined_primary_') and r['ratio']==ratio]
        check(f'growth sign r={ratio:g}', vals, 'all below/above zero according to threshold', all(v*(ratio-1)>0 for v in vals))
    neutral = max(abs(r['fitted_rate']) for r in results if r['label'].startswith('refined_primary_') and r['ratio']==1.)
    check('critical case is within declared neutral band', neutral, a['neutral_rate_band'], neutral < a['neutral_rate_band'])
    for field, key in [('relative_divergence','divergence_relative'), ('relative_wall','wall_relative'),
                       ('relative_kinetic_budget','energy_budget_relative'), ('relative_thermal_budget','energy_budget_relative')]:
        value = max(r[field] for r in results if r['dt'] <= .005)
        check('refined '+field, value, a[key], value <= a[key])
    spatial = []; temporal = []; amplitude = []
    for ratio in PLAN['rayleigh_over_critical']:
        base = by[f'primary_r{ratio:g}_s410']
        spatial.append(abs(by[f'refined_primary_r{ratio:g}_s410']['fitted_rate']-by[f'refined_grid32_r{ratio:g}']['fitted_rate']))
        temporal.append(abs(base['fitted_rate']-by[f'halfstep_r{ratio:g}']['fitted_rate']))
    check('grid 24/32 small-signal rates', max(spatial), a['refined_rate_difference'], max(spatial) < a['refined_rate_difference'])
    check('dt .01/.005 small-signal rates', max(temporal), a['refined_rate_difference'], max(temporal) < a['refined_rate_difference'])
    # Independent Fourier interpolation onto a common grid uses scipy.signal.
    from scipy.signal import resample
    def field_difference(left, right):
        l = arrays[left][1][-1]; r = arrays[right][1][-1]
        target = r.shape[-1]
        l = resample(resample(l, target, axis=1), target, axis=2)
        return float(np.linalg.norm(l-r)/np.linalg.norm(r))
    grid_fields = [field_difference('nonlinear_grid16','nonlinear_grid24'), field_difference('nonlinear_grid24','nonlinear_grid32')]
    check('larger-amplitude final field grid 24/32', grid_fields[-1], a['refined_field_relative_error'], grid_fields[-1] < a['refined_field_relative_error'])
    for ratio in PLAN['quarter_step_ratios']:
        full = arrays[f'primary_r{ratio:g}_s410'][0][:, 1]
        half = 2*arrays[f'halfamplitude_r{ratio:g}'][0][:, 1]
        amplitude.append(float(np.max(abs(full-half))/np.max(abs(full))))
    check('half-amplitude linear-regime control', max(amplitude), a['amplitude_control_relative'], max(amplitude) < a['amplitude_control_relative'])
    absolute_temperature=max(abs(by[f'fixed_absolute_temperature_r{ratio:g}']['fitted_rate']-by[f'refined_primary_r{ratio:g}_s410']['fitted_rate']) for ratio in [.95,1.05])
    check('fixed absolute temperature perturbation growth rate',absolute_temperature,a['growth_rate_absolute_error'],absolute_temperature<a['growth_rate_absolute_error'])
    check('no perturbation remains exactly zero', float(np.max(abs(arrays['zero_disturbance'][0][:, 1:]))), 0., not np.any(arrays['zero_disturbance'][0][:, 1:]))
    check('removing buoyancy removes growth', by['buoyancy_off']['fitted_rate'], '<0', by['buoyancy_off']['fitted_rate'] < 0)
    coarse = by['coarse_negative_control']
    check('coarse negative control detected: erases genuine instability', coarse['fitted_rate'], 'wrong sign, rejected as inaccurate', coarse['fitted_rate'] < 0 and coarse['rate_error'] > a['growth_rate_absolute_error'])
    repeated = all(np.array_equal(x,y) for x,y in zip(arrays['refined_grid32_r1.05'],arrays['repeat_refined_grid32']))
    check('selected complete run reproduced exactly', repeated, True, repeated)
    groups = []
    for ratio in PLAN['rayleigh_over_critical']:
        r = [v for v in results if v['label'].startswith('refined_primary_') and v['ratio']==ratio]
        vals = np.array([v['fitted_rate'] for v in r])
        half = float(student_t.ppf(.975,len(vals)-1)*vals.std(ddof=1)/np.sqrt(len(vals)))
        groups.append(dict(ratio=ratio, rayleigh=ratio*RAC, analytic_rate=leading_rate(ratio*RAC),
                           mean_rate=float(vals.mean()), seed_95_t_interval=[float(vals.mean()-half),float(vals.mean()+half)],
                           amplification_mean=float(np.mean([v['mode_amplification'] for v in r]))))

    inertia = []
    for ratio in PLAN['inertia_comparison']['rayleigh_ratios']:
        times = np.linspace(0, .8 if ratio < 0 else 6, 801)
        ra = ratio*RAC
        for pr in PLAN['inertia_comparison']['prandtl']:
            # Matched theta and well-prepared Stokes velocity at time zero.
            initial = np.array([Q*ratio, 1.])
            matrix = modal_matrix(ra,pr)
            full = np.array([expm(matrix*t)@initial for t in times])
            reduced_theta = np.exp(Q*(ratio-1)*times)
            reduced = np.column_stack((Q*ratio*reduced_theta,reduced_theta))
            name=f'inertia_r{ratio:g}_pr{pr:g}'
            np.savez_compressed(DATA/(name+'.npz'), t=times, full=full, stokes=reduced)
            eig=np.linalg.eigvals(matrix)
            inertia.append(dict(label=name, ratio=ratio, prandtl=pr,
                full_leading_rate=float(eig.real.max()), full_frequency=float(abs(eig.imag).max()),
                stokes_rate=Q*(ratio-1), theta_zero_crossings=int(np.sum(full[1:,1]*full[:-1,1]<0)),
                theta_relative_curve_error=float(np.linalg.norm(full[:,1]-reduced_theta)/np.linalg.norm(full[:,1]))))
    check('stable fluid mode oscillates with inertia at Pr=1', next(r['theta_zero_crossings'] for r in inertia if r['ratio']==-4 and r['prandtl']==1), '>0', next(r['theta_zero_crossings'] for r in inertia if r['ratio']==-4 and r['prandtl']==1)>0)
    dump('inertia.json', inertia)
    with (DATA/'cases.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
    dump('record_columns.json', COLUMNS)
    dump('checks.json', checks)
    dump('summary.json', dict(question=PLAN['question'], evidence=PLAN['evidence'], critical_rayleigh=RAC,
         critical_horizontal_wavenumber=float(np.pi/np.sqrt(2)), pde_runs=len(results), modal_comparisons=len(inertia),
         primary_groups=groups, nonlinear_field_grid_differences=grid_fields,
         max_temporal_rate_change=max(temporal), max_spatial_rate_change=max(spatial),
         coarse_control=coarse, checks_passed=sum(c['passed'] for c in checks), checks_total=len(checks),
         full_repeat_scope='Selected grid32, ratio1.05 full runs at dt0.01 and dt0.005; refined repeat compared here. Not the whole matrix.',
         repeat_identical=repeated, elapsed_seconds=time.perf_counter()-started,
         environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
         source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['plan.json','convection.py','run_study.py']}))
    print(json.dumps(dict(pde_runs=len(results),checks_passed=sum(c['passed'] for c in checks),checks_total=len(checks),elapsed_seconds=time.perf_counter()-started)),flush=True)
    if not all(c['passed'] for c in checks):
        raise SystemExit('Some declared checks failed; retain and report them before interpreting the study.')

if __name__=='__main__':
    main()
