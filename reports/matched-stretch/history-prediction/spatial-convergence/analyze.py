"""Frozen operational spatial screens, paired effects and whole-run uncertainty."""
from pathlib import Path
import sys
import json
import csv
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_refinement import sha, save_json, source_matches, wrap, cloud_target, Flow
from features import future_target, index_at
import importlib.util
spec = importlib.util.spec_from_file_location('mirror_effects', HERE.parent/'mirror-direction'/'analyze.py')
effects = importlib.util.module_from_spec(spec); spec.loader.exec_module(effects)
compare, bootstrap, factorial = effects.prior.compare, effects.prior.bootstrap, effects.factorial
CELLS = ['original_normal', 'original_reversed', 'mirror_normal', 'mirror_reversed']
FAMILIES = ['global', 'cloud', 'pair']


def difference(a, b, cutoff=None):
    out = compare(a, b, cutoff)
    out['relative_RMSE'] = out['same_label_RMSE']/max(float(np.sqrt(np.mean(b*b))), 1e-12)
    return out


def targets(z, p):
    out = {}
    for horizon in p['horizons']:
        i = index_at(z['normal_time'], horizon)
        cells = {k: {} for k in FAMILIES}
        for label in ['normal', 'reversed']:
            center = z[label+'_positions'][i]; idx = z['neighbor_indices']
            for arrangement, prefix in [('original', ''), ('mirror', 'global_')]:
                key = arrangement+'_'+label
                cells['global'][key] = future_target(z[label+'_'+prefix+'tangent'][0], z[label+'_'+prefix+'tangent'][i], horizon)
                r0 = z['initial_offsets'] if arrangement == 'original' else z['initial_local_mirror_offsets']
                r1 = wrap(center[idx]-center[:, None]) if arrangement == 'original' else wrap(z[label+'_local_neighbor_positions'][i]-center[:, None])
                cells['cloud'][key] = cloud_target(r0, r1, horizon)
                cells['pair'][key] = np.log(np.maximum(np.linalg.norm(r1, axis=-1), 1e-12)/np.maximum(np.linalg.norm(r0, axis=-1), 1e-12)).mean(1)/horizon
        out[str(horizon)] = cells
    return out


def effect_values(cells):
    out = {family: factorial(cells[family]) for family in ['cloud', 'pair']}
    out['global_direction'] = {}
    for arrangement in ['original', 'mirror']:
        a = cells['global'][arrangement+'_normal']; b = cells['global'][arrangement+'_reversed']
        out['global_direction'][arrangement] = {'RMS': float(np.sqrt(np.mean((a-b)**2)))}
    return out


def passed_screen(value, threshold, lower=False):
    return value is not None and (value >= threshold if lower else value <= threshold)


def main():
    p = json.loads((HERE/'protocol.json').read_text()); execution = json.loads((HERE/'execution.json').read_text())
    if not execution['complete'] or len(execution['jobs']) != 18: raise ValueError('Initial matrix incomplete')
    jobs = execution['jobs']; grids = p['new_grids'].copy()
    extension = HERE/'extension-execution.json'
    if extension.exists():
        e = json.loads(extension.read_text())
        if not e['complete'] or len(e['jobs']) != 5: raise ValueError('Extension incomplete')
        jobs = jobs+e['jobs']; grids.append(p['extension_grid'])
    data = {}; meta = {}; spectra = {}; ts = {}
    for job in jobs:
        name = job['name']; path = HERE/'recorded-data'/(name+'.npz')
        m = json.loads(path.with_suffix('.json').read_text()); fp = path.with_name(name+'_endpoint-spectrum.npz')
        if sha(path) != m['sha256'] or sha(fp) != m['spectrum_sha256']: raise ValueError('Output hash mismatch')
        if sha(HERE.parent/m['checkpoint_path']) != m['checkpoint_sha256']: raise ValueError('Reference checkpoint changed')
        for rel, digest in m['sources'].items():
            if not source_matches(HERE.parent.parent/rel, digest): raise ValueError('Executed source changed '+rel)
        with np.load(path, allow_pickle=False) as z: data[name] = dict(z)
        with np.load(fp, allow_pickle=False) as z: spectra[name] = dict(z)
        if not all(np.isfinite(a).all() for a in data[name].values()): raise ValueError('Nonfinite data')
        ts[name] = targets(data[name], p); meta[name] = m
    result = {'protocol': p, 'grids': grids, 'analysis_sha256': sha(HERE/'analyze.py'),
              'runs': {}, 'spatial_comparisons': {}, 'half_step_comparisons': {},
              'coarse_time_comparisons': {}, 'secondary_horizon_spatial_comparisons': {},
              'summary': {}, 'screens': [], 'all_screens_pass': None}
    table = []; threshold = p['event_threshold']; primary = str(p['primary_horizon'])
    for name, value in ts.items():
        result['runs'][name] = {'seed': meta[name]['seed'], 'n': meta[name]['n'], 'dt': meta[name]['dt'],
            'data_sha256': meta[name]['sha256'], 'diagnostics': meta[name]['diagnostics'], 'effects': {}, 'global_comparisons': {}}
        for horizon, cells in value.items():
            result['runs'][name]['effects'][horizon] = effect_values(cells)
            comparisons = {}
            for label, a, b in [('direction_original', 'original_normal', 'original_reversed'),
                                ('direction_mirror', 'mirror_normal', 'mirror_reversed'),
                                ('mirror_normal', 'original_normal', 'mirror_normal'),
                                ('mirror_reversed', 'original_reversed', 'mirror_reversed')]:
                comparisons[label] = compare(cells['global'][a], cells['global'][b], threshold if horizon == primary else None)
            result['runs'][name]['global_comparisons'][horizon] = comparisons
    def screen(kind, identifier, value, limit, lower=False):
        strict = kind == 'difference_contracts'
        result['screens'].append({'kind': kind, 'identifier': identifier, 'value': value, 'limit': limit,
                                  'comparison': '>=' if lower else ('<' if strict else '<='),
                                  'pass': value < limit if strict else passed_screen(value, limit, lower)})
    s = p['screens']; finest = grids[-1]; previous = grids[-2]
    for horizon in p['horizons']:
        if str(horizon) == primary: continue
        h = str(horizon); rec = {}
        for seed in p['seeds']:
            rec[str(seed)] = {}
            for na, nb in zip(grids[:-1], grids[1:]):
                a = f's{seed}_n{na}_dt{p["dt"]}'; b = f's{seed}_n{nb}_dt{p["dt"]}'
                rec[str(seed)][f'{na}-{nb}'] = {family: {cell: difference(ts[a][h][family][cell], ts[b][h][family][cell]) for cell in CELLS} for family in FAMILIES}
        result['secondary_horizon_spatial_comparisons'][h] = rec
    for seed in p['seeds']:
        result['spatial_comparisons'][str(seed)] = {}; differences = {}; vectors = {}
        for na, nb in zip(grids[:-1], grids[1:]):
            a = f's{seed}_n{na}_dt{p["dt"]}'; b = f's{seed}_n{nb}_dt{p["dt"]}'
            key = f'{na}-{nb}'; record = {}
            for family in FAMILIES:
                record[family] = {}
                for cell in CELLS:
                    d = difference(ts[a][primary][family][cell], ts[b][primary][family][cell], threshold if family == 'global' else None)
                    record[family][cell] = d
                    identifier = f's{seed}/{family}/{cell}'
                    table.append({'seed': seed, 'pair': key, 'family': family, 'cell': cell,
                                  'RMSE': d['same_label_RMSE'], 'relative_RMSE': d['relative_RMSE'],
                                  'Jaccard': d.get('event_Jaccard')})
                    if identifier in differences:
                        contraction = d['same_label_RMSE']/max(differences[identifier], 1e-12)
                        screen('difference_contracts', key+'/'+identifier, contraction, 1.)
                        delta = ts[a][primary][family][cell]-ts[b][primary][family][cell]
                        prior_delta = vectors[identifier]
                        d['successive_difference_alignment'] = float(np.dot(delta, prior_delta)/max(np.linalg.norm(delta)*np.linalg.norm(prior_delta), 1e-24))
                    differences[identifier] = d['same_label_RMSE']
                    vectors[identifier] = ts[a][primary][family][cell]-ts[b][primary][family][cell]
                    if nb == finest:
                        screen('relative_target_RMSE', identifier, d['relative_RMSE'], s['finest_pair_relative_RMSE_max'])
                        if family == 'global': screen('event_Jaccard', identifier, d['event_Jaccard'], s['pointwise_event_Jaccard_min'], True)
            f = Flow(56, p['nu'], workers=1)
            record['shared_band_Eulerian_relative_RMS'] = {}
            for label in ['normal', 'reversed']:
                delta = spectra[a][label]-spectra[b][label]
                record['shared_band_Eulerian_relative_RMS'][label] = float(np.sqrt(f.inner(delta, delta)/f.inner(spectra[b][label], spectra[b][label])))
            result['spatial_comparisons'][str(seed)][key] = record
        a = f's{seed}_n{previous}_dt{p["dt"]}'; b = f's{seed}_n{finest}_dt{p["dt"]}'
        ea = result['runs'][a]['effects'][primary]; eb = result['runs'][b]['effects'][primary]
        for family in eb:
            for key in eb[family]:
                drift = abs(ea[family][key]['RMS']-eb[family][key]['RMS'])/max(eb[family][key]['RMS'], 1e-12)
                screen('effect_relative_drift', f's{seed}/{family}/{key}', drift, s['effect_RMS_relative_drift_max'])
        # Old grid56 dt.005 starts from precisely the same checkpoint/positions.
        name = f's{seed}_n56_dt0.005'; mirrorpath = HERE.parent/'mirror-direction'/'recorded-data'/(name+'.npz')
        oldpath = HERE.parent/'velocity-reversal'/'recorded-data'/(name+'.npz')
        mm = json.loads(mirrorpath.with_suffix('.json').read_text())
        if sha(mirrorpath) != mm['sha256'] or sha(oldpath) != mm['reused_archive_sha256']: raise ValueError('Coarse archived data changed')
        with np.load(mirrorpath, allow_pickle=False) as z: old = dict(z)
        with np.load(oldpath, allow_pickle=False) as z:
            for label in ['normal', 'reversed']:
                old[label+'_positions'] = z[label+'_positions']; old[label+'_tangent'] = z[label+'_tangent']
        ref = f's{seed}_n56_dt{p["dt"]}'
        if not np.array_equal(old['normal_positions'][0], data[ref]['normal_positions'][0]): raise ValueError('Changed fixed centers')
        old_targets = targets(old, p)
        result['coarse_time_comparisons'][str(seed)] = {fam: {cell: difference(old_targets[primary][fam][cell], ts[ref][primary][fam][cell], threshold if fam == 'global' else None) for cell in CELLS} for fam in FAMILIES}
    for n in p['time_control_grids']+([p['extension_grid']] if extension.exists() else []):
        a = f's{p["time_control_seed"]}_n{n}_dt{p["dt"]}'; b = f's{p["time_control_seed"]}_n{n}_dt{p["half_step_dt"]}'
        rec = {fam: {cell: difference(ts[a][primary][fam][cell], ts[b][primary][fam][cell], threshold if fam == 'global' else None) for cell in CELLS} for fam in FAMILIES}
        result['half_step_comparisons'][str(n)] = rec
        if n == finest:
            spatial = result['spatial_comparisons'][str(p['time_control_seed'])][f'{previous}-{finest}']
            for fam in FAMILIES:
                for cell in CELLS:
                    ratio = rec[fam][cell]['same_label_RMSE']/max(spatial[fam][cell]['same_label_RMSE'], 1e-12)
                    screen('temporal_fraction', fam+'/'+cell, ratio, s['time_error_fraction_of_finest_spatial_error_max'])
    for n in grids:
        names = [f's{seed}_n{n}_dt{p["dt"]}' for seed in p['seeds']]
        ef = [result['runs'][name]['effects'][primary] for name in names]
        result['summary'][str(n)] = {'effects': {fam: {key: bootstrap([e[fam][key]['RMS'] for e in ef]) for key in ef[0][fam]} for fam in ef[0]},
            'global_comparisons': {key: {metric: bootstrap([result['runs'][name]['global_comparisons'][primary][key][metric] for name in names]) for metric in ['same_label_RMSE', 'event_disagreement_fraction', 'event_Jaccard', 'correlation']} for key in result['runs'][names[0]]['global_comparisons'][primary]}}
    result['all_screens_pass'] = all(x['pass'] for x in result['screens'])
    result['failures'] = [x for x in result['screens'] if not x['pass']]
    save_json(HERE/'results.json', result)
    np.savez_compressed(HERE/'targets.npz', **{name+'_'+h+'_'+fam+'_'+cell: value for name, hs in ts.items() for h, fams in hs.items() for fam, cells in fams.items() for cell, value in cells.items()})
    with (HERE/'grid-comparisons.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0])); writer.writeheader(); writer.writerows(table)
    if not extension.exists():
        save_json(HERE/'extension-decision.json', {'extend': not result['all_screens_pass'], 'frozen_protocol_sha256': sha(HERE/'protocol.json'), 'screen_failures': result['failures']})
    plot(result)
    print(json.dumps({'grids': grids, 'screens': len(result['screens']), 'passed': sum(x['pass'] for x in result['screens']), 'failures': result['failures'], 'summary': result['summary']}, indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    grids = result['grids']; fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout='constrained')
    colors = ['#2166ac', '#d6604d', '#7b3294']
    pairs = list(zip(grids[:-1], grids[1:]))
    for family, color in zip(FAMILIES, colors):
        means = []
        for a, b in pairs:
            values = [result['spatial_comparisons'][str(seed)][f'{a}-{b}'][family][cell]['relative_RMSE'] for seed in result['protocol']['seeds'] for cell in CELLS]
            means.append(100*np.mean(values))
        axes[0, 0].plot([b for a, b in pairs], means, 'o-', label=family, color=color)
    axes[0, 0].axhline(1, color='gray', ls='--'); axes[0, 0].set(xlabel='Finer grid size', ylabel='Mean adjacent relative RMS (%)', title='Target differences; screens use every run/cell'); axes[0, 0].legend()
    for a, b in pairs:
        vals = [result['spatial_comparisons'][str(seed)][f'{a}-{b}']['global'][cell]['event_Jaccard'] for seed in result['protocol']['seeds'] for cell in CELLS]
        axes[0, 1].scatter([b]*len(vals), vals, color='#2166ac', alpha=.6)
    axes[0, 1].axhline(.95, color='gray', ls='--'); axes[0, 1].set(xlabel='Finer grid size', ylabel='High-event intersection / union', title='Fixed cutoff; all seeds and four cells')
    keys = ['arrangement_average', 'direction_average', 'interaction_difference_in_differences']
    for key, color in zip(keys, colors):
        vals = [result['summary'][str(n)]['effects']['cloud'][key]['mean'] for n in grids]
        axes[1, 0].plot(grids, vals, 'o-', label=key.replace('_', ' '), color=color)
    axes[1, 0].set(xlabel='Grid size', ylabel='Mean per-run contrast RMS', title='Anchored finite-cloud effects; different scales'); axes[1, 0].legend(fontsize=8)
    for arrangement, color in [('original', colors[0]), ('mirror', colors[1])]:
        vals = [100*result['summary'][str(n)]['global_comparisons']['direction_'+arrangement]['event_disagreement_fraction']['mean'] for n in grids]
        axes[1, 1].plot(grids, vals, 'o-', label=arrangement+' cloud', color=color)
    axes[1, 1].set(xlabel='Grid size', ylabel='Mean event disagreement (%)', title='Effect of reversing physical motion'); axes[1, 1].legend()
    fig.savefig(HERE/'convergence.png', dpi=160); fig.savefig(HERE/'convergence.svg'); plt.close(fig)


if __name__ == '__main__': main()
