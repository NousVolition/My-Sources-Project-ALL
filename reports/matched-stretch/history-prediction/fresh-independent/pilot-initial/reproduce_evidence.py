"""Read-only reproduction of the initial pilot failure; optional plot export."""
from pathlib import Path
import sys, json, argparse, hashlib
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY.parent))
from features import dataset, HISTORY_NAMES, neighbors, wrap, index_at


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot-directory', type=Path)
    args = parser.parse_args()
    protocol = json.loads((STUDY/'protocol.json').read_text())
    audit = json.loads((HERE/'pilot-audit-80.json').read_text())
    verification = json.loads((HERE/'verification.json').read_text())
    diagnosis = json.loads((HERE/'feature-diagnosis.json').read_text())
    for row in verification['initial_pilot_archives']:
        path = STUDY/'recorded-data'/(row['name']+'.npz')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256']
    for screen in audit['screens']:
        passed = (screen['value'] < screen['bound'] if screen['metric'] == 'target_RMS_ratio'
                  else screen['value'] <= screen['bound'] if screen['direction'] == 'max'
                  else screen['value'] >= screen['bound'])
        assert screen['pass'] == passed
    old = json.loads((STUDY.parent/'results/primary_models.json').read_text())
    scale = np.array(old['models']['B_history']['continuous']['scale'])[len(old['models']['A_full']['continuous']['scale']):]
    options = dict(lookback=protocol['primary_lookback'], horizon=protocol['primary_horizon'],
                   k=protocol['neighbors'], count=protocol['base']['markers'], stride=1)
    paths = [STUDY/f'recorded-data/s50000_n{n}_nu0.02_dt0.0025.npz' for n in [56,80]]
    data = [dataset([path], protocol, **options) for path in paths]
    for key in ['group','anchor','particle']:
        np.testing.assert_array_equal(data[0][key], data[1][key])
    delta = (data[0]['history']-data[1]['history'])/scale
    value = float(np.sqrt(np.mean(delta**2)))
    np.testing.assert_allclose(value, diagnosis['saved_RMS'], rtol=1e-12, atol=1e-15)
    feature_rms = np.sqrt(np.mean(delta**2,axis=0))
    for row in diagnosis['top_feature_RMS']:
        np.testing.assert_allclose(feature_rms[HISTORY_NAMES.index(row['feature'])],
                                   row['standardized_RMS'], rtol=1e-12, atol=1e-15)
    row_rms = np.sqrt(np.mean(delta**2,axis=1))
    largest = int(np.argmax(row_rms))
    assert int(data[0]['particle'][largest]) == 109
    assert np.isclose(data[0]['anchor'][largest],.6)
    fraction = float(row_rms[largest]**2/(len(row_rms)*value**2))
    np.testing.assert_allclose(fraction, diagnosis['largest_row_fraction_of_total_squared_history_difference'])
    for expected in diagnosis['current_neighbor_selection_for_largest_row']:
        path = STUDY/f"recorded-data/s50000_n{expected['grid']}_nu0.02_dt0.0025.npz"
        with np.load(path,allow_pickle=False) as z:
            x = z['positions'][index_at(z['time'],.6)]
        assert neighbors(x,8)[109].tolist() == expected['eight_neighbors']
        distance = np.linalg.norm(wrap(x-x[109]),axis=1)
        order = np.argsort(distance)
        np.testing.assert_allclose(distance[order[9]]-distance[order[8]],
                                   expected['distance_margin'], rtol=1e-10, atol=1e-14)
    print(json.dumps(dict(archives_verified=16, screens_pass=115, screens_total=116,
                         reconstructed_history_RMS=value, bound=.005,
                         largest_row_fraction=fraction, frozen_failure_reproduced=True)))
    if args.plot_directory:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        args.plot_directory.mkdir(parents=True,exist_ok=True)
        fig,axes = plt.subplots(1,2,figsize=(11.2,4.4),layout='constrained')
        locations = np.arange(4)
        for ci,label in enumerate(['56\u219280','80\u2192112']):
            for ax,metric,scale_value in [(axes[0],'target_relative_RMS',100),
                                          (axes[1],'history_original_train_standardized_RMS',1)]:
                values = [next(s['value'] for s in audit['screens']
                               if s['comparison']==label and s['seed']==seed and s['metric']==metric)*scale_value
                          for seed in audit['seeds']]
                ax.bar(locations+(ci-.5)*.34,values,width=.32,label=label,color=['#5B80A5','#C68647'][ci])
        for ax in axes:
            ax.set_xticks(locations,[str(s) for s in audit['seeds']])
            ax.set_xlabel('Independent pilot seed'); ax.legend(frameon=False)
            ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True)
        axes[0].axhline(.25,color='#9F4242',linestyle='--')
        axes[0].set_ylabel('Future target relative RMS difference (%)')
        axes[0].set_title('Future target screens pass'); axes[0].set_ylim(0,.28)
        axes[0].text(2.05,.254,'limit 0.25%',color='#9F4242',fontsize=9)
        axes[1].axhline(.005,color='#9F4242',linestyle='--')
        axes[1].set_ylabel('History RMS difference / old training scales')
        axes[1].set_title('One history stability screen fails')
        axes[1].text(2.0,.0052,'limit 0.005',color='#9F4242',fontsize=9)
        fig.savefig(args.plot_directory/'initial-accuracy.png',dpi=180)
        fig.savefig(args.plot_directory/'initial-accuracy.svg')
        plt.close(fig)


if __name__ == '__main__':
    main()
