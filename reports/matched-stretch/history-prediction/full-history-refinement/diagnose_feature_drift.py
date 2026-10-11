"""Descriptive additive-score diagnosis; does not refit or change any screen."""
from pathlib import Path
import argparse
import json
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    cli=argparse.ArgumentParser();cli.add_argument('--stage',type=int);args=cli.parse_args()
    suffix=f'-{args.stage}' if args.stage else ''
    result=json.loads((HERE/f'results{suffix}.json').read_text())
    models=json.loads((HERE.parent/'results'/'primary_models.json').read_text())
    import hashlib
    assert hashlib.sha256((HERE.parent/'results'/'primary_models.json').read_bytes()).hexdigest()==result['frozen_models_sha256']
    lo,hi=result['grids'][-2:]
    out={'grids':[lo,hi],'models_refitted':False,
         'interpretation':'Exact decomposition of a fitted linear score change, not causal attribution. Group RMS values do not add; covariances matter. Event decomposition is of log odds, not probabilities.',
         'models':{}}
    with np.load(HERE/f'predictions-and-features{suffix}.npz',allow_pickle=False) as z:
        for name in ['A_full','B_history']:
            fields=['present']+(['history'] if name=='B_history' else [])
            a=np.column_stack([z[f'n{lo}_{field}'] for field in fields])
            b=np.column_stack([z[f'n{hi}_{field}'] for field in fields])
            groups={'velocity_and_speed':[0,1,2,17],
                    'velocity_gradient_and_summaries':list(range(3,17)),
                    'present_neighbor_geometry':[i for i in range(18,31) if i!=28],
                    'present_strain_geometry_alignment':[28]}
            if name=='B_history':groups['past_descriptors']=list(range(31,54))
            assert sorted(i for indices in groups.values() for i in indices)==list(range(a.shape[1]))
            record={}
            for kind in ['continuous','event']:
                model=models['models'][name][kind]
                terms=(b-a)/np.asarray(model['scale'])*np.asarray(model['coefficient']).reshape(-1)
                parts=np.column_stack([terms[:,indices].sum(1) for indices in groups.values()])
                total=terms.sum(1)
                np.testing.assert_allclose(parts.sum(1),total,atol=1e-14,rtol=1e-12)
                if kind=='continuous':
                    saved=z[f'n{hi}_{name}_prediction']-z[f'n{lo}_{name}_prediction']
                    np.testing.assert_allclose(total,saved,atol=1e-12,rtol=1e-10)
                record[kind]={'score':'log_odds' if kind=='event' else 'continuous_prediction',
                    'total_RMS':float(np.sqrt(np.mean(total**2))),
                    'group_order':list(groups),'group_indices':groups,
                    'group_RMS':{label:float(np.sqrt(np.mean(parts[:,i]**2))) for i,label in enumerate(groups)},
                    'uncentered_group_second_moment_matrix':(parts.T@parts/len(parts)).tolist(),
                    'feature_contribution_RMS':dict(zip(models['models'][name]['feature_names'],np.sqrt(np.mean(terms**2,axis=0)).tolist())),
                    'exact_reconstruction_max_error':float(abs(parts.sum(1)-total).max())}
            out['models'][name]=record
    path=HERE/f'feature-drift{suffix}.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'grids':[lo,hi],'models':{name:{kind:{k:v for k,v in r.items() if k in ['total_RMS','group_RMS','score']} for kind,r in rec.items()} for name,rec in out['models'].items()}},indent=2))


if __name__=='__main__':main()
