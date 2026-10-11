"""Verify completed archives, fixed models, saved features/predictions and stage bounds."""
from pathlib import Path
import sys
import json
import argparse
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent
sys.path.insert(0,str(HERE))
from analyze_history import checked
from run_history import sha
from export_models import predict_model


def main():
    cli=argparse.ArgumentParser();cli.add_argument('--stage',type=int);args=cli.parse_args()
    suffix=f'-{args.stage}' if args.stage else ''
    p=json.loads((HERE/'protocol.json').read_text())
    r=json.loads((HERE/f'results{suffix}.json').read_text());execution=json.loads((HERE/f'execution{suffix}.json').read_text())
    if args.stage:
        decision=json.loads((HERE/f'decision{suffix}.json').read_text())
        assert decision['protocol_sha256']==sha(HERE/'protocol.json')
        assert decision['results_sha256']==sha(HERE/f'results{suffix}.json')
        done={'complete':True,'results_sha256':decision['results_sha256'],
              'report_sha256':sha(HERE/f'README{suffix}.md'),'analysis_sha256':r['analysis_sha256'],
              'frozen_models_sha256':r['frozen_models_sha256'],'grids':r['grids'],
              'new_simulations':sum(not rec['reused'] for rec in execution['records']),
              'all_screens_pass':decision['screens_pass']}
        assert r['grids'][-1]==args.stage
    else:
        done=json.loads((HERE/'completion.json').read_text())
    assert done['complete'] and execution['complete']
    assert done['results_sha256']==sha(HERE/f'results{suffix}.json')
    assert done['report_sha256']==sha(HERE/f'README{suffix}.md')
    assert done['analysis_sha256']==sha(HERE/'analyze_history.py')
    assert done['frozen_models_sha256']==p['models']['sha256']==sha(BASE/'results'/'primary_models.json')
    assert r['models_refitted'] is False and r['protocol']==p
    assert execution['grids']==r['grids']==done['grids']
    expected={(seed,n,p['base']['dt'],'main') for n in r['grids'] for seed in p['test_seeds']}
    expected|={(seed,n,p['half_step_dt'],'half') for n in p['half_step_grids'] if n in r['grids'] for seed in p['half_step_seeds']}
    actual={(rec['seed'],rec['config']['n'],rec['config']['dt'],rec['kind']) for rec in execution['records']}
    assert actual==expected and len(actual)==len(execution['records'])
    assert done['new_simulations']<=p['max_new_simulations']
    assert done['all_screens_pass']==r['all_screens_pass']==all(s['pass'] for s in r['screens'])
    records=[]
    for rec in execution['records']:
        path=HERE/rec['path'];m=checked(path,rec['config'])
        assert m['sha256']==rec['sha256'] and m['seed']==rec['seed']
        with np.load(path,allow_pickle=False) as z:
            assert z['positions'].shape==(49,512,3)
            np.testing.assert_array_equal(z['tangent'][0],np.broadcast_to(np.eye(3),(512,3,3)))
        records.append({'path':rec['path'],'sha256':rec['sha256'],'kind':rec['kind'],'reused':rec['reused']})
    model=json.loads((BASE/'results'/'primary_models.json').read_text())
    observed_rows=len(p['test_seeds'])*len(p['primary']['anchors'])*p['base']['markers'];errors={}
    with np.load(HERE/f'predictions-and-features{suffix}.npz',allow_pickle=False) as z:
        for n in r['grids']:
            target=z[f'n{n}_target'];group=z[f'n{n}_group'];anchor=z[f'n{n}_anchor'];particle=z[f'n{n}_particle']
            assert target.shape==(observed_rows,) and np.isfinite(target).all()
            assert set(group)==set(p['test_seeds']) and set(anchor)==set(p['primary']['anchors'])
            assert set(particle)==set(range(p['base']['markers']))
            for seed in p['test_seeds']:
                for t in p['primary']['anchors']:
                    ix=(group==seed)&(anchor==t)
                    assert ix.sum()==p['base']['markers'] and len(set(particle[ix]))==p['base']['markers']
            for name,fields in [('A_full',[f'n{n}_present']),('B_history',[f'n{n}_present',f'n{n}_history'])]:
                X=np.column_stack([z[field] for field in fields])
                for kind,output_suffix in [('continuous','prediction'),('event','probability')]:
                    pred=predict_model(model['models'][name][kind],X)
                    saved=z[f'n{n}_{name}_{output_suffix}'];error=float(abs(pred-saved).max())
                    assert error<1e-12;errors[f'n{n}/{name}/{kind}']=error
    out={'complete':True,'stage_only':args.stage,'bounded_matrix_complete':args.stage is None,
         'trajectories_verified':len(records),'new_simulations':done['new_simulations'],
         'whole_independent_reused_groups':8,'model_refits':0,'all_seed_grid_step_combinations_match':True,
         'saved_predictor_reproduction_max_errors':errors,'records':records,
         'results_sha256':sha(HERE/f'results{suffix}.json'),'verifier_sha256':sha(Path(__file__))}
    (HERE/f'verification{suffix}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'}))


if __name__=='__main__':main()
