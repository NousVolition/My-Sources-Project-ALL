"""Fits train/validation only. It cannot open any test-flow observation."""
import argparse,json
from datetime import datetime,timezone
from pathlib import Path
from analyze import dataset,fit_model,shuffled
from run import ROOT,P,dump,digest

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'data');ap.add_argument('--out',type=Path,default=ROOT/'results');args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True);path=args.out/'frozen-models.json'
    if path.exists():raise RuntimeError('Refusing to overwrite frozen models')
    models={}
    for horizon in P['horizons']:
        train=dataset(args.data,P['train_seeds'],horizon);val=dataset(args.data,P['validation_seeds'],horizon)
        for column,target in enumerate(P['prediction']['targets']):
            key=f'{target}:T{horizon:g}';models[key]={}
            for family in ['current','history','capacity_current','shuffled_history','across_flow_history']:
                tr,va=train,val
                if family=='shuffled_history':tr,va=shuffled(train),shuffled(val)
                if family=='across_flow_history':tr,va=shuffled(train,True),shuffled(val,True)
                models[key][family]=fit_model(tr,va,column,family)
    dump(path,models)
    dump(args.out/'model-freeze-provenance.json',dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),models_sha256=digest(path),
        allowed_input_seeds=P['train_seeds']+P['validation_seeds'],forbidden_input_seeds=P['test_seeds'],test_observations_opened=False,
        source_sha256={f:digest(ROOT/f) for f in ['freeze_models.py','analyze.py','protocol.json','production-amendment.json']},
        note='Test trajectories may already exist, but this fitting routine never reads them. No test performance examined before freezing.'))
    print('Frozen 30 models from training and validation only.',flush=True)

if __name__=='__main__':main()
