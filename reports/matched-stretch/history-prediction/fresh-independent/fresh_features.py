"""Reuse the causal parent descriptors; add explicitly diagnostic mirror controls."""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from features import dataset,history,index_at,wrap
from evaluate import matrices


def make_dataset(paths,p,options=None):
    options=options or {}
    options={'lookback':p['primary_lookback'],'horizon':p['primary_horizon'],
             'k':p['neighbors'],'count':p['base']['markers'],'stride':1,**options}
    d=dataset(paths,p,**options)
    mirrored=[];reversed_mirror=[]
    for path in paths:
        with np.load(path,allow_pickle=False) as z:
            time=z['time'];positions=z['positions'][:,:options['count']]
            velocity=z['velocity'][:,:options['count']];gradient=z['gradient'][:,:options['count']]
        for t in p['anchors']:
            i=index_at(time,t);b=index_at(time,t-options['lookback'])
            # Reflect the recorded history about the periodic origin, keeping
            # current fluid inputs and target fixed. This is not fluid evolution.
            window=positions[b:i+1:options['stride']].copy();window[:,:,0]*=-1
            window=wrap(window)
            for reverse,output in [(False,mirrored),(True,reversed_mirror)]:
                output.append(history(window,velocity[i],gradient[i],options['k'],options['lookback'],reverse=reverse))
    d['mirrored']=np.concatenate(mirrored);d['mirrored_reverse']=np.concatenate(reversed_mirror)
    return d


def make_matrices(d):
    ret=matrices(d)
    ret['B_mirrored_history']=np.column_stack([d['present'],d['mirrored']])
    ret['B_mirrored_reversed']=np.column_stack([d['present'],d['mirrored_reverse']])
    return ret


def cases(p):
    result=[{'name':'primary','options':{},'quantile':p['primary_event_quantile']}]
    for config in p['observation_sensitivities']:
        result.append({'name':config['name'],'options':{k:v for k,v in config.items() if k not in ['name','quantile']},
                       'quantile':config.get('quantile',p['primary_event_quantile'])})
    return result
