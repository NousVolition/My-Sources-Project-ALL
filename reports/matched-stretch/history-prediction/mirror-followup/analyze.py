"""Reproduce a bounded mirror/chirality follow-up without running fluid dynamics."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse
import json
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parent))
from chirality import signed_features
from features import history,index_at,future_target,shuffle_blocks
from evaluate import paths_for,cached_dataset,split_masks,fit_pair,summarize,paired
from simulate import sha,save_json


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',type=Path,default=HERE.parent/'recorded-data')
    args=ap.parse_args()
    protocol=json.loads((HERE/'protocol.json').read_text())
    p=json.loads((HERE.parent/'protocol.json').read_text())
    paths=paths_for(args.data,p['train_seeds']+p['validation_seeds']+p['test_seeds'])
    provenance={x.name:sha(x) for path in paths for x in (path,path.with_suffix('.json'))}
    source_paths=[HERE/'analyze.py',HERE/'chirality.py',HERE/'protocol.json',HERE.parent/'features.py',HERE.parent/'evaluate.py',HERE.parent/'protocol.json',HERE.parent/'simulate.py',HERE.parent.parent/'numerics.py']
    sources={str(x.relative_to(HERE.parent.parent)).replace('\\','/'):sha(x) for x in source_paths}
    d=cached_dataset(args.data,p,{'lookback':.2,'horizon':.2,'k':8,'count':512,'stride':1},paths)
    cs=[]; hs=[]
    audit={k:0. for k in ['history_mirror','reverse_history_mirror','current_chirality_mirror','past_chirality_mirror','past_chirality_reverse','particle_permutation','target_mirror']}
    R=np.diag([-1,1,1]); block=0
    def update(name,a,b): audit[name]=max(audit[name],float(np.max(abs(a-b))))
    for path in paths:
        with np.load(path,allow_pickle=False) as z:
            for t in p['anchors']:
                i=index_at(z['time'],t); b=index_at(z['time'],t-.2)
                x=z['positions'][b:i+1]; u=z['velocity'][i]; J=z['gradient'][i]
                c,h=signed_features(x); cs.append(c); hs.append(h)
                ix=slice(block*512,(block+1)*512); block+=1
                update('history_mirror',history(x@R,u@R,R@J@R,8,.2),d['history'][ix])
                update('reverse_history_mirror',history(x@R,u@R,R@J@R,8,.2,reverse=True),d['reverse'][ix])
                cm,hm=signed_features(x@R)
                update('current_chirality_mirror',cm,-c)
                update('past_chirality_mirror',hm,h*[-1,1,1])
                _,hr=signed_features(x,reverse=True)
                update('past_chirality_reverse',hr,h*[-1,1,-1])
                perm=np.random.default_rng(block+900).permutation(512)
                cp,hp=signed_features(x[:,perm])
                update('particle_permutation',np.column_stack([cp,hp]),np.column_stack([c,h])[perm])
                a=index_at(z['time'],t+p['gap']); e=index_at(z['time'],t+p['gap']+.2)
                update('target_mirror',future_target(R@z['tangent'][a]@R,R@z['tangent'][e]@R,.2),d['target'][ix])
        print('Audited '+path.stem,flush=True)
    c=np.concatenate(cs); h=np.concatenate(hs)
    B=np.column_stack([d['present'],d['history']]); BC=np.column_stack([B,c])
    HR=h*[-1,1,-1]
    BR=np.column_stack([d['present'],d['reverse'],c])
    Xs={'B_original':B,'B_current_chirality':BC,
        'C_history_chirality':np.column_stack([BC,h]),
        'C_current_capacity':np.column_stack([BC,c**2,c**3,c**4]),
        'C_shuffled':np.column_stack([BC,shuffle_blocks(h,d['group'],d['anchor'],604)]),
        'C_mirrored_history_refit':np.column_stack([BC,h*[-1,1,-1]]),
        'C_reversed':np.column_stack([BR,HR]),
        'C_reversed_mirrored_history_refit':np.column_stack([BR,HR*[-1,1,-1]])}
    masks=split_masks(d,p); tr,va,te=masks
    y=d['target']; cutoff=float(np.quantile(y[tr],.9)); event=(y>=cutoff).astype(int)
    result={'protocol':protocol,'data_sha256':provenance,'source_sha256':sources,'audit_max_abs_error':audit,
            'audited_windows':block,'audited_marker_windows':len(y),'cutoff':cutoff,'models':{},'paired':{},'equivalence':{}}
    saved={'target':y[te],'event':event[te],'group':d['group'][te],'anchor':d['anchor'][te],'particle':d['particle'][te]}
    for name,X in Xs.items():
        fit,_=fit_pair(X,y,event,masks,d['group'])
        pred=fit.pop('prediction'); prob=fit.pop('probability'); fit.pop('validation_prediction')
        fit.update(summarize(y[te],event[te],pred,prob,d['group'][te],fit['decision_threshold']))
        fit['input_dimensions']=X.shape[1]; result['models'][name]=fit
        saved[name+'_prediction']=pred; saved[name+'_probability']=prob
        print(name,fit['macro']['RMSE'],fit['macro']['AUPRC_AP'],flush=True)
    for base in ['B_current_chirality','C_current_capacity','C_shuffled']:
        result['paired']['C_history_chirality_vs_'+base]=paired(result['models'][base],result['models']['C_history_chirality'],p)
    for original,mirrored in [('C_history_chirality','C_mirrored_history_refit'),('C_reversed','C_reversed_mirrored_history_refit')]:
        result['equivalence'][original]={key:float(np.max(abs(saved[original+'_'+key]-saved[mirrored+'_'+key]))) for key in ['prediction','probability']}
    result['new_feature_ranges']={'current_chirality':[float(c.min()),float(c.max())],
                                  'history_min':h.min(axis=0).tolist(),'history_max':h.max(axis=0).tolist()}
    assert max(audit.values())<1e-9, audit
    assert all(max(v.values())<1e-8 for v in result['equivalence'].values())
    save_json(HERE/'results.json',result)
    np.savez_compressed(HERE/'predictions.npz',**saved)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=['B_original','B_current_chirality','C_history_chirality','C_current_capacity','C_shuffled','C_reversed']
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    labels=['Original history','+ present handedness','+ past handedness','Current capacity control','Shuffled past handedness','Reversed history']
    for ax,metric,title in zip(axes,['RMSE','AUPRC_AP'],['Future stretch error (lower better)','Strong-event average precision (higher better)']):
        for i,name in enumerate(names):
            vals=[row[metric] for row in result['models'][name]['per_run'].values()]
            ax.scatter(vals,[i]*len(vals),s=16,alpha=.55,color='#447799')
            ax.scatter([np.mean(vals)],[i],marker='D',s=40,color='#c04a24')
        ax.set_yticks(range(len(names)),labels if ax==axes[0] else ['']*len(names))
        ax.invert_yaxis(); ax.set_title(title,fontsize=10); ax.grid(axis='x',alpha=.2)
    fig.suptitle('Exploratory mirror follow-up: eight previously inspected test runs\nDots: individual runs; diamonds: run means; grid 26 remains unresolved',fontsize=11)
    fig.savefig(HERE/'comparison.png',dpi=160); fig.savefig(HERE/'comparison.svg'); plt.close(fig)
    print('COMPLETE',flush=True)


if __name__=='__main__': main()
