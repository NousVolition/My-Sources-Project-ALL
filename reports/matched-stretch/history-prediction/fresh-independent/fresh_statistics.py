"""Whole-run inference for spatial event ranking; no exponential sign enumeration."""
import numpy as np
from sklearn.metrics import average_precision_score


def spatial_scores(event,prob,group,anchor):
    out={}
    for seed in np.unique(group):
        rows=[]
        for t in np.unique(anchor):
            use=(group==seed)&(anchor==t);y=event[use];q=prob[use]
            n=int(y.sum());k=int(np.ceil(.1*len(y)))
            # Same fixed marker budget for both models; stable ordering only
            # settles exactly tied scores and never uses target labels.
            chosen=np.argsort(-q,kind='stable')[:k]
            hits=int(y[chosen].sum())
            rows.append({'anchor':float(t),'events':n,'markers':len(y),
                'AP':float(average_precision_score(y,q)) if n else 0.,
                'precision_at_10percent':hits/k,'recall_at_10percent':hits/n if n else 0.,
                'predicted_marker_budget':k})
        out[str(int(seed))]={'mean_AP':float(np.mean([r['AP'] for r in rows])),
                           'mean_precision_at_10percent':float(np.mean([r['precision_at_10percent'] for r in rows])),
                           'mean_recall_at_10percent':float(np.mean([r['recall_at_10percent'] for r in rows])),
                           'anchors':rows}
    return out


def paired_values(values,p):
    values=np.asarray(values,dtype=float)
    rng=np.random.default_rng(p['bootstrap']['seed'])
    ix=rng.integers(0,len(values),(p['bootstrap']['replicates'],len(values)))
    return {'mean':float(values.mean()),'ci95':np.quantile(values[ix].mean(1),[.025,.975]).tolist(),
            'per_run':values.tolist(),'independent_runs':len(values)}


def compare(a,b,p):
    keys=sorted(a['per_run'],key=int);assert keys==sorted(b['per_run'],key=int)
    out={metric:paired_values([b['per_run'][k][metric]-a['per_run'][k][metric] for k in keys],p)
         for metric in a['macro']}
    out['spatial_AP']=paired_values([b['spatial'][k]['mean_AP']-a['spatial'][k]['mean_AP'] for k in keys],p)
    out['relative_RMSE_reduction']=paired_values([1-b['per_run'][k]['RMSE']/a['per_run'][k]['RMSE'] for k in keys],p)
    out['group_order']=keys
    return out


def classify(interval,useful,numerically_audited):
    if not numerically_audited:return 'numerically unresolved; no usefulness or exclusion claim'
    lo,hi=interval
    if lo>useful:return 'useful benefit supported for the predeclared comparison'
    if hi<useful:
        return 'useful-size benefit ruled out for the predeclared comparison; small positive effect remains possible' if hi>=0 else 'useful-size benefit ruled out; evidence of worse average prediction'
    if lo>0:return 'positive effect detected; useful size remains unresolved'
    return 'unresolved: interval permits no benefit and a useful-size benefit'
