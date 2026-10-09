"""Run-held-out predictors, controls, calibration and paired run-level uncertainty."""
from pathlib import Path
import argparse
import itertools
import json
import os
import time
os.environ.setdefault("OMP_NUM_THREADS","1")
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.metrics import mean_absolute_error, r2_score, average_precision_score, brier_score_loss, log_loss, precision_score, recall_score, precision_recall_curve
from features import dataset,shuffle_blocks,PRESENT_NAMES,HISTORY_NAMES
from simulate import ROOT,sha,save_json


def split_masks(data,p):
    masks=[np.isin(data["group"],p[key]) for key in ["train_seeds","validation_seeds","test_seeds"]]
    assert all(m.any() for m in masks)
    assert np.all(sum(m.astype(int) for m in masks)==1)
    assert not set(data["group"][masks[0]]) & set(data["group"][masks[2]])
    return masks


def matrices(d):
    A,H=d["present"],d["history"]
    weak=np.column_stack([A[:,:3],A[:,14],A[:,15],A[:,16],A[:,17]])
    # Exactly 23 additional current-only terms: random fixed selection of degree-2
    # terms, generated without fitting or observing targets. Equal B input count.
    poly=PolynomialFeatures(2,include_bias=False).fit_transform(A)
    ix=np.random.default_rng(513).choice(np.arange(A.shape[1],poly.shape[1]),H.shape[1],replace=False)
    past=np.array([i for i,n in enumerate(HISTORY_NAMES) if n.startswith("past_")])
    no_past=np.array([i for i in range(H.shape[1]) if i not in past])
    no_exchange=np.array([i for i,n in enumerate(HISTORY_NAMES) if not n.startswith("neighbor_")])
    ret={"A_full":A,"B_history":np.column_stack([A,H]),"A_weak":weak,"B_weak":np.column_stack([weak,H]),
         "A_matched_current":np.column_stack([A,poly[:,ix]]),"history_only":H,
         "B_past_only":np.column_stack([A,H[:,past]]),"B_no_past":np.column_stack([A,H[:,no_past]]),
         "B_no_exchange":np.column_stack([A,H[:,no_exchange]]),
         "B_reversed":np.column_stack([A,d["reverse"]]),
         "B_irrelevant":np.column_stack([A,d["irrelevant"]]),
         "B_permuted_labels":np.column_stack([A,H])}
    for seed in [91,92,93]:
        ret[f"B_shuffled_{seed}"]=np.column_stack([A,shuffle_blocks(H,d["group"],d["anchor"],seed)])
    return ret


def macro_loss(y,p,g,classification=False):
    return float(np.mean([log_loss(y[g==v],p[g==v],labels=[0,1]) if classification else
                          np.sqrt(np.mean((y[g==v]-p[g==v])**2)) for v in np.unique(g)]))


def fit_pair(X,y,event,masks,groups,family="linear",permuted=None):
    tr,va,te=masks
    output={}
    fitted=[]
    for classification,target in [(False,y),(True,event)]:
        train_target=target[tr] if permuted is None else permuted[classification][tr]
        candidates=([.01,.1,1.,10.] if classification else [.1,1.,10.,100.,1000.]) if family=="linear" else [1.,10.]
        best=None; records=[]
        for param in candidates:
            if family=="linear":
                est=LogisticRegression(C=param,max_iter=500,solver="lbfgs",tol=1e-6) if classification else Ridge(alpha=param)
                model=make_pipeline(StandardScaler(),est)
            else:
                cls=HistGradientBoostingClassifier if classification else HistGradientBoostingRegressor
                model=cls(max_iter=80,max_leaf_nodes=7,max_depth=3,learning_rate=.08,min_samples_leaf=40,
                          l2_regularization=param,early_stopping=False,random_state=26)
            model.fit(X[tr],train_target)
            predict=lambda xx: model.predict_proba(xx)[:,1] if classification else model.predict(xx)
            vp=predict(X[va])
            loss=macro_loss(target[va],vp,groups[va],classification)
            records.append({"regularization":param,"validation_loss":loss})
            if best is None or loss<best[0]: best=(loss,param,model,vp)
        loss,param,model,vp=best
        pred=model.predict_proba(X[te])[:,1] if classification else model.predict(X[te])
        key="event" if classification else "continuous"
        output[key]={"selected_regularization":param,"validation_loss":loss,"candidates":records}
        output["probability" if classification else "prediction"]=pred
        if classification:
            precision,recall,thresholds=precision_recall_curve(event[va],vp)
            f1=2*precision[:-1]*recall[:-1]/np.maximum(precision[:-1]+recall[:-1],1e-12)
            output["decision_threshold"]=float(thresholds[np.argmax(f1)])
        else: output["validation_prediction"]=vp
        if family=="linear":
            co=model.steps[-1][1].coef_.reshape(-1)
            output[key]["standardized_coefficients"]=co.tolist()
        fitted.append(model)
    return output,fitted


def calibration(y,p):
    rows=[]; ece=0.
    bins=np.minimum((p*10).astype(int),9)
    for b in range(10):
        use=bins==b
        if use.any():
            predicted=float(p[use].mean()); observed=float(y[use].mean()); n=int(use.sum())
            ece+=n/len(y)*abs(predicted-observed)
            rows.append({"bin":b,"count":n,"probability":predicted,"observed":observed})
    return float(ece),rows


def metrics(y,event,pred,prob,threshold):
    ece,_=calibration(event,prob)
    return {"RMSE":float(np.sqrt(np.mean((y-pred)**2))),"MAE":float(mean_absolute_error(y,pred)),
            "R2":float(r2_score(y,pred)),"AUPRC_AP":float(average_precision_score(event,prob)),
            "Brier":float(brier_score_loss(event,prob)),"log_loss":float(log_loss(event,prob,labels=[0,1])),
            "ECE10":ece,"precision":float(precision_score(event,prob>=threshold,zero_division=0)),
            "recall":float(recall_score(event,prob>=threshold,zero_division=0)),"prevalence":float(event.mean())}


def summarize(y,event,pred,prob,group,threshold):
    by={str(g):metrics(y[group==g],event[group==g],pred[group==g],prob[group==g],threshold) for g in np.unique(group)}
    return {"pooled":metrics(y,event,pred,prob,threshold),"per_run":by,
            "macro":{k:float(np.mean([v[k] for v in by.values()])) for k in next(iter(by.values()))},
            "calibration":calibration(event,prob)[1]}


def paired(a,b,p):
    keys=sorted(a["per_run"])
    assert keys==sorted(b["per_run"])
    rng=np.random.default_rng(p["bootstrap"]["seed"])
    draw=rng.integers(0,len(keys),(p["bootstrap"]["replicates"],len(keys)))
    out={}
    for metric in a["macro"]:
        delta=np.array([b["per_run"][key][metric]-a["per_run"][key][metric] for key in keys])
        boot=delta[draw].mean(1)
        value={"mean_B_minus_A":float(delta.mean()),"ci95":np.quantile(boot,[.025,.975]).tolist(),"paired_run_deltas":delta.tolist()}
        if metric in ["RMSE","AUPRC_AP"]:
            flips=np.array(list(itertools.product([-1,1],repeat=len(keys))))
            value["two_sided_signflip_p_unadjusted"]=float(np.mean(abs((flips*delta).mean(1))>=abs(delta.mean())-1e-14))
        out[metric]=value
    return out


def paths_for(data,seeds,n=26,nu=.02,dt=.005):
    paths=[Path(data)/f"s{s}_n{n}_nu{nu}_dt{dt}.npz" for s in seeds]
    if not all(p.exists() for p in paths): raise FileNotFoundError("Missing required simulation(s)")
    return paths


def cached_dataset(data,p,options,paths=None):
    paths=paths or paths_for(data,p["train_seeds"]+p["validation_seeds"]+p["test_seeds"])
    # Content-derived cache includes all simulation hashes and feature source hash.
    import hashlib
    ident=json.dumps({"files":{str(x):sha(x) for x in paths},"options":options,"protocol":p,"features":sha(ROOT/"features.py")},sort_keys=True)
    digest=hashlib.sha256(ident.encode()).hexdigest()
    folder=ROOT/"cache"; folder.mkdir(exist_ok=True)
    path=folder/(digest+".npz")
    if path.exists():
        with np.load(path,allow_pickle=False) as d: return dict(d)
    d=dataset(paths,p,**options)
    np.savez_compressed(path,**d)
    return d


def experiment(d,p,quantile=.9,all_models=True):
    masks=split_masks(d,p); tr,va,te=masks
    y=d["target"]; cutoff=float(np.quantile(y[tr],quantile)); event=(y>=cutoff).astype(int)
    Xs=matrices(d)
    names=list(Xs) if all_models else ["A_full","B_history"]
    results={"train_event_cutoff":cutoff,"quantile":quantile,"n_train":int(tr.sum()),"n_validation":int(va.sum()),
             "n_test":int(te.sum()),"models":{},"paired":{}}
    predictions={"target":y[te],"event":event[te],"group":d["group"][te],"anchor":d["anchor"][te],
                 "particle":d["particle"][te],"position":d["position"][te]}
    fitted={}
    for name in names+(["tree_A_full","tree_B_history"] if all_models else []):
        tree=name.startswith("tree_"); xname=name[5:] if tree else name
        perm=None
        if name=="B_permuted_labels":
            yp=shuffle_blocks(y,d["group"],d["anchor"],61)
            ep=shuffle_blocks(event,d["group"],d["anchor"],61)
            perm={False:yp,True:ep}
        fit,model=fit_pair(Xs[xname],y,event,masks,d["group"],family="tree" if tree else "linear",permuted=perm)
        pred=fit.pop("prediction"); prob=fit.pop("probability"); fit.pop("validation_prediction")
        fit.update(summarize(y[te],event[te],pred,prob,d["group"][te],fit["decision_threshold"]))
        fit["input_dimensions"]=Xs[xname].shape[1]
        results["models"][name]=fit
        predictions[name+"_prediction"]=pred; predictions[name+"_probability"]=prob
        fitted[name]=model
        print(f"  {name}: RMSE={fit['macro']['RMSE']:.5f}, AP={fit['macro']['AUPRC_AP']:.5f}",flush=True)
    for name in results["models"]:
        if name!="A_full": results["paired"][name+"_vs_A_full"]=paired(results["models"]["A_full"],results["models"][name],p)
    if all_models:
        results["paired"]["tree_B_vs_tree_A"]=paired(results["models"]["tree_A_full"],results["models"]["tree_B_history"],p)
        results["paired"]["B_vs_matched_current"]=paired(results["models"]["A_matched_current"],results["models"]["B_history"],p)
        results["paired"]["B_weak_vs_A_weak"]=paired(results["models"]["A_weak"],results["models"]["B_weak"],p)
    return results,predictions,fitted


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--data",type=Path,default=ROOT/"data"); ap.add_argument("--quick",action="store_true")
    args=ap.parse_args(); p=json.loads((ROOT/"protocol.json").read_text())
    out=ROOT/"results"; out.mkdir(exist_ok=True)
    options={"lookback":.2,"horizon":.2,"k":8,"count":512,"stride":1}
    start=time.perf_counter()
    d=cached_dataset(args.data,p,options)
    main_result,pred,models=experiment(d,p)
    save_json(out/"primary.json",main_result); np.savez_compressed(out/"predictions.npz",**pred)
    save_json(out/"feature_names.json",{"present":PRESENT_NAMES,"history":HISTORY_NAMES})
    if not args.quick:
        sensitivity={}
        configurations=[("lookback_0.1",{"lookback":.1}), ("lookback_0.3",{"lookback":.3}),
                        ("horizon_0.1",{"horizon":.1}), ("horizon_0.4",{"horizon":.4}),
                        ("neighbors_4",{"k":4}), ("neighbors_16",{"k":16}),
                        ("markers_256",{"count":256}), ("history_stride_2",{"stride":2})]
        for name,change in configurations:
            print(name,flush=True)
            ds=cached_dataset(args.data,p,{**options,**change})
            res,_,_=experiment(ds,p,all_models=False); sensitivity[name]=res
            save_json(out/"sensitivity.json",sensitivity)
        for q in [.8,.95]:
            print("threshold "+str(q),flush=True)
            res,_,_=experiment(d,p,quantile=q,all_models=False); sensitivity["quantile_"+str(q)]=res
        # All grid/timestep/viscosity changes are evaluated with frozen primary fits.
        for name,change in [("grid_38",{"n":38}),("step_half",{"dt":.0025}),
                            ("viscosity_0.01",{"nu":.01}),("viscosity_0.04",{"nu":.04})]:
            print(name,flush=True)
            ds=cached_dataset(args.data,p,options,paths_for(args.data,p["test_seeds"],**change))
            xs=matrices(ds); yy=ds["target"]; ee=(yy>=main_result["train_event_cutoff"]).astype(int)
            res={"frozen_primary_models":True,"n_test":len(yy),"models":{}}
            for m in ["A_full","B_history"]:
                pp=models[m][0].predict(xs[m]); prob=models[m][1].predict_proba(xs[m])[:,1]
                res["models"][m]=summarize(yy,ee,pp,prob,ds["group"],main_result["models"][m]["decision_threshold"])
            res["paired"]={"B_history_vs_A_full":paired(res["models"]["A_full"],res["models"]["B_history"],p)}
            res["paired_target_RMSE_from_base"]=float(np.sqrt(np.mean((yy-pred["target"])**2)))
            sensitivity[name]=res
        save_json(out/"sensitivity.json",sensitivity)
    save_json(out/"execution.json",{"seconds":time.perf_counter()-start,"source_hashes":{p.name:sha(p) for p in ROOT.glob("*.py")},
                                      "protocol_sha256":sha(ROOT/"protocol.json")})
    print("Evaluation complete",flush=True)


if __name__=="__main__":main()
