"""Export primary linear fits as inspectable JSON; verify saved test predictions."""
import os
os.environ.setdefault("OMP_NUM_THREADS","1")
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
import json
import numpy as np
from scipy.special import expit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge,LogisticRegression
from simulate import ROOT,save_json
from features import PRESENT_NAMES,HISTORY_NAMES
from evaluate import cached_dataset,split_masks


def predict_model(model,X):
    z=(X-np.array(model["mean"]))/np.array(model["scale"])
    score=z@np.array(model["coefficient"])+model["intercept"]
    return expit(score) if model["kind"]=="logistic" else score


def main():
    p=json.loads((ROOT/"protocol.json").read_text())
    r=json.loads((ROOT/"results"/"primary.json").read_text())
    d=cached_dataset(ROOT/"data",p,{"lookback":.2,"horizon":.2,"k":8,"count":512,"stride":1})
    tr,va,te=split_masks(d,p)
    event=(d["target"]>=r["train_event_cutoff"]).astype(int)
    output={"lookback":.2,"horizon":.2,"gap":.025,"neighbors":8,"trained_marker_count":512,
            "event_target_cutoff":r["train_event_cutoff"],"models":{},"prediction_verification":{}}
    with np.load(ROOT/"results"/"predictions.npz",allow_pickle=False) as saved:
        for name,X,names in [("A_full",d["present"],PRESENT_NAMES),
                              ("B_history",np.column_stack([d["present"],d["history"]]),PRESENT_NAMES+HISTORY_NAMES)]:
            output["models"][name]={"feature_names":names,"decision_threshold":r["models"][name]["decision_threshold"]}
            for key,target,kind in [("continuous",d["target"],"ridge"),("event",event,"logistic")]:
                param=r["models"][name][key]["selected_regularization"]
                estimator=Ridge(alpha=param) if kind=="ridge" else LogisticRegression(C=param,max_iter=500,solver="lbfgs",tol=1e-6)
                model=make_pipeline(StandardScaler(),estimator).fit(X[tr],target[tr])
                scaler=model.steps[0][1];est=model.steps[-1][1]
                serial={"kind":kind,"regularization":param,"mean":scaler.mean_.tolist(),"scale":scaler.scale_.tolist(),
                        "coefficient":est.coef_.reshape(-1).tolist(),"intercept":float(np.asarray(est.intercept_).reshape(-1)[0])}
                predicted=predict_model(serial,X[te])
                original=saved[name+("_prediction" if kind=="ridge" else "_probability")]
                err=float(abs(predicted-original).max())
                assert err<1e-10
                output["models"][name][key]=serial
                output["prediction_verification"][name+"_"+key]=err
    save_json(ROOT/"results"/"primary_models.json",output)
    print(json.dumps(output["prediction_verification"]))


if __name__=="__main__":main()
