"""Reuse the original marker-stress trajectories without inventing independent runs.

These files provide geometry-only compatibility and resolution diagnostics. They
do NOT supply an independent-run predictive test or true tangent FTLE labels.
"""
import json
import numpy as np
from features import history,neighbors,wrap,index_at
from simulate import ROOT,sha,save_json


def local_stretch(a,b,k,horizon):
    idx=neighbors(a,k)
    r=wrap(a[idx]-a[:,None]); q=wrap(b[idx]-b[:,None])
    gram=np.einsum("nki,nkj->nij",r,r)
    cross=np.einsum("nki,nkj->nij",r,q)
    reg=1e-10*np.maximum(np.trace(gram,axis1=1,axis2=2),1e-12)[:,None,None]
    B=np.linalg.solve(gram+reg*np.eye(3),cross)
    return np.log(np.linalg.svd(B,compute_uv=False)[:,0]) / horizon


def main():
    folder=ROOT.parent/"marker-stress"
    configs=["reference","tracking_coarse","tracking_medium","fluid_half","grid128_medium","grid128","grid128_refined","saved_time_coarse"]
    out={"independent_initial_conditions":1,"predictive_train_test_comparison_valid":False,
         "reason":"All eight files share one physical initial-condition family. Companion markers/clones are not independent fluid runs. Grid variants also sample the raw nonperiodic start differently.",
         "reused_particles":"First 250 primary markers; exclude 250 displaced companions and 5 exact clones",
         "history_window":.04,"gap":.02,"horizon":.04,"anchors":[.10,.16,.22],"neighbors":8,
         "target":"Finite-neighborhood affine logarithmic stretch, not tangent FTLE; current start-neighbor set fixed through target window",
         "historical_features":"Distance, angle, exchange and least-squares deformation evaluated using the new feature extractor. Alignment requiring unavailable gradient is omitted from summary.",
         "records":{}}
    targets={}
    for name in configs:
        path=folder/(name+".npz")
        with np.load(path,allow_pickle=False) as d:
            ts=d["time"]; x=d["positions"][:,:250]; u=d["velocities"][:,:250]
        hs=[]; ys=[]
        for t in out["anchors"]:
            i=index_at(ts,t); b=index_at(ts,t-.04)
            a=index_at(ts,t+.02); e=index_at(ts,t+.06)
            H=history(x[b:i+1],u[i],np.zeros((250,3,3)),8,.04)
            hs.append(H); ys.append(local_stretch(x[a],x[e],8,.04))
        H=np.concatenate(hs); y=np.concatenate(ys); targets[name]=y
        out["records"][name]={"file":str(path),"sha256":sha(path),"n_rows":len(y),
                               "target_mean":float(y.mean()),"target_std":float(y.std()),
                               "mean_neighbor_endpoint_exchange":float(H[:,9].mean()),
                               "mean_distance_log_change":float(H[:,0].mean())}
    for name in configs:
        out["records"][name]["target_RMSE_vs_reference"]=float(np.sqrt(np.mean((targets[name]-targets["reference"])**2)))
    out["source_verification_sha256"]=sha(folder/"verification.json")
    (ROOT/"results").mkdir(exist_ok=True)
    save_json(ROOT/"results"/"existing_data_audit.json",out)
    print(json.dumps({"reused_trajectories":len(configs),"independent_initial_conditions":1,
                      "grid128_refined_target_RMSE":out["records"]["grid128_refined"]["target_RMSE_vs_reference"]}))


if __name__=="__main__":main()
