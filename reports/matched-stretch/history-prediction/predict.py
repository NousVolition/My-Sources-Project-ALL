"""Apply saved A/B models to a past-only observation NPZ, producing marker CSV.

Input keys: time (T,), positions (T,N,3), velocity (N,3), gradient (N,3,3).
Positions must contain ONLY the .2 lookback window up to the forecast instant.
Neither future positions nor marker identifiers are accepted or needed.
"""
import argparse
import csv
import json
import numpy as np
from features import present,history
from export_models import predict_model
from simulate import ROOT


def forecast(observation,models):
    ts=np.asarray(observation["time"]); x=np.asarray(observation["positions"])
    u=np.asarray(observation["velocity"]); J=np.asarray(observation["gradient"])
    if len(ts)<3 or not np.all(np.diff(ts)>0) or not np.isclose(ts[-1]-ts[0],models["lookback"]):
        raise ValueError("Supply strictly increasing times for the complete trained lookback")
    if len(ts)!=9 or not np.allclose(np.diff(ts),.025):
        raise ValueError("The saved primary model requires 9 observations spaced by 0.025")
    if x.shape!=(len(ts),len(u),3) or J.shape!=(len(u),3,3) or u.shape!=(len(u),3):
        raise ValueError("Observation shapes do not agree")
    if not all(np.isfinite(a).all() for a in [ts,x,u,J]): raise ValueError("Nonfinite observations")
    A=present(x[-1],u,J,models["neighbors"])[0]
    H=history(x,u,J,models["neighbors"],models["lookback"])
    out={}
    for name,X in [("A_full",A),("B_history",np.column_stack([A,H]))]:
        model=models["models"][name]
        out[name+"_ftle"]=predict_model(model["continuous"],X)
        out[name+"_probability"]=predict_model(model["event"],X)
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument("observation");ap.add_argument("output")
    ap.add_argument("--models",default=str(ROOT/"results"/"primary_models.json"));args=ap.parse_args()
    models=json.loads(open(args.models,encoding="utf-8").read())
    with np.load(args.observation,allow_pickle=False) as data:
        output=forecast(data,models);pos=data["positions"][-1];t=data["time"][-1]
    with open(args.output,"w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["row","prediction_time","target_start","target_end","x","y","z"]+list(output))
        for i in range(len(pos)):
            w.writerow([i,t,t+models["gap"],t+models["gap"]+models["horizon"],*pos[i],*[v[i] for v in output.values()]])


if __name__=="__main__":main()
