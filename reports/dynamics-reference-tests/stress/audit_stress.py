"""Save a reference-derivative diagnostic and reproduce one finest-step long run."""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
from run_stress import PLAN, transport_exact, fingerprint
from core import lyapunov

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"


def derivative_diagnostic():
    theta=np.random.default_rng(44).uniform(0,2*np.pi,32)
    h=1e-5
    t=.37
    angular=(transport_exact(theta-2*h,t)-8*transport_exact(theta-h,t)
             +8*transport_exact(theta+h,t)-transport_exact(theta+2*h,t))/(12*h)
    settings=PLAN["transport"]
    source=np.full(len(theta),settings["forcing_mean"])
    for k,amp in settings["forcing_cosine_amplitudes"].items(): source+=amp*np.cos(int(k)*theta)
    rhs=source-settings["leak"]*transport_exact(theta,t)-settings["omega"]*angular
    rows=[]
    for dt in (1e-5,5e-6,2.5e-6):
        numerical=(transport_exact(theta,t+dt)-transport_exact(theta,t-dt))/(2*dt)
        rows.append(dict(time_difference_step=dt,max_residual=float(np.max(abs(numerical-rhs)))))
    result=dict(initial_test_tolerance=1e-7,unchanged_final_tolerance=1e-7,values=rows,
                explanation="The first independent check failed marginally. Refining its second-order time difference reduces the residual; the tested model and acceptance tolerance were unchanged. The initial failure XML is retained.")
    (DATA/"reference-derivative-diagnostic.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2),flush=True)


def selected_reproduction():
    started=time.perf_counter()
    settings=PLAN["long_run"]
    seed=settings["seeds"][0]
    h=settings["steps"][-1]
    path=DATA/f"long_seed{seed}_h{h}.npz"
    with np.load(path,allow_pickle=False) as saved:
        if saved["fingerprint"].item()!=fingerprint()[0]: raise ValueError("Original code/plan fingerprint changed")
        spectrum,blocks,states=lyapunov(saved["initial"],h=h,transient=settings["transient"],
                                       duration=settings["duration"],block=settings["block"])
        comparisons={name:bool(np.array_equal(value,saved[name])) for name,value in
                     [("spectrum",spectrum),("block_rates",blocks),("states",states)]}
        sizes={name:list(saved[name].shape) for name in comparisons}
    result=dict(seed=seed,step=h,duration=settings["duration"],arrays_exactly_equal=comparisons,array_shapes=sizes,
                passed=all(comparisons.values()),elapsed_seconds=time.perf_counter()-started,
                raw_file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope="One of 24 extended trajectories independently rerun with the same code and environment, including 10000 sampled states and 20 tangent-rate blocks. This is not a repeat of all extended runs.")
    (DATA/"selected-reproduction.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2),flush=True)
    if not result["passed"]: raise RuntimeError("Selected reproduction failed")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--derivative-only",action="store_true")
    args=parser.parse_args()
    derivative_diagnostic()
    if not args.derivative_only: selected_reproduction()
