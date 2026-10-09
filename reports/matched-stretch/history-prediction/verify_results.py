"""Integrity, numerical sensitivity and identity audits on generated trajectories."""
import json
import numpy as np
from simulate import ROOT,sha,save_json
from features import present,history,index_at,future_target,wrap


def main():
    p=json.loads((ROOT/"protocol.json").read_text())
    rows=[]
    for path in sorted((ROOT/"data").glob("*.npz")):
        m=json.loads(path.with_suffix(".json").read_text())
        assert m["sha256"]==sha(path)
        assert m["sources"]["numerics.py"]==sha(ROOT.parent/"numerics.py")
        assert m["sources"]["simulate.py"]==sha(ROOT/"simulate.py")
        assert m["sources"]["protocol.json"]==sha(ROOT/"protocol.json")
        group=m["seed"]
        partition=next(key for key in ["train_seeds","validation_seeds","test_seeds"] if group in p[key])
        diag=m["diagnostics"]
        rows.append({"file":path.name,"sha256":m["sha256"],"independent_group":group,"split":partition,
                     "config":m["config"],"max_CFL":m["max_cfl"],"seconds":m["seconds"],
                     "max_abs_energy_residual":max(abs(r["energy_budget_residual"]) for r in diag),
                     "max_divergence":max(r["divergence_max"] for r in diag),
                     "max_high_band_energy":max(r["high_band_energy_fraction"] for r in diag),
                     "max_volume_error_from_time_zero":max(r["max_det_error"] for r in diag)})
    assert len(rows)==64 and len({r["independent_group"] for r in rows})==32
    comparisons={}
    for label,n,dt in [("grid",38,.005),("time",26,.0025)]:
        summary=[]
        for seed in p["test_seeds"]:
            a=np.load(ROOT/"data"/f"s{seed}_n26_nu0.02_dt0.005.npz",allow_pickle=False)
            b=np.load(ROOT/"data"/f"s{seed}_n{n}_nu0.02_dt{dt}.npz",allow_pickle=False)
            diff=np.linalg.norm(wrap(a["positions"]-b["positions"]),axis=-1)
            dy=[]; det_error=[]
            for t in p["anchors"]:
                i=index_at(a["time"],t+p["gap"]); j=index_at(a["time"],t+p["gap"]+.2)
                dy.append(future_target(a["tangent"][i],a["tangent"][j],.2)-future_target(b["tangent"][i],b["tangent"][j],.2))
                det_error.extend(abs(np.linalg.det(b["tangent"][j])/np.linalg.det(b["tangent"][i])-1))
            summary.append({"seed":seed,"path_RMS":float(np.sqrt(np.mean(diff**2))),"path_max":float(diff.max()),
                            "target_RMSE":float(np.sqrt(np.mean(np.array(dy)**2))),
                            "refined_target_volume_error_p95":float(np.quantile(det_error,.95)),
                            "refined_target_volume_error_max":float(np.max(det_error))})
            a.close();b.close()
        comparisons[label]=summary
    invariance=[]
    for seed in p["test_seeds"][:2]:
        with np.load(ROOT/"data"/f"s{seed}_n26_nu0.02_dt0.005.npz",allow_pickle=False) as d:
            i=index_at(d["time"],.5); b=index_at(d["time"],.3)
            x=d["positions"][b:i+1]; u=d["velocity"][i]; J=d["gradient"][i]
        perm=np.random.default_rng(seed).permutation(len(u))
        A=present(x[-1],u,J,8)[0]; H=history(x,u,J,8,.2)
        Ap=present(x[-1,perm],u[perm],J[perm],8)[0]; Hp=history(x[:,perm],u[perm],J[perm],8,.2)
        error=float(max(abs(Ap-A[perm]).max(),abs(Hp-H[perm]).max()))
        assert error<1e-10
        invariance.append({"seed":seed,"max_feature_error_after_ID_permutation":error})
    out={"status":"passed","simulation_files":len(rows),"independent_seeds":32,
         "source_commit":"17c14e89940d4c96ec9df315b1a689430d8c5ddf","parent_solver_sha256":sha(ROOT.parent/"numerics.py"),
         "files":rows,"refinement":comparisons,"identity_invariance":invariance,
         "limits":"Passed checks verify implementation/integrity, not continuum convergence. Trilinear interpolation is not exactly divergence-free and numerical FTLE inherits that error."}
    save_json(ROOT/"results"/"verification.json",out)
    print(json.dumps({k:out[k] for k in ["status","simulation_files","independent_seeds","identity_invariance"]}))


if __name__=="__main__":main()
