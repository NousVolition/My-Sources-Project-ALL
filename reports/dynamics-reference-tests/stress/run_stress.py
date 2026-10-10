"""Run the predeclared local stress plan, retaining results and failed criteria."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import sys
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.stats import t as student_t

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from core import advance, integrate, lorenz, jacobian, lyapunov

PLAN=json.loads((ROOT/"stress_plan.json").read_text())
DATA=ROOT/"data"
DATA.mkdir(exist_ok=True)
CHECKS=[]


def plain(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    raise TypeError(type(x))


def write_json(name,value):
    (DATA/name).write_text(json.dumps(value,indent=2,default=plain,allow_nan=False)+"\n",encoding="utf-8")


def csv_out(name,rows):
    with (DATA/name).open("w",newline="",encoding="utf-8") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def check(name,condition,kind="accuracy",**metrics):
    CHECKS.append(dict(name=name,passed=bool(condition),kind=kind,**metrics))
    print(("PASS " if condition else "FAIL ")+name,flush=True)
    write_json("checks.json",CHECKS)


def fingerprint():
    paths=[ROOT/"run_stress.py",ROOT/"stress_plan.json",ROOT.parent/"core.py"]
    source={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    key=hashlib.sha256(json.dumps(source,sort_keys=True).encode()).hexdigest()
    return key,source


def lorenz_batch(y,rhos):
    x,v,z=y.T
    return np.column_stack((10*(v-x),x*(rhos-z)-v,x*v-(8/3)*z))


def spectrum_stats(rows,step,prefix):
    group=[r for r in rows if r["step"]==step and r["duration"]==prefix]
    rates=np.array([r["lambda1"] for r in group])
    half=student_t.ppf(.975,len(rates)-1)*rates.std(ddof=1)/np.sqrt(len(rates))
    return dict(step=step,duration=prefix,n=len(group),mean_largest_rate=rates.mean(),
                between_start_sd=rates.std(ddof=1),approximate_95_t_interval=[rates.mean()-half,rates.mean()+half],
                mean_z=np.mean([r["mean_z"] for r in group]),
                positive_x_fraction=np.mean([r["positive_x_fraction"] for r in group]),
                min_largest_rate=rates.min(),max_largest_rate=rates.max())


def long_runs(resume=False):
    settings=PLAN["long_run"]
    run_rows,prefix_rows=[],[]
    key,_=fingerprint()
    run_start=time.perf_counter()
    count=len(settings["steps"])*len(settings["seeds"])
    for h in settings["steps"]:
        for seed in settings["seeds"]:
            path=DATA/f"long_seed{seed}_h{h}.npz"
            initial=np.random.default_rng(seed).uniform([-.5,-.5,20.],[.5,.5,30.])
            started=time.perf_counter()
            reused=False
            if resume and path.exists():
                with np.load(path,allow_pickle=False) as saved:
                    if saved["fingerprint"].item()!=key:
                        raise ValueError(f"Cannot resume changed code/plan: {path.name}")
                    spectrum=saved["spectrum"].copy()
                    blocks=saved["block_rates"].copy()
                    states=saved["states"].copy()
                reused=True
            else:
                spectrum,blocks,states=lyapunov(initial,h=h,transient=settings["transient"],
                                               duration=settings["duration"],block=settings["block"])
                np.savez_compressed(path,spectrum=spectrum,block_rates=blocks,states=states,initial=initial,
                                    step=h,transient=settings["transient"],duration=settings["duration"],
                                    sample_interval=settings["qr_interval"],block_duration=settings["block"],fingerprint=key)
            trace_error=abs(spectrum.sum()+41/3)
            maximum=np.max(abs(states))
            row=dict(seed=seed,step=h,lambda1=spectrum[0],lambda2=spectrum[1],lambda3=spectrum[2],
                     trace_error=trace_error,max_abs_state=maximum,mean_z=states[:,2].mean(),
                     positive_x_fraction=np.mean(states[:,0]>0),wall_seconds=time.perf_counter()-started,reused=reused)
            run_rows.append(row)
            for duration in settings["prefix_durations"]:
                block_count=round(duration/settings["block"])
                state_count=round(duration/settings["qr_interval"])
                rates=blocks[:block_count].mean(axis=0)
                selected=states[:state_count]
                prefix_rows.append(dict(seed=seed,step=h,duration=duration,lambda1=rates[0],lambda2=rates[1],lambda3=rates[2],
                                        mean_z=selected[:,2].mean(),positive_x_fraction=np.mean(selected[:,0]>0)))
            check(f"Long-run invariants seed={seed}, h={h}",
                  np.all(np.isfinite(states)) and spectrum[0]>settings["per_run_minimum_largest_rate"]
                  and trace_error<settings["per_run_maximum_trace_error"]
                  and maximum<settings["per_run_maximum_absolute_state"],kind="long-run guard",
                  largest_rate=spectrum[0],trace_error=trace_error,max_abs_state=maximum)
            csv_out("long_runs.csv",run_rows)
            csv_out("duration_prefixes.csv",prefix_rows)
            elapsed=time.perf_counter()-run_start
            progress=dict(completed=len(run_rows),total=count,elapsed_seconds=elapsed,
                          simple_remaining_estimate_seconds=(elapsed/len(run_rows))*(count-len(run_rows)),last_seed=seed,last_step=h)
            write_json("progress.json",progress)
            print(f"PROGRESS {len(run_rows)}/{count}; {elapsed:.1f}s elapsed; last {row['wall_seconds']:.1f}s",flush=True)
    groups=[spectrum_stats(prefix_rows,h,duration) for h in settings["steps"] for duration in settings["prefix_durations"]]
    fine,finest=settings["steps"][-2:]
    at=lambda step,duration:next(x for x in groups if x["step"]==step and x["duration"]==duration)
    full=settings["duration"]
    a,b=at(fine,full),at(finest,full)
    tests=[("Rate stable across finest steps","mean_largest_rate",settings["finest_two_step_mean_rate_difference_tolerance"]),
           ("Mean z stable across finest steps","mean_z",settings["finest_two_step_mean_z_difference_tolerance"]),
           ("Lobe occupancy stable across finest steps","positive_x_fraction",settings["finest_two_step_positive_x_fraction_difference_tolerance"])]
    for name,metric,tolerance in tests:
        difference=abs(a[metric]-b[metric])
        check(name,difference<=tolerance,kind="predeclared robustness criterion",difference=difference,tolerance=tolerance)
    shorter=at(finest,500.)
    for name,metric,tolerance in [
        ("Rate stable from 500 to 1000 units","mean_largest_rate",settings["finest_step_500_to_1000_mean_rate_difference_tolerance"]),
        ("Mean z stable from 500 to 1000 units","mean_z",settings["finest_step_500_to_1000_mean_z_difference_tolerance"])]:
        difference=abs(shorter[metric]-b[metric])
        check(name,difference<=tolerance,kind="predeclared robustness criterion",difference=difference,tolerance=tolerance)
    write_json("long_run_statistics.json",groups)
    return dict(runs=len(run_rows),independent_starts=len(settings["seeds"]),duration=full,transient=settings["transient"],
                groups=groups,wall_seconds=time.perf_counter()-run_start,
                note="The guard limits and comparison tolerances are engineering criteria, not proofs or universal accuracy guarantees.")


def parameter_scan():
    settings=PLAN["parameter_scan"]
    initials,rhos,labels=[],[],[]
    local_rows=[]
    for rho in settings["rhos"]:
        if rho>1:
            a=np.sqrt((8/3)*(rho-1))
            eq=np.array([a,a,rho-1])
            eigen=np.linalg.eigvals(jacobian(eq,rho))
            largest=eigen.real.max()
            local_rows.append(dict(rho=rho,largest_real_eigenvalue=largest,
                                   linear_efold_time=1/abs(largest),local_branch_stable=largest<0,
                                   five_efold_time=5/abs(largest)))
            near=[eq+np.array([.001,0,0]),eq*np.array([-1,-1,1])+np.array([-.001,0,0])]
        else:
            near=[np.array([.001,0,0]),np.array([-.001,0,0])]
        random=np.random.default_rng(831).uniform([-15,-15,0],[15,15,45],size=(4,3))
        for index,y in enumerate([*near,*random]):
            initials.append(y); rhos.append(rho); labels.append("near_equilibrium" if index<2 else "broad_start")
    initial=np.asarray(initials)
    rho_array=np.asarray(rhos)
    results=[]
    for h in settings["steps"]:
        y=initial.copy()
        stride=round(settings["sample_interval"]/h)
        sampled=[y.copy()]
        max_state=float(np.max(abs(y)))
        for i in range(round(settings["duration"]/h)):
            y=advance(lambda a:lorenz_batch(a,rho_array),y,h)
            max_state=max(max_state,float(np.max(abs(y))))
            if (i+1)%stride==0: sampled.append(y.copy())
        states=np.asarray(sampled)
        sample_times=np.arange(len(states))*settings["sample_interval"]
        np.savez_compressed(DATA/f"parameter_scan_h{h}.npz",states=states,time=sample_times,rhos=rho_array,
                            initial=initial,labels=np.array(labels),step=h)
        check(f"Parameter scan remains finite h={h}",np.all(np.isfinite(states)) and max_state<300,
              kind="finite-window guard",maximum_absolute_state=max_state,trajectories=len(initial))
        terminal=states[sample_times>=settings["duration"]-settings["terminal_window"]]
        for column,rho in enumerate(rho_array):
            equilibria=[np.zeros(3)]
            if rho>1:
                a=np.sqrt((8/3)*(rho-1))
                equilibria.extend([np.array([a,a,rho-1]),np.array([-a,-a,rho-1])])
            eq=min(equilibria,key=lambda q:np.linalg.norm(states[-1,column]-q))
            distance=np.linalg.norm(terminal[:,column]-eq,axis=1)
            maximum=float(distance.max())
            results.append(dict(rho=rho,step=h,start_index=column%settings["starts_per_rho"],start_type=labels[column],
                                terminal_max_distance_to_equilibrium=maximum,
                                near_equilibrium_in_sampled_terminal_window=maximum<=settings["equilibrium_distance_tolerance"],
                                terminal_mean_z=terminal[:,column,2].mean()))
    csv_out("parameter_scan.csv",results)
    csv_out("local_relaxation_times.csv",local_rows)
    return dict(trajectories=len(results),settings=settings,local_times=local_rows,
                classification_warning="Not near an equilibrium within 200 units does not imply chaos or absence of an attracting equilibrium.")


def transport_exact(theta,t):
    s=PLAN["transport"]
    out=np.full_like(theta,s["forcing_mean"]/s["leak"])
    out+=(s["initial_uniform_density"]-s["forcing_mean"]/s["leak"])*np.exp(-s["leak"]*t)
    for k,amplitude in s["forcing_cosine_amplitudes"].items():
        k=int(k)
        rate=s["leak"]+1j*k*s["omega"]
        coefficient=(amplitude/2)/rate*(1-np.exp(-rate*t))
        out+=2*np.real(coefficient*np.exp(1j*k*theta))
    return out


def transport_rhs(n):
    s=PLAN["transport"]
    theta=2*np.pi*np.arange(n)/n
    force=np.full(n,s["forcing_mean"],dtype=float)
    for k,amp in s["forcing_cosine_amplitudes"].items(): force+=amp*np.cos(int(k)*theta)
    wave=np.fft.fftfreq(n,d=1/n)
    def rhs(y):
        derivative=np.fft.ifft(1j*wave*np.fft.fft(y)).real
        return force-s["leak"]*y-s["omega"]*derivative
    return theta,rhs


def transport_stress():
    s=PLAN["transport"]
    rows=[]
    for n in s["grids"]:
        theta,rhs=transport_rhs(n)
        exact=transport_exact(theta,s["duration"])
        for h in s["steps"]:
            t,state=integrate(rhs,np.full(n,s["initial_uniform_density"]),s["duration"],h)
            error=float(np.max(abs(state[-1]-exact)))
            mass_error=float(np.max(abs(2*np.pi*state.mean(axis=1)-120*np.pi)))
            rows.append(dict(grid=n,step=h,maximum_endpoint_error=error,rms_endpoint_error=np.sqrt(np.mean((state[-1]-exact)**2)),
                             maximum_mass_error=mass_error,minimum_density=np.min(state)))
            # Endpoints suffice here; exact solution is available for all times.
            np.savez_compressed(DATA/f"transport_n{n}_h{h}.npz",theta=theta,numerical=state[-1],exact=exact,
                                time=s["duration"],step=h,grid=n)
        finest=rows[-1]
        if n<64:
            check(f"Aliasing detected on underresolved grid n={n}",finest["maximum_endpoint_error"]>s["underresolved_error_minimum"],
                  kind="intentional negative control",error=finest["maximum_endpoint_error"],mass_error=finest["maximum_mass_error"])
        else:
            own=rows[-len(s["steps"]):]
            order=np.log2(own[-2]["maximum_endpoint_error"]/own[-1]["maximum_endpoint_error"])
            check(f"Resolved transport accuracy n={n}",finest["maximum_endpoint_error"]<s["fine_resolved_profile_maximum_error"]
                  and s["resolved_time_order_range"][0]<order<s["resolved_time_order_range"][1],
                  error=finest["maximum_endpoint_error"],observed_time_order=order,mass_error=finest["maximum_mass_error"])
    csv_out("transport_convergence.csv",rows)
    # Isolate the stability limit of a Fourier mode: y'=-i*376*y.
    amplification=[]
    for h in (.002,.004,.008):
        z=-1j*376*h
        factor=1+z+z*z/2+z**3/6+z**4/24
        value=np.array([1.+0j])
        for _ in range(round(1/h)): value=advance(lambda a:-1j*376*a,value,h)
        observed=float(abs(value[0]))
        predicted=float(abs(factor)**round(1/h))
        amplification.append(dict(step=h,steps=round(1/h),per_step_amplitude_factor=abs(factor),
                                  predicted_final_amplitude=predicted,observed_final_amplitude=observed,exact_amplitude=1.))
        check(f"High-frequency RK4 amplification h={h}",np.isclose(observed,predicted,rtol=1e-10)
              and ((observed>1e3) if h==.008 else observed<=1.),kind="intentional numerical stress",
              observed_amplitude=observed,exact_amplitude=1.,per_step_factor=abs(factor))
    csv_out("rk4_frequency_stress.csv",amplification)
    return dict(rows=rows,frequency_stress=amplification,
                note="Conservation alone can pass while spatial aliases or time-integration damping corrupt a profile.")


def independent_reference():
    s=PLAN["independent_reference"]
    initial=np.array([1.,1.,1.])
    t,y=integrate(lorenz,initial,s["duration"],s["step"])
    ref=solve_ivp(lambda t,y:lorenz(y),(0,s["duration"]),initial,method=s["method"],t_eval=t,rtol=s["rtol"],atol=s["atol"])
    error=float(np.max(np.linalg.norm(y-ref.y.T,axis=1)))
    check("Extended short-time independent solver comparison",ref.success and error<s["maximum_state_error"],
          maximum_state_error=error,tolerance=s["maximum_state_error"])
    np.savez_compressed(DATA/"independent_reference.npz",time=t,rk4=y,dop853=ref.y.T)
    return dict(maximum_state_error=error,settings=s)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--resume",action="store_true",help="Reuse matching completed long-run checkpoints")
    args=parser.parse_args()
    started=time.perf_counter()
    key,sources=fingerprint()
    write_json("run_provenance.json",dict(fingerprint=key,source_sha256=sources,python=sys.version,
                                         numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),plan=PLAN))
    # Cheap diagnostics first, followed by the requested long ensemble.
    print("RUN independent_reference",flush=True)
    reference=independent_reference()
    print("RUN transport_stress",flush=True)
    transport=transport_stress()
    print("RUN parameter_scan",flush=True)
    scan=parameter_scan()
    print("RUN long_runs",flush=True)
    long=long_runs(args.resume)
    result=dict(elapsed_seconds=time.perf_counter()-started,checks_passed=sum(c["passed"] for c in CHECKS),
                checks_total=len(CHECKS),failed_checks=[c for c in CHECKS if not c["passed"]],
                independent_reference=reference,transport=transport,parameter_scan=scan,long_run=long,
                completed=True,meaning="Completed bounded pilot; failed criteria remain failures. No new empirical measurements.")
    write_json("summary.json",result)
    print(json.dumps({k:result[k] for k in ("elapsed_seconds","checks_passed","checks_total","failed_checks")},default=plain),flush=True)
    if result["failed_checks"]: sys.exit(1)


if __name__=="__main__":main()
