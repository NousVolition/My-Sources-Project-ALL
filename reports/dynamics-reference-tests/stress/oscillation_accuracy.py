"""Quantify amplitude and accumulated phase errors, including unitary methods."""
from pathlib import Path
import cmath
import csv
import hashlib
import json
import math
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"
FIG=ROOT/"figures"
PLAN=json.loads((ROOT/"oscillation_accuracy_plan.json").read_text())


def factor(method,z):
    if method=="rk4": return 1+z+z*z/2+z**3/6+z**4/24
    if method=="implicit_midpoint": return (1+z/2)/(1-z/2)
    if method=="gauss_legendre4": return (1+z/2+z*z/12)/(1-z/2+z*z/12)
    if method=="exact_rotation": return cmath.exp(z)
    raise ValueError(method)


def direct_step(method,y,h,omega):
    rate=-1j*omega
    if method=="rk4":
        a=rate*y;b=rate*(y+h*a/2);c=rate*(y+h*b/2);d=rate*(y+h*c)
        return y+h*(a+2*b+2*c+d)/6
    if method=="implicit_midpoint":
        # Solve y_next = y + h*rate*(y+y_next)/2.
        return (y+h*rate*y/2)/(1-h*rate/2)
    if method=="gauss_legendre4":
        s=math.sqrt(3)/6
        a=np.array([[.25,.25-s],[.25+s,.25]])
        stages=np.linalg.solve(np.eye(2)-h*rate*a,np.full(2,rate*y))
        return y+h*(stages[0]+stages[1])/2
    if method=="exact_rotation": return y*cmath.exp(rate*h)
    raise ValueError(method)


def measure(method,h,duration,omega=376.):
    n=round(duration/h)
    if not math.isclose(n*h,duration,abs_tol=1e-12): raise ValueError("fractional step count")
    r=factor(method,-1j*omega*h)
    samples=2*math.pi/(omega*h)
    # Structural identities avoid magnifying floating-point norm error across
    # millions of hypothetical steps. Direct stepping separately measures roundoff.
    if method in ("implicit_midpoint","gauss_legendre4","exact_rotation"):
        log_amplitude=0.
    else:
        # For pure imaginary z, |R|^2=1-w^6/72+w^8/576. log1p
        # preserves the small attenuation lost by rounding abs(R) to one.
        w=omega*h
        log_amplitude=.5*n*math.log1p(-w**6/72+w**8/576)
    amplitude_error=abs(math.expm1(log_amplitude)) if log_amplitude<600 else None
    phase_resolved=samples>=PLAN["minimum_samples_per_exact_period"]
    if method=="exact_rotation":
        phase_error=0. # analytic oracle; accumulated roundoff assessed in direct audit
    elif phase_resolved:
        phase_error=n*cmath.phase(r)+omega*duration
    else: phase_error=None # too few samples; don't mistake wrapped phase for accuracy
    phase_degrees=None if phase_error is None else abs(phase_error)*180/math.pi
    stable=abs(r)<=1+2e-15
    accepted=stable and phase_resolved and amplitude_error is not None and amplitude_error<=PLAN["maximum_relative_amplitude_error"] and phase_degrees is not None and phase_degrees<=PLAN["maximum_absolute_phase_error_degrees"]
    return dict(method=method,step=h,duration=duration,step_count=n,samples_per_exact_period=samples,
                numerically_stable=stable,log10_final_amplitude=log_amplitude/math.log(10),
                relative_amplitude_error=amplitude_error,accumulated_phase_error_degrees=phase_degrees,
                phase_sampling_adequate=phase_resolved,accuracy_accepted=accepted)


def csv_out(name,rows):
    with (DATA/name).open("w",newline="",encoding="utf-8") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def run():
    started=time.perf_counter()
    DATA.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
    rows=[measure(method,h,duration,PLAN["omega"]) for method in PLAN["methods"] for duration in PLAN["durations"] for h in PLAN["steps"]]
    csv_out("oscillation_accuracy.csv",rows)
    choices=[]
    for method in PLAN["methods"]:
        for duration in PLAN["durations"]:
            accepted=[r for r in rows if r["method"]==method and r["duration"]==duration and r["accuracy_accepted"]]
            choice=max(accepted,key=lambda x:x["step"]) if accepted else None
            choices.append(dict(method=method,duration=duration,largest_tested_accepted_step=None if choice is None else choice["step"],
                                required_steps=None if choice is None else choice["step_count"],
                                amplitude_error=None if choice is None else choice["relative_amplitude_error"],
                                phase_error_degrees=None if choice is None else choice["accumulated_phase_error_degrees"]))
    csv_out("oscillation_accepted_steps.csv",choices)
    audits=[];trace=[]
    h=PLAN["direct_audit_step"];duration=PLAN["direct_audit_duration"];omega=PLAN["omega"]
    for method in PLAN["methods"]:
        y=1.+0j;angle=0.
        for i in range(round(duration/h)):
            next_y=direct_step(method,y,h,omega)
            angle+=cmath.phase(next_y/y)
            y=next_y
            if (i+1)%40==0:
                t=(i+1)*h
                trace.append(dict(method=method,time=t,real=y.real,imag=y.imag,amplitude=abs(y),
                                  accumulated_phase_error_radians=angle+omega*t))
        predicted=measure(method,h,duration,omega)
        observed_amp=abs(abs(y)-1)
        observed_phase=abs(angle+omega*duration)*180/math.pi
        passed=abs(observed_amp-predicted["relative_amplitude_error"])<1e-9 and abs(observed_phase-predicted["accumulated_phase_error_degrees"])<1e-5
        audits.append(dict(method=method,direct_steps=round(duration/h),amplitude_error=observed_amp,
                           phase_error_degrees=observed_phase,formula_amplitude_error=predicted["relative_amplitude_error"],
                           formula_phase_error_degrees=predicted["accumulated_phase_error_degrees"],formula_agreement_passed=passed))
    csv_out("oscillation_direct_traces.csv",trace)
    csv_out("oscillation_direct_audit.csv",audits)
    # Plot only points whose phase is resolved. The missing coarse points are rejected.
    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    for method,label in (("rk4","RK4"),("implicit_midpoint","Implicit midpoint"),("gauss_legendre4","Gauss-Legendre order 4")):
        selected=[r for r in rows if r["method"]==method and r["duration"]==100 and r["phase_sampling_adequate"]]
        x=[r["step"] for r in selected]
        # Zero norm error displayed as a labeled numerical floor for visibility.
        amp=[max(r["relative_amplitude_error"],1e-14) for r in selected]
        axes[0].loglog(x,amp,"o-",label=label)
        axes[1].loglog(x,[r["accumulated_phase_error_degrees"] for r in selected],"o-",label=label)
    axes[0].axhline(.001,color=".3",ls="--",label="0.1% amplitude limit")
    axes[1].axhline(1,color=".3",ls="--",label="1 degree phase limit")
    axes[0].set(xlabel="Time step",ylabel="Relative amplitude error",title="Amplitude after 100 model units")
    axes[1].set(xlabel="Time step",ylabel="Accumulated phase error (degrees)",title="Preserving amplitude can still miss the phase")
    for ax in axes: ax.legend(fontsize=7)
    fig.tight_layout()
    for extension in ("png","svg"):fig.savefig(FIG/f"06_amplitude_and_phase.{extension}",dpi=160,bbox_inches="tight")
    plt.close(fig)
    provenance={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/"oscillation_accuracy.py",ROOT/"oscillation_accuracy_plan.json"]}
    summary=dict(elapsed_seconds=time.perf_counter()-started,combinations=len(rows),accepted_combinations=sum(r["accuracy_accepted"] for r in rows),
                 choices=choices,direct_audits=audits,all_direct_audits_passed=all(r["formula_agreement_passed"] for r in audits),plan=PLAN,
                 source_sha256=provenance,
                 implication="Add amplitude and accumulated-phase gates. A conservative method still needs a phase check. No automatic solver replacement in the Lorenz or coupled-waterwheel runs is justified by this isolated linear-mode example.")
    (DATA/"oscillation-summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    if not summary["all_direct_audits_passed"]: raise RuntimeError("Direct audit disagrees with discrete-map analysis")


if __name__=="__main__":run()
