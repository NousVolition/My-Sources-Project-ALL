"""Run the complete deterministic numerical study and write auditable outputs."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.stats import t as student_t
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator, ScalarFormatter

from core import (integrate, advance, logistic_exact, lorenz, jacobian,
                  variational, lyapunov, delayed_branch, wheel_rhs, wheel_project)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
checks, summary = [], {}
plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 160, "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    raise TypeError(type(value))


def save_json(name, value):
    (DATA / name).write_text(json.dumps(value, indent=2, default=plain, allow_nan=False) + "\n", encoding="utf-8")


def csv_out(name, rows):
    with (DATA / name).open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def check(name, condition, **values):
    checks.append(dict(name=name, passed=bool(condition), **values))
    print(("PASS " if condition else "FAIL ") + name, flush=True)


def figure(name, fig):
    fig.tight_layout()
    for extension in ("png", "svg"):
        fig.savefig(FIG / (name + "." + extension), bbox_inches="tight")
    plt.close(fig)


def methods():
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))
    for problem in ("logistic", "RC"):
        for method, order in (("euler", 1), ("heun", 2), ("rk4", 4)):
            errors = []
            for h in (.2, .1, .05, .025):
                if problem == "logistic":
                    fun, initial = lambda x: x*(1-x), .1
                    exact = lambda t: logistic_exact(t, initial)
                else:
                    fun, initial = lambda x: (2-x)/.75, 0.
                    exact = lambda t: 2*(1-np.exp(-t/.75))
                t, state = integrate(fun, [initial], 4, h, method)
                error = np.max(np.abs(state[:, 0]-exact(t)))
                errors.append(error)
                rows.append(dict(problem=problem, method=method, step=h, max_error=error))
            observed = np.log2(errors[-2]/errors[-1])
            check(f"{problem}: {method} convergence", abs(observed-order) < .22,
                  observed_order=observed, expected_order=order, finest_max_error=errors[-1])
            if problem == "logistic":
                axes[0].loglog([.2, .1, .05, .025], errors, "o-", label=f"{method}: p={observed:.2f}")
    csv_out("method_convergence.csv", rows)
    t = np.linspace(0, 8, 500)
    slope_t, slope_x = np.meshgrid(np.linspace(0, 8, 17), np.linspace(0, 2, 17))
    slope = slope_x*(1-slope_x)
    length = np.sqrt(1+slope*slope)
    axes[1].quiver(slope_t, slope_x, 1/length, slope/length, color=".8", angles="xy")
    trajectories = []
    for x0 in (.05, .2, .5, 1., 1.5, 2.):
        exact = logistic_exact(t, x0)
        axes[1].plot(t, exact, label=f"x0={x0}")
        trajectories.extend(dict(time=a, initial=x0, x=b) for a, b in zip(t, exact))
    csv_out("logistic_exact_trajectories.csv", trajectories)
    for h in (.25, 1.5, 2.8):
        t, y = integrate(lambda x: x*(1-x), [.1], 80*h, h, "euler")
        axes[2].plot(t, y[:, 0], ".-", markersize=3, label=f"Euler h={h}")
        csv_out(f"logistic_euler_h{h}.csv", [dict(time=a, x=b) for a,b in zip(t, y[:,0])])
        if h == 2.8:
            oscillations = np.count_nonzero(np.diff(np.sign(y[-30:, 0]-1)))
            check("Coarse Euler logistic instability is detected", oscillations > 10,
                  crossings_last_30_steps=int(oscillations), equilibrium_multiplier=1-h)
    axes[0].set(xlabel="Step size", ylabel="Maximum absolute error", title="Known orders recovered")
    axes[0].set_xticks([.025,.05,.1,.2])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].xaxis.set_major_formatter(ScalarFormatter())
    axes[1].set(xlabel="Time", ylabel="x", title="Logistic ODE: smooth approach to 1")
    axes[2].set(xlabel="Time", ylabel="x", xlim=(0,40), title="Large steps invent oscillation")
    for ax in axes: ax.legend(fontsize=8)
    figure("01_numerical_methods", fig)
    decay_rows = []
    for h in (.5, 1.5, 2.2):
        t, y = integrate(lambda x: -x, [1.], 20*h, h, "euler")
        decay_rows.extend(dict(step=h, time=a, numerical=b, exact=np.exp(-a)) for a,b in zip(t,y[:,0]))
        check(f"Euler linear stability h={h}", np.allclose(y[:,0], (1-h)**np.arange(21), atol=1e-10),
              final=y[-1,0], growth_factor=1-h)
    csv_out("euler_stability_control.csv", decay_rows)
    summary["RC"] = dict(equilibrium=2., time_constant=.75,
                         time_to_95_percent=.75*np.log(20), units="illustrative consistent units")


def scalar_examples():
    kinds = {"stable": (lambda x: -x**3, lambda t,x: x/np.sqrt(1+2*x*x*t)),
             "unstable": (lambda x: x**3, lambda t,x: x/np.sqrt(1-2*x*x*t)),
             "half_stable": (lambda x: x*x, lambda t,x: x/(1-x*t)),
             "neutral": (lambda x: x*0, lambda t,x: np.full_like(t,x))}
    fig, axes = plt.subplots(2, 2, figsize=(9,6))
    rows = []
    for ax, (name, (fun, exact)) in zip(axes.flat, kinds.items()):
        max_error = 0.
        for initial in (-1., -.3, .3, 1.):
            t, y = integrate(fun, [initial], .4, .0005)
            max_error = max(max_error, float(np.max(np.abs(y[:,0]-exact(t,initial)))))
            ax.plot(t, y[:,0], label=f"x0={initial}")
            rows.extend(dict(kind=name, initial=initial, time=a, numerical=b, exact=c)
                        for a,b,c in zip(t,y[:,0],exact(t,initial)))
        check(f"Phase line: {name}", max_error < 2e-10, max_absolute_error=max_error)
        ax.set(title=name.replace("_", " "), xlabel="Time", ylabel="x")
        ax.axhline(0, color=".6", lw=.7)
    figure("02_scalar_stability", fig)
    csv_out("scalar_stability.csv", rows)
    t = np.linspace(0, 1.8, 361)
    inverse_rows, inverse_error = [], 0.
    for initial in (-1., .4, 1., 1.1):
        tn, yn = integrate(lambda x: x*(x-1), [initial], 1.8, .005)
        exact = initial/(initial+(1-initial)*np.exp(tn))
        inverse_error = max(inverse_error, np.max(abs(yn[:,0]-exact)))
        inverse_rows.extend(dict(initial=initial,time=a,numerical=b,exact=c) for a,b,c in zip(tn,yn[:,0],exact))
    csv_out("inverse_phase_portrait.csv", inverse_rows)
    check("Inverse portrait: candidate x'=x(x-1)", inverse_error < 1e-8,
          max_absolute_error=inverse_error, note="One consistent equation; the sketch does not uniquely identify f.")


def nonuniqueness():
    times = np.linspace(0, 3, 601)
    rows, residual = [], 0.
    fig, axes = plt.subplots(1, 2, figsize=(10.5,3.8))
    for delay in (0., 1., 2.):
        for sign in (-1., 1.):
            x, dx = delayed_branch(times, delay, sign)
            residual = max(residual, np.max(abs(dx-np.cbrt(x))))
            rows.extend(dict(time=t,delay=delay,sign=sign,x=a,derivative=b,residual=b-np.cbrt(a))
                        for t,a,b in zip(times,x,dx))
            axes[0].plot(times, x, label=f"delay {delay:g}, sign {sign:g}")
    zero = solve_ivp(lambda t,y: np.cbrt(y), (0,3), [0.], t_eval=times, rtol=1e-10, atol=1e-12)
    axes[0].plot(times, zero.y[0], "k--", lw=2, label="solver from exact zero")
    check("Nonunique exact branches satisfy the ODE", residual < 5e-15,
          max_residual=residual, common_initial_value=0, distinct_branches=7)
    check("Zero-initialized solver stays on zero branch", np.max(abs(zero.y)) == 0.,
          implication="Finding one solution does not establish uniqueness.")
    errors = []
    for epsilon in (1e-3, 1e-6, 1e-9):
        t, y = integrate(lambda x: np.cbrt(x), [epsilon], 1., .0001)
        exact = (epsilon**(2/3)+2*t/3)**1.5
        errors.append(dict(epsilon=epsilon, fixed_step=.0001, final_numerical=y[-1,0],
                           final_exact=exact[-1], final_error=abs(y[-1,0]-exact[-1])))
        axes[1].plot(t,y[:,0], label=f"RK4 x0={epsilon:g}")
    axes[1].plot(times[:201], (2*times[:201]/3)**1.5, "k--", label="exact immediate departure")
    axes[0].set(title="Same initial value, different exact solutions", xlabel="Time",ylabel="x")
    axes[1].set(title="Tiny perturbations select a branch",xlabel="Time",ylabel="x")
    for ax in axes: ax.legend(fontsize=7)
    figure("03_nonuniqueness", fig)
    csv_out("nonunique_branches.csv", rows)
    csv_out("nonunique_perturbations.csv", errors)
    summary["nonuniqueness"] = dict(max_residual=residual, branches_checked=7,
                                    warning="The vector field is not locally Lipschitz at zero. This is not chaotic sensitivity.")


def overdamped():
    times = np.linspace(0, 12, 1201)
    fig, axes = plt.subplots(1, 2, figsize=(10.5,3.8))
    rows, errors = [], []
    for mass in (1., .1, .01, .001):
        matrix = np.array([[0.,1.],[-1/mass,-1/mass]])
        def rhs(t,y):
            x,v,dissipated = y
            return [v,(-v-x)/mass,v*v]
        sol = solve_ivp(rhs, (0,12), [1.,0.,0.], method="Radau", t_eval=times, rtol=2e-10,atol=2e-12)
        exact = np.array([expm(matrix*t) @ [1.,0.] for t in times])
        error = np.max(abs(sol.y[:2].T-exact))
        energy = .5*sol.y[0]**2+.5*mass*sol.y[1]**2
        budget = np.max(abs(energy+sol.y[2]-.5))
        reduced_error = np.max(abs(sol.y[0]-np.exp(-times)))
        crossings = int(np.count_nonzero(np.diff(np.sign(sol.y[0]))))
        check(f"Damped oscillator m={mass}", sol.success and error < 3e-8 and budget < 3e-9,
              exact_state_error=error, energy_budget_error=budget, reduced_model_max_error=reduced_error,
              zero_crossings=crossings)
        errors.append(dict(mass=mass, model_error=reduced_error, numerical_error=error, energy_budget_error=budget,crossings=crossings))
        axes[0].plot(times, sol.y[0], label=f"mass={mass:g}")
        axes[1].plot(times, energy, label=f"mass={mass:g}")
        rows.extend(dict(mass=mass,time=t,x=x,v=v,energy=e,dissipated=d,reduced=np.exp(-t))
                    for t,x,v,e,d in zip(times,sol.y[0],sol.y[1],energy,sol.y[2]))
    check("Overdamped approximation improves as inertia falls",
          errors[1]["model_error"] > errors[2]["model_error"] > errors[3]["model_error"])
    axes[0].plot(times,np.exp(-times),"k--",label="zero-inertia x=exp(-t)")
    axes[0].set(xlabel="Time",ylabel="Displacement",title="Fixed drag=1, spring=1; reduce mass")
    axes[1].set(xlabel="Time",ylabel="Mechanical energy",title="Energy lost equals integrated drag loss")
    for ax in axes: ax.legend(fontsize=8)
    figure("04_overdamped_limit",fig)
    csv_out("overdamped_trajectories.csv",rows)
    csv_out("overdamped_errors.csv",errors)
    summary["overdamped"] = errors


def saddle_and_reaction():
    fig, axes = plt.subplots(1, 2, figsize=(10.5,3.8))
    rows = []
    for mu in np.logspace(-6,-1,51):
        # expm1 avoids cancellation in r-x-exp(-x) near the saddle-node.
        f = lambda x: mu-x-np.expm1(-x)
        roots = [brentq(f,-3,0,xtol=1e-14),brentq(f,0,3,xtol=1e-14)]
        for root in roots:
            approx = np.sign(root)*np.sqrt(2*mu)
            rows.append(dict(mu=mu,root=root,normal_form_root=approx,residual=f(root),
                             eigenvalue=np.expm1(-root),relative_root_error=abs(root-approx)/abs(root)))
    residual = max(abs(r["residual"]) for r in rows)
    near_error = max(r["relative_root_error"] for r in rows if r["mu"] <= 1.000001e-6)
    check("Saddle-node roots and local approximation", residual < 2e-13 and near_error < .0005,
          max_root_residual=residual, relative_normal_form_error_at_mu_1e_6=near_error)
    for sign,style,label in ((-1,"--","unstable"),(1,"-","stable")):
        branch = [r for r in rows if np.sign(r["root"])==sign]
        axes[0].plot([r["mu"] for r in branch],[r["root"] for r in branch],style,label=label)
    axes[0].set(xlabel="mu = r - 1",ylabel="Equilibrium x", title="r - x - exp(-x): two roots for r>1")
    axes[0].legend()
    csv_out("saddle_node_roots.csv",rows)
    reaction_rows = []
    # Two explicitly assumed mass-action scenarios, since the screenshot is cropped.
    kf, kr, total, initial = 1., .5, 1., .05
    for scenario,rate,capacity in (("closed",kf*total,kf*total/(kf+kr)),
                                   ("fixed_A",kf*.8,kf*.8/kr)):
        if scenario == "closed":
            def rhs(y):
                a,x = y
                reaction = kf*a*x-kr*x*x
                return np.array([-reaction,reaction])
            t, y = integrate(rhs,[total-initial,initial],10.,.01)
            concentration = y[:,1]
            mass_error = np.max(abs(np.sum(y,axis=1)-total))
        else:
            t,y = integrate(lambda x: rate*x-kr*x*x,[initial],10.,.01)
            concentration, mass_error = y[:,0], None
        exact = logistic_exact(t,initial,rate,capacity)
        error = np.max(abs(concentration-exact))
        check(f"Autocatalysis under stated assumption: {scenario}",
              error < 2e-10 and (mass_error is None or mass_error < 1e-12),
              exact_error=error, conserved_A_plus_X_error=mass_error, positive_equilibrium=capacity)
        axes[1].plot(t,concentration,label=scenario)
        reaction_rows.extend(dict(scenario=scenario,time=a,x=b,exact=c) for a,b,c in zip(t,concentration,exact))
    axes[1].set(xlabel="Time",ylabel="X concentration",title="Autocatalysis depends on the A supply")
    axes[1].legend()
    figure("05_bifurcation_and_reaction",fig)
    csv_out("autocatalysis.csv",reaction_rows)
    summary["saddle_node"] = dict(critical_r=1., normal_form="x_dot = (r-1) - x^2/2 + O(x^3)",
                                 relative_root_error_at_mu_1e_6=near_error)


def lorenz_checks():
    sigma,beta = 10.,8/3
    hopf = sigma*(sigma+beta+3)/(sigma-beta-1)
    rows = []
    for rho in (.5,1.,1.1,10.,24.,hopf-.001,hopf,hopf+.001,28.):
        eqs = [("origin",np.zeros(3))]
        if rho > 1:
            for sign in (-1,1):
                a = sign*np.sqrt(beta*(rho-1))
                eqs.append((f"C{sign:+d}",np.array([a,a,rho-1])))
        for name,eq in eqs:
            vals = np.linalg.eigvals(jacobian(eq,rho))
            for ev in vals:
                rows.append(dict(rho=rho,equilibrium=name,real=ev.real,imag=ev.imag,
                                 residual=np.max(abs(lorenz(eq,rho)))))
    csv_out("lorenz_equilibrium_eigenvalues.csv",rows)
    at = [r["real"] for r in rows if r["rho"]==hopf and r["equilibrium"]=="C+1"]
    check("Lorenz local Hopf threshold", abs(max(at)) < 2e-13,
          rho_H=hopf, largest_real_eigenvalue=max(at))
    before = max(r["real"] for r in rows if r["rho"]==hopf-.001 and r["equilibrium"]=="C+1")
    after = max(r["real"] for r in rows if r["rho"]==hopf+.001 and r["equilibrium"]=="C+1")
    check("Lorenz equilibrium stability changes across Hopf", before < 0 < after,
          largest_real_before=before,largest_real_after=after)
    rng = np.random.default_rng(631)
    reflection = np.array([-1.,-1.,1.])
    symmetry_error = max(np.max(abs(lorenz(y*reflection)-lorenz(y)*reflection)) for y in rng.normal(size=(100,3)))
    check("Lorenz reflection symmetry", symmetry_error == 0., max_residual=symmetry_error)
    initial = np.array([1.,1.,1.])
    reference = solve_ivp(lambda t,y:lorenz(y),(0,2),initial,method="DOP853",rtol=2e-13,atol=2e-14,dense_output=True)
    convergence = []
    for h in (.01,.005,.0025,.00125):
        t,y = integrate(lorenz,initial,2.,h)
        error = np.max(np.linalg.norm(y-reference.sol(t).T,axis=1))
        convergence.append(dict(step=h,max_norm_error=error))
    order = np.log2(convergence[-2]["max_norm_error"]/convergence[-1]["max_norm_error"])
    check("Lorenz short-time trajectory convergence", 3.8 < order < 4.2, observed_order=order,
          finest_max_error=convergence[-1]["max_norm_error"])
    csv_out("lorenz_time_convergence.csv",convergence)
    t,v = integrate(variational,np.r_[initial,np.eye(3).ravel()],.5,.00025)
    logvolume = np.array([np.linalg.slogdet(mat.reshape(3,3))[1] for mat in v[:,3:]])
    exact = -(sigma+1+beta)*t
    vol_error = np.max(abs(logvolume-exact))
    check("Phase volume follows exact divergence identity", vol_error < 1e-7,
          max_log_volume_error=vol_error, exact_contraction_rate=-(sigma+1+beta))
    csv_out("lorenz_phase_volume.csv",[dict(time=a,log_volume=b,exact=c) for a,b,c in zip(t,logvolume,exact)])
    fig,axes = plt.subplots(1,2,figsize=(10.5,3.8))
    axes[0].plot(t,logvolume,label="integrated tangent matrix")
    axes[0].plot(t,exact,"k--",label="exact divergence")
    axes[0].set(xlabel="Time",ylabel="log(V/V0)",title="State-space volume contracts")
    t,pair = integrate(lambda y: np.r_[lorenz(y[:3]),lorenz(y[3:])],np.r_[initial,initial+[1e-8,0,0]],40.,.0025)
    separation=np.linalg.norm(pair[:,:3]-pair[:,3:],axis=1)
    axes[1].semilogy(t,separation)
    axes[1].set(xlabel="Time",ylabel="Distance between states",title="Nearby trajectories can separate")
    axes[0].legend(fontsize=8)
    figure("06_lorenz_volume_and_separation",fig)
    np.savez_compressed(DATA/"lorenz_nearby_pair.npz",time=t,states=pair,separation=separation)
    summary["lorenz_local"] = dict(hopf=hopf,phase_volume_rate=-(sigma+1+beta),
                                   phase_volume_log_error=vol_error,time_convergence=convergence)


def lorenz_ensemble():
    rows,block_rows = [],[]
    seeds=(101,202,303,404)
    for h in (.005,.0025):
        for seed in seeds:
            initial=np.random.default_rng(seed).uniform([-.5,-.5,20.],[.5,.5,30.])
            spectrum,blocks,states=lyapunov(initial,h=h)
            rows.append(dict(seed=seed,step=h,initial_x=initial[0],initial_y=initial[1],initial_z=initial[2],
                             lambda1=spectrum[0],lambda2=spectrum[1],lambda3=spectrum[2],
                             sum=spectrum.sum(),mean_z=states[:,2].mean()))
            block_rows.extend(dict(seed=seed,step=h,block=i,lambda1=a,lambda2=b,lambda3=c)
                              for i,(a,b,c) in enumerate(blocks))
            np.savez_compressed(DATA/f"lorenz_seed{seed}_h{h}.npz",states=states,block_rates=blocks,
                                sample_interval=.1,transient=30.,duration=180.,initial=initial)
            check(f"Lorenz Lyapunov run seed={seed}, h={h}", spectrum[0]>.5 and abs(spectrum.sum()+41/3)<.001,
                  spectrum=spectrum, sum_error=spectrum.sum()+41/3)
    csv_out("lorenz_ensemble.csv",rows)
    csv_out("lorenz_lyapunov_blocks.csv",block_rows)
    estimates=[]
    for h in (.005,.0025):
        rates=np.array([r["lambda1"] for r in rows if r["step"]==h])
        half=student_t.ppf(.975,len(rates)-1)*rates.std(ddof=1)/np.sqrt(len(rates))
        estimates.append(dict(step=h,mean_largest_rate=rates.mean(),sd_between_initial_conditions=rates.std(ddof=1),
                              approximate_95_t_interval=[rates.mean()-half,rates.mean()+half],min=rates.min(),max=rates.max()))
    check("Lorenz largest growth rate survives step refinement",
          abs(estimates[0]["mean_largest_rate"]-estimates[1]["mean_largest_rate"])<.1,
          difference=abs(estimates[0]["mean_largest_rate"]-estimates[1]["mean_largest_rate"]))
    summary["lorenz_ensemble"]=dict(transient=30.,measurement_duration=180.,seeds=seeds,qr_interval=.1,
                                    estimates=estimates,interval_caveat="Approximate t intervals across four independent initial conditions; not a proof, not physical uncertainty, and not a bound on finite-time bias. Ten-time-unit blocks are not treated as independent samples.")
    fig,axes=plt.subplots(1,2,figsize=(10.5,3.8))
    for h,shift in ((.005,-.07),(.0025,.07)):
        rates=[r["lambda1"] for r in rows if r["step"]==h]
        axes[0].plot(np.arange(4)+shift,rates,"o",label=f"h={h}")
    axes[0].set(xticks=np.arange(4),xticklabels=seeds,xlabel="Independent initial-condition seed",ylabel="Largest finite-time Lyapunov rate",title="Positive growth persists across runs")
    axes[0].legend(fontsize=8)
    item=np.load(DATA/"lorenz_nearby_pair.npz")
    orbit=item["states"][item["time"]>=10,:3]
    axes[1].plot(orbit[:,0],orbit[:,2],lw=.35)
    axes[1].set(xlabel="x",ylabel="z",title="Classic Lorenz model\nsigma=10, beta=8/3, rho=28")
    figure("07_lorenz_ensemble",fig)


def waterwheel():
    rows=[]
    for n in (16,32,64):
        theta,rhs=wheel_rhs(n)
        initial=np.r_[45+np.sin(theta)+27*np.cos(theta)+2*np.cos(3*theta),1.]
        for h in (.002,.001,.0005):
            t,y=integrate(rhs,initial,2.,h)
            projected=wheel_project(y,theta)
            # Independent high-order reference to the b=1 Lorenz reduction.
            ref=solve_ivp(lambda t,z:lorenz(z,beta=1.),(0,2),[1.,1.,1.],method="DOP853",t_eval=t,rtol=2e-13,atol=2e-14)
            projected_error=np.max(np.linalg.norm(projected-ref.y.T,axis=1))
            mass=2*np.pi*np.mean(y[:,:-1],axis=1)
            exact_mass=2*np.pi*(60+(45-60)*np.exp(-t))
            budget_error=np.max(abs(mass-exact_mass))
            a3=2*np.mean(y[:,:-1]*np.sin(3*theta),axis=1)
            b3=2*np.mean(y[:,:-1]*np.cos(3*theta),axis=1)
            higher_error=np.max(abs(np.hypot(a3,b3)-2*np.exp(-t)))
            rows.append(dict(grid=n,step=h,lorenz_projection_error=projected_error,mass_error=budget_error,
                             third_mode_amplitude_error=higher_error,minimum_mass_density=np.min(y[:,:-1])))
            check(f"Waterwheel budget and positivity n={n}, h={h}",budget_error<1e-8 and np.min(y[:,:-1])>0,
                  mass_budget_error=budget_error,minimum_mass_density=np.min(y[:,:-1]))
            np.savez_compressed(DATA/f"waterwheel_n{n}_h{h}.npz",time=t,theta=theta,state=y,projected=projected,
                                reference_lorenz=ref.y.T,mass=mass,exact_mass=exact_mass,third_amplitude=np.hypot(a3,b3))
    fine=[r for r in rows if r["grid"]==64]
    order=np.log2(fine[-2]["lorenz_projection_error"]/fine[-1]["lorenz_projection_error"])
    check("Waterwheel first-mode mapping converges to Lorenz beta=1", 3.7<order<4.3,
          observed_order=order,finest_projection_error=fine[-1]["lorenz_projection_error"])
    amplitude_errors=[r["third_mode_amplitude_error"] for r in fine]
    check("Unforced higher waterwheel mode decays at leakage rate",
          amplitude_errors[-1]<amplitude_errors[-2]<amplitude_errors[0] and amplitude_errors[-1]<1e-5,
          errors=amplitude_errors)
    # A bandwidth-limited solution: adequate grids should agree up to floating-point error.
    finalstates=[np.load(DATA/f"waterwheel_n{n}_h0.0005.npz")["projected"] for n in (16,32,64)]
    grid_error=max(np.max(abs(a-finalstates[-1])) for a in finalstates[:-1])
    check("Waterwheel resolved Fourier grids agree",grid_error<1e-9,max_grid_difference=grid_error)
    csv_out("waterwheel_convergence.csv",rows)
    item=np.load(DATA/"waterwheel_n64_h0.0005.npz")
    fig,axes=plt.subplots(1,3,figsize=(13.5,3.7))
    axes[0].plot(item["time"],item["mass"],label="computed")
    axes[0].plot(item["time"],item["exact_mass"],"k--",label="exact source/leak budget")
    axes[0].set(xlabel="Time",ylabel="Total mass",title="Waterwheel mass balance")
    axes[0].legend(fontsize=8)
    axes[1].plot(item["time"],item["projected"][:,0])
    axes[1].set(xlabel="Time",ylabel="Angular velocity",title="Rotation coupled to water distribution")
    axes[2].plot(item["time"],item["third_amplitude"],label="computed third mode")
    axes[2].plot(item["time"],2*np.exp(-item["time"]),"k--",label="2 exp(-t)")
    axes[2].set(xlabel="Time",ylabel="Harmonic amplitude",title="Higher mode decays; no torque feedback")
    axes[2].legend(fontsize=8)
    figure("08_waterwheel_budget",fig)
    summary["waterwheel"]=dict(beta=1.,sigma=10.,rho=28.,source_mean=60.,source_cosine=28.,leak=1.,drag=10.,torque=10.,
                               time_order=order,grid_difference=grid_error,finest=fine[-1],
                               limitation="Only Fourier modes 0,1,3 are populated. Grid agreement demonstrates resolution of this smooth case, not arbitrary nonsmooth inflow.")


def main():
    started=time.perf_counter()
    for task in (methods,scalar_examples,nonuniqueness,overdamped,saddle_and_reaction,lorenz_checks,lorenz_ensemble,waterwheel):
        print("RUN " + task.__name__,flush=True)
        task()
    save_json("checks.json",checks)
    summary["checks"]=dict(passed=sum(c["passed"] for c in checks),total=len(checks))
    save_json("summary.json",summary)
    environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,
                     matplotlib=matplotlib.__version__,elapsed_seconds=time.perf_counter()-started,
                     source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("*.py")})
    save_json("environment.json",environment)
    print(json.dumps(summary["checks"]),flush=True)
    if not all(c["passed"] for c in checks): sys.exit(1)


if __name__ == "__main__": main()
