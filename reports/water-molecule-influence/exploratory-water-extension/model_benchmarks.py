"""Mathematical mechanism tests. None of these variables are water molecules.

These explicitly chosen models are positive/negative reference cases for the
measurements. They do not infer water physics from a matching-looking plot.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;D=HERE/'data'


def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()),encoding='utf-8')


def solve(f,y0,end,step=.1,rtol=1e-9):
    t=np.linspace(0,end,round(end/step)+1)
    s=solve_ivp(f,(0,end),y0,t_eval=t,method='DOP853',rtol=rtol,atol=rtol*1e-2,max_step=step*4)
    assert s.success and np.isfinite(s.y).all(),s.message
    return s.t,s.y.T


def rk4(f,x,dt):
    a=f(x);b=f(x+dt*a/2);c=f(x+dt*b/2);d=f(x+dt*c)
    return x+dt*(a+2*b+2*c+d)/6


def main():
    D.mkdir(exist_ok=True);res={};arrays={};checks={}
    # Adler/firefly/overdamped-pendulum normal form.
    rows=[]
    for mu in [0.,.5,.95,1.,1.01,1.05,1.2,2.]:
        t,z=solve(lambda t,z:mu-np.sin(z),[.3],500.,.1)
        late=t>=250;rate=float(np.polyfit(t[late],z[late,0],1)[0]);pred=np.sqrt(mu*mu-1) if mu>1 else 0.
        crossings=np.flatnonzero(np.diff(np.floor(z[:,0]/(2*np.pi)))>0)
        periods=np.diff(t[crossings]);measured=float(2*np.pi/np.mean(periods[-5:])) if len(periods)>=5 else None
        rows.append(dict(mu=mu,late_regression_rate=rate,analytic_mean_rate=pred,cycle_rate=measured,phase_excursion=float(np.ptp(z[late,0]))))
        arrays[f'adler_{mu:g}']=np.column_stack([t,z])
        if mu>1 and measured is not None:assert abs(measured-pred)<.01
        if mu<1:assert abs(z[-1,0]-np.arcsin(mu))<1e-6
    res['adler']=rows
    # Josephson RCSJ, dimensionless beta*phi''+phi'+sin(phi)=I.
    rows=[]
    for beta in [.2,2.,10.]:
        for current in [.4,.6,.8,1.2]:
            for initial in ['rest','running']:
                y0=[np.arcsin(min(current,1.)),0.] if initial=='rest' else [0.,current]
                f=lambda t,z,b=beta,I=current:np.array([z[1],(I-z[1]-np.sin(z[0]))/b])
                t,z=solve(f,y0,400.,.1);late=t>=300
                rate=float(np.mean(z[late,1]));rows.append(dict(beta=beta,current=current,initial=initial,mean_phase_rate=rate,voltage_proportional_to_phase_rate=True))
                if beta==10:arrays[f'josephson_{current:g}_{initial}']=np.column_stack([t,z])
                if beta==.2 and current==1.2 and initial=='running':
                    _,ref=solve(f,y0,400.,.1,1e-11);checks['josephson_tolerance_max_difference']=float(np.max(abs(z-ref)));assert checks['josephson_tolerance_max_difference']<1e-6
    res['josephson']=rows
    # Protein feedback: alpha=.1, K=delta=1. n=1 versus cooperative n=2.
    rows=[]
    for n in [1,2]:
        for beta in np.linspace(.5,3.5,61):
            f=lambda p:.1+beta*p**n/(1+p**n)-p
            grid=np.linspace(0,5,1001);roots=[];candidates=[]
            candidates.extend(float(q) for q in grid if abs(f(q))<1e-12)
            for left,right in zip(grid[:-1],grid[1:]):
                if f(left)*f(right)<0:
                    candidates.append(brentq(f,left,right))
            for root in sorted(candidates):
                if roots and abs(root-roots[-1]['p'])<1e-8:continue
                deriv=beta*n*root**(n-1)/(1+root**n)**2-1
                roots.append(dict(p=root,slope=deriv,stable=deriv<0))
            rows.append(dict(n=n,beta=beta,roots=roots))
    assert all(len(row['roots'])==1 for row in rows if row['n']==1)
    assert any(len(row['roots'])==3 for row in rows if row['n']==2)
    res['protein_feedback']=rows
    for p0 in [.01,4.]:
        t,z=solve(lambda t,z:.1+2.3*z*z/(1+z*z)-z,[p0],100.,.1)
        arrays[f'protein_start_{p0:g}']=np.column_stack([t,z])
    # Cusp normal form (not a reconstruction of the cropped insect equations).
    cusp=[]
    for k in np.linspace(-1,2,61):
        for r in np.linspace(-1.5,1.5,61):
            roots=np.roots([1.,0.,-k,-r]);real=sorted(float(q.real) for q in roots if abs(q.imag)<1e-7)
            stable=sum(k-3*q*q< -1e-8 for q in real)
            cusp.append([k,r,len(real),stable])
    arrays['cusp_stability']=np.array(cusp);res['cusp']=dict(equation='x_dot=r+k*x-x^3',grid_points=len(cusp),two_stable_points=sum(q[3]==2 for q in cusp))
    # Threshold delay: exact rest vs a small seed after a parameter quench.
    delays=[]
    for epsilon in [.05,.1,.2,.5,1.]:
        a0=1e-5;end=400.;t,z=solve(lambda t,z:epsilon*z-z**3,[a0],end,.05)
        threshold=.5*np.sqrt(epsilon);hit=np.flatnonzero(z[:,0]>=threshold);bracket=[float(t[hit[0]-1]),float(t[hit[0]])]
        exact=np.log((epsilon/a0**2-1)/3)/(2*epsilon)
        assert bracket[0]<=exact<=bracket[1]
        delays.append(dict(epsilon=epsilon,threshold=threshold,first_passage_bracket=bracket,analytic_time=exact))
        arrays[f'quench_{epsilon:g}']=np.column_stack([t,z])
    _,zero=solve(lambda t,z:z-z**3,[0.],10.);assert np.max(abs(zero))==0
    res['quench_delay']=delays
    # Attraction without Lyapunov stability on a circle. The distance is chord
    # distance to theta=0 modulo 2*pi, not the unwrapped angular coordinate.
    attraction=[]
    for theta0 in [.02,.05,.1]:
        t,z=solve(lambda t,z:1-np.cos(z),[theta0],8/theta0,.1,1e-11)
        exact=2*np.arctan2(np.ones_like(t),1/np.tan(theta0/2)-t)
        error=float(np.max(abs(z[:,0]-exact)));assert error<1e-7
        distance=abs(np.exp(1j*z[:,0])-1);assert distance.max()>1.999
        attraction.append(dict(initial_angle=theta0,initial_distance=float(distance[0]),maximum_distance=float(distance.max()),final_distance=float(distance[-1]),time_to_opposite_point=float(1/np.tan(theta0/2)),analytic_max_error=error))
        arrays[f'attraction_{theta0:g}']=np.column_stack([t,z,distance])
    res['attraction_vs_stability']=dict(equation='theta_dot=1-cos(theta) on the unit circle',exact_solution='cot(theta(t)/2)=cot(theta0/2)-t',rows=attraction,
        conclusion='All points approach theta=0 modulo 2*pi, but arbitrarily small positive initial angles eventually reach the opposite point. Attraction therefore does not imply Lyapunov stability.')
    # Maxwell--Bloch laser equation supplied by the user, and its adiabatic limit.
    laser=[]
    for lam in [-.5,.1,1.]:
        for gamma in [20.,100.]:
            f=lambda t,z,lam=lam,gamma=gamma:np.array([z[1]-z[0],gamma*(z[0]*z[2]-z[1]),gamma*(lam+1-z[2]-lam*z[0]*z[1])])
            t,z=solve(f,[.02,.02*(1+lam),1+lam],80.,.1)
            _,reduced=solve(lambda t,z,lam=lam:lam*z*(1-z*z)/(1+lam*z*z),[.02],80.,.1)
            laser.append(dict(pump_parameter=lam,fast_decay=gamma,final_E=float(z[-1,0]),reduction_RMS_error=float(np.sqrt(np.mean((z[:,0]-reduced[:,0])**2)))))
            arrays[f'laser_{lam:g}_{gamma:g}']=np.column_stack([t,z,reduced])
    res['laser']=laser
    # Standard SIR is explicitly specified because the screenshot omits its definitions.
    epidemic=[]
    for beta in [.5,1.,2.]:
        def f(t,z):
            s,i,r=z;return [-beta*s*i,beta*s*i-i,i]
        t,z=solve(f,[.999,.001,0.],100.,.05)
        err=float(np.max(abs(z.sum(axis=1)-1)));assert err<1e-10 and z.min()>-1e-12
        peak=int(z[:,1].argmax());epidemic.append(dict(beta=beta,gamma=1.,initial_reproduction_ratio=beta*.999,peak_time=float(t[peak]),peak_infected=float(z[peak,1]),conservation_error=err))
        arrays[f'sir_{beta:g}']=np.column_stack([t,z])
    res['epidemic_threshold']=epidemic
    # Identical phase oscillators on equal-degree rings versus global coupling.
    oscillator=[]
    for topology in ['ring','all']:
        w=np.ones((9,9))-np.eye(9) if topology=='all' else np.roll(np.eye(9),1,axis=0)+np.roll(np.eye(9),-1,axis=0)
        w/=w.sum(axis=1,keepdims=True)
        for coupling in [-1.,0.,.5,2.]:
            for seed in range(3):
                theta=np.random.default_rng(413+seed).uniform(-np.pi,np.pi,9)
                def f(t,z):
                    diff=z[None,:]-z[:,None];return 1+coupling*np.sum(w*np.sin(diff),axis=1)
                t,z=solve(f,theta,80.,.1)
                order=np.abs(np.exp(1j*z).mean(axis=1));r2=np.abs(np.exp(2j*z).mean(axis=1))
                oscillator.append(dict(topology=topology,coupling=coupling,seed=seed,late_order=float(order[-201:].mean()),late_second_harmonic_order=float(r2[-201:].mean())))
                if seed==0:arrays[f'phase_{topology}_{coupling:g}']=np.column_stack([t,z,order,r2])
                if coupling==0:assert np.ptp(order)<1e-10
    res['identical_phase_networks']=oscillator
    # Three saddles: exact local eigenvalues and noisy winnerless competition.
    rho=np.eye(3)+1.6*np.roll(np.eye(3),-1,axis=1)+.5*np.roll(np.eye(3),1,axis=1)
    hetero=[]
    for eta in [0.,1e-4,.01]:
        for seed in range(3):
            dt=.01;end=300.;rng=np.random.default_rng(910+seed);signal=abs(rng.normal(size=(round(end/.1),3)))
            x=np.array([.7,.15,.05])+rng.uniform(0,.01,3);states=[];ts=[];noise_index=-1
            for step in range(round(end/dt)+1):
                if step%10==0:states.append(x.copy());ts.append(step*dt)
                if step==round(end/dt):break
                noise=eta*signal[min(int(round(step*dt,8)/.1+1e-7),len(signal)-1)]
                x=rk4(lambda q:q*(1-rho@q)+noise,x,dt)
                assert x.min()>=-1e-12 and np.isfinite(x).all()
            z=np.array(states);winner=z.argmax(axis=1);changes=np.flatnonzero(np.diff(winner)!=0)+1;dwell=np.diff(np.r_[0,changes,len(winner)])*.1
            seq=winner[np.r_[0,changes]]
            hetero.append(dict(noise=eta,seed=seed,switches=len(changes),median_dwell=float(np.median(dwell)),early_median_dwell=float(np.median(dwell[:max(1,len(dwell)//2)])),late_median_dwell=float(np.median(dwell[max(1,len(dwell)//2):])) if len(dwell)>1 else None,winner_sequence=seq.tolist()))
            if seed==0:arrays[f'heteroclinic_{eta:g}']=np.column_stack([ts,z])
    res['heteroclinic']=dict(rho=rho,axis_eigenvalues=[-1.,-.6,.5],rows=hetero,noise_definition='Rectified Gaussian rate held for a fixed 0.1 time units; eta=0 is deterministic. No absolute Brownian increments.')
    # Nine-state modular competition: hierarchy is deliberately encoded in rho.
    matrix=np.eye(9)
    for i in range(9):
        for j in range(9):
            gi,ki=divmod(i,3);gj,kj=divmod(j,3)
            if i==j:continue
            matrix[i,j]=(1.6 if (kj-ki)%3==1 else .5) if gi==gj else (1.06 if (gj-gi)%3==1 else .97)
    hierarchy=[]
    for eta in [0.,1e-4,.001]:
        rng=np.random.default_rng(770);dt=.01;end=600.;noise=abs(rng.normal(size=(round(end/.1),9)))
        x=rng.uniform(.01,.1,9);states=[];ts=[]
        for step in range(round(end/dt)+1):
            if step%20==0:states.append(x.copy());ts.append(step*dt)
            if step==round(end/dt):break
            drive=eta*noise[min(int(round(step*dt,8)/.1+1e-7),len(noise)-1)]
            x=rk4(lambda q:q*(1-matrix@q)+drive,x,dt);assert x.min()>=-1e-12
        z=np.array(states);chunks=z.reshape(-1,3,3).sum(axis=2);sw=int(np.count_nonzero(np.diff(z.argmax(axis=1))));csw=int(np.count_nonzero(np.diff(chunks.argmax(axis=1))))
        share=chunks/np.maximum(chunks.sum(axis=1,keepdims=True),1e-20)
        hierarchy.append(dict(noise=eta,activity_switches=sw,chunk_switches=csw,ratio=sw/max(csw,1),max_chunk_share=float(np.max(share)),late_fraction_with_one_chunk_above_80pct=float(np.mean(share[len(share)//2:].max(axis=1)>.8))))
        arrays[f'hierarchy_{eta:g}']=np.column_stack([ts,z,chunks])
    res['engineered_hierarchy']=dict(rho=matrix,rows=hierarchy,warning='Three groups of three and their coupling are imposed in this model; this is not discovered water chunking.')
    # Diffusive versus directed chains; no pacemaker versus explicitly assigned leader.
    spatial=[]
    for mode,delta in [('none',0.),('diffusive',.2),('pacemaker',.2),('pacemaker',1.)]:
        rng=np.random.default_rng(661);x=rng.uniform(.01,.5,(16,3));dt=.01;end=200.;states=[];ts=[]
        def f(q):
            growth=q*(1-q@rho.T);coupling=np.zeros_like(q)
            if mode=='diffusive':
                coupling[1:]+=q[:-1]-q[1:];coupling[:-1]+=q[1:]-q[:-1]
            elif mode=='pacemaker':coupling[1:]=q[:-1]-q[1:]
            return growth+delta*coupling
        for step in range(round(end/dt)+1):
            if step%20==0:states.append(x.copy());ts.append(step*dt)
            if step==round(end/dt):break
            x=rk4(f,x,dt);assert x.min()>=-1e-12
        z=np.array(states);late=z[len(z)//2:];err=float(np.sqrt(np.mean((late[:,1:]-late[:,:1])**2)))
        spatial.append(dict(mode=mode,coupling=delta,late_RMS_difference_from_unit0=err));arrays[f'chain_{mode}_{delta:g}']=np.concatenate([np.array(ts)[:,None],z.reshape(len(z),-1)],axis=1)
    res['spatial_activity_networks']=spatial
    # Same physical regularized noise path with RK4 step refinement.
    eta=.001;end=10.;rng=np.random.default_rng(77);noise=abs(rng.normal(size=(100,3)))
    ends=[]
    for dt in [.01,.005]:
        x=np.array([.7,.15,.05])
        for step in range(round(end/dt)):
            drive=eta*noise[min(int(round(step*dt,8)/.1+1e-7),99)];x=rk4(lambda q:q*(1-rho@q)+drive,x,dt)
        ends.append(x)
    checks['regularized_noise_step_refinement_error']=float(np.max(abs(ends[1]-ends[0])));assert checks['regularized_noise_step_refinement_error']<1e-7
    res['checks']=dict(**checks,passed=True)
    res['scope']='Mechanism reference models, not evidence that water obeys these equations. Prespecified illustrative parameters; negative and weak findings retained.'
    np.savez_compressed(D/'model_trajectories.npz',**arrays);save(D/'model_results.json',res)
    print(json.dumps(dict(model_families=len([k for k in res if k not in ['checks','scope']]),arrays=len(arrays),checks=res['checks'])),flush=True)


if __name__=='__main__':main()
