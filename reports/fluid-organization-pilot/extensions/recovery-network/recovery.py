"""3D periodic Navier--Stokes forcing-shutoff experiment, in model units."""
import argparse
import hashlib
import time
import numpy as np
from core import PilotFlow, initial, energy, L
from common import ROOT,P,save,arrivals,interval,prediction_comparison

def strength(t,arm,ramp):
    if arm=='unchanged':return 1.
    if arm=='abrupt':return 0.
    if arm=='gradual':return max(0.,1.-t/ramp)
    raise ValueError(arm)

def rk4(f,h,t,dt,force,arm,ramp):
    ks=[];rates=[];speed=0.
    for j,a in enumerate((0.,.5,.5,1.)):
        z=h if j==0 else h+a*dt*ks[-1]
        u=f.real(z);forcing=strength(t+a*dt,arm,ramp)*force
        ks.append(f.rhs_only(z,u)+forcing)
        rates.append(np.array([f.inner(z,forcing),f.nu*f.inner(f.curl(z),f.curl(z))])/L**3)
        speed=max(speed,float(np.sqrt((u*u).sum(0)).max()))
    w=(1,2,2,1)
    return f.project(h+dt*sum(a*k for a,k in zip(w,ks))/6),dt*sum(a*b for a,b in zip(w,rates))/6,speed

def observe(f,h,force):
    e=energy(f,h);w=f.curl(h);ens=.5*f.inner(w,w)/L**3
    u=f.real(h);speed=np.sqrt((u*u).sum(0))
    corr=f.inner(h,force)/max(np.sqrt(f.inner(h,h)*f.inner(force,force)),1e-30)
    helicity=f.inner(h,w)/L**3
    low=f.k2<=(2*np.pi/L)**2*2.01
    lowfrac=f.inner(h*low,h)/max(f.inner(h,h),1e-30)
    rate=f.rhs_only(h)+force
    rmsrate=np.sqrt(f.inner(rate,rate)/L**3)
    return np.array([e,ens,helicity,corr,lowfrac,speed.max(),np.mean(speed**4),rmsrate])

def fingerprint():
    names=('recovery.py','core.py','vendor/numerics.py','protocol.json')
    return {n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}

def run(seed,n=None,dt=None,label='base'):
    c=P['fluid'];n=n or c['n'];dt=dt or c['dt'];dest=ROOT/'data'/'fluid';dest.mkdir(parents=True,exist_ok=True)
    path=dest/f'{label}_{seed}.npz';meta=path.with_suffix('.json');hashes=fingerprint()
    if path.exists():
        import json
        old=json.loads(meta.read_text(encoding='utf-8'))
        if old['sources']!=hashes or old['n']!=n or old['dt']!=dt:raise ValueError('Stale '+str(path))
        if old['sha256']!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('Corrupt cached data')
        return path
    start=time.perf_counter();f=PilotFlow(n,c['nu']);force=initial(f,7788,kind='beltrami')*c['forcing_rms']
    if f.nu*f.k2[f.keep].max()*dt>.6:raise ValueError('Diffusion stability gate')
    speed=np.random.default_rng(seed+8000).uniform(*c['initial_rms_range'])
    h=initial(f,seed,speed);stride=round(c['save_dt']/dt)
    history=[];bt=[];budget=np.zeros(2);burn_e0=energy(f,h);maxcfl=0.
    for j in range(round(c['burn']/dt)+1):
        if j%stride==0:history.append(observe(f,h,force));bt.append(j*dt-c['burn'])
        if j<round(c['burn']/dt):
            h,b,v=rk4(f,h,j*dt,dt,force,'unchanged',c['ramp_time']);budget+=b;maxcfl=max(maxcfl,v*dt/f.dx)
    burn_res=abs(energy(f,h)-burn_e0-budget[0]+budget[1])/burn_e0
    h0=h.copy();e0=energy(f,h0);hs=[h0.copy() for _ in c['arms']]
    budget=np.zeros((3,2));times=[];obs=[];energies=[];diff=[];budgets=[];divmax=0.
    clone_error=max(float(np.max(abs(q-h0))) for q in hs)
    for j in range(round(c['end']/dt)+1):
        t=j*dt
        if j%stride==0:
            times.append(t);energies.append([energy(f,q) for q in hs])
            obs.append([observe(f,q,force*strength(t,arm,c['ramp_time'])) for q,arm in zip(hs,c['arms'])])
            diff.append([energy(f,q-hs[0]) for q in hs]);budgets.append(budget.copy())
            divmax=max(divmax,max(float(np.max(abs(f.real(1j*sum(k*v for k,v in zip(f.k,q)))))) for q in hs))
        if j<round(c['end']/dt):
            for a,arm in enumerate(c['arms']):
                hs[a],b,v=rk4(f,hs[a],t,dt,force,arm,c['ramp_time']);budget[a]+=b;maxcfl=max(maxcfl,v*dt/f.dx)
    es=np.array(energies);bs=np.array(budgets);res=float(np.max(abs(es-e0-bs[:,:,0]+bs[:,:,1]))/e0)
    np.savez_compressed(path,time=times,energy=es,difference_energy=diff,budget=bs,
                        observables=obs,history_times=bt,history=history,initial_h=h0,final_h=np.stack(hs))
    save(meta,{'seed':seed,'n':n,'dt':dt,'sources':hashes,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
               'initial_rms':speed,'burn_budget_residual':float(burn_res),'energy_budget_residual':res,
               'divergence_max':divmax,'clone_error':clone_error,'max_cfl':maxcfl,'seconds':time.perf_counter()-start})
    print(f'fluid {label} {seed}: {time.perf_counter()-start:.1f}s; balance {res:.2g}',flush=True)
    return path

def analyze():
    import json,csv
    c=P['fluid'];rows=[];parts={};curves=[];meta=[];t=None
    for split in ('train','validation','test'):
        data={'current':[],'history':[],'y':[],'groups':[],'slot':[]}
        for seed in c[split+'_seeds']:
            path=ROOT/'data'/'fluid'/f'base_{seed}.npz';z=np.load(path);t=z['time'];e=z['energy'];e0=e[0,0]
            h=z['history'];ht=z['history_times'];b=z['budget'];meta.append(json.loads(path.with_suffix('.json').read_text()))
            curve=e/e0;curves.append(curve);ix=t<=8+1e-9
            integrals=np.trapezoid(curve[ix],t[ix],axis=0)
            for a,arm in enumerate(c['arms']):
                first,settle=arrivals(t,curve[:,a],c['settling_energy_fraction'],c['residence_time'])
                row={'seed':seed,'split':split,'arm':arm,'first_arrival':first,'sustained_arrival':settle,
                     'right_censored':settle is None,'integrated_energy_0_8':float(integrals[a]),
                     'integrated_difference_0_8':float(integrals[a]-integrals[0]),
                     'peak_above_initial':float(max(0,curve[:,a].max()-1)),
                     'injected_work_over_initial_energy':float(b[-1,a,0]/e0),
                     'dissipated_over_initial_energy':float(b[-1,a,1]/e0)}
                rows.append(row)
                if a>0:
                    data['current'].append(np.r_[h[-1],float(a==1),float(a==2)])
                    data['history'].append(np.concatenate([(h[-1]-h[np.argmin(abs(ht+lag))])/lag for lag in c['history_lags']]))
                    data['y'].append([integrals[a],integrals[a]-integrals[0]]);data['groups'].append(seed);data['slot'].append(a)
        parts[split]={k:np.array(v) for k,v in data.items()}
    output=ROOT/'results';output.mkdir(exist_ok=True)
    with (output/'recovery_runs.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={}
    for arm in c['arms']:
        rr=[r for r in rows if r['arm']==arm]
        valid=[r['sustained_arrival'] for r in rr if r['sustained_arrival'] is not None]
        summary[arm]={'n':len(rr),'settled':len(valid),'censored':len(rr)-len(valid),
                      'settling_time_among_observed':interval(valid) if valid else None,
                      'integrated_energy_0_8':interval([r['integrated_energy_0_8'] for r in rr]),
                      'peak_above_initial':interval([r['peak_above_initial'] for r in rr])}
    paired_delays=[b['sustained_arrival']-a['sustained_arrival'] for a,b in zip(rows[1::3],rows[2::3])
                   if a['sustained_arrival'] is not None and b['sustained_arrival'] is not None]
    summary['gradual_minus_abrupt_settling']=interval(paired_delays) if paired_delays else None
    predictions,preds=prediction_comparison(parts,c['prediction_targets'])
    save(output/'recovery_prediction.json',predictions)
    np.savez_compressed(output/'recovery_predictions.npz',truth=parts['test']['y'],groups=parts['test']['groups'],**preds)
    sensitivity={}
    for frac in c['settling_sensitivity_fractions']:
        for residence in c['residence_sensitivity']:
            sensitivity[f'fraction{frac}_residence{residence}']={}
            for a,arm in enumerate(c['arms']):
                vals=[arrivals(t,q[:,a],frac,residence)[1] for q in curves];v=[x for x in vals if x is not None]
                sensitivity[f'fraction{frac}_residence{residence}'][arm]={'settled':len(v),'censored':len(vals)-len(v),'mean_observed':float(np.mean(v)) if v else None}
    refinement=[]
    for seed in c['refinement_seeds']:
        base=np.load(ROOT/'data'/'fluid'/f'base_{seed}.npz')
        half=np.load(ROOT/'data'/'fluid'/f'half_{seed}.npz')
        for label,refpath,reference in [('time',f'half_{seed}.npz',base),('grid',f'fine_{seed}.npz',half)]:
            ref=np.load(ROOT/'data'/'fluid'/refpath);e0=reference['energy'][0,0]
            curve_error=float(np.max(abs(ref['energy']-reference['energy']))/e0)
            delta=[]
            for a in (1,2):
                ts=[arrivals(z['time'],z['energy'][:,a]/z['energy'][0,0],.1,1)[1] for z in (reference,ref)]
                delta.append(abs(ts[0]-ts[1]) if None not in ts else None)
            refinement.append({'seed':seed,'comparison':label,'max_energy_difference_over_initial':curve_error,'settling_differences':delta,
                               'pass':curve_error<c['gates']['relative_energy_curve'] and all(x is not None and x<=c['gates']['settling_time_absolute'] for x in delta)})
    summary.update({'validation':{'max_energy_budget_residual':max(m['energy_budget_residual'] for m in meta),
                                  'max_burn_budget_residual':max(m['burn_budget_residual'] for m in meta),
                                  'max_divergence':max(m['divergence_max'] for m in meta),'max_cfl':max(m['max_cfl'] for m in meta),
                                  'clone_error_max':max(m['clone_error'] for m in meta),'refinement':refinement},
                    'settling_sensitivity':sensitivity})
    save(output/'recovery_summary.json',summary)
    np.savez_compressed(output/'recovery_curves.npz',time=t,energy=np.array(curves))
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--screen',action='store_true');ap.add_argument('--analyze',action='store_true');args=ap.parse_args()
    if args.analyze:analyze()
    elif args.screen:run(999,label='screen')
    else:
        c=P['fluid']
        for seed in c['train_seeds']+c['validation_seeds']+c['test_seeds']:run(seed)
        for seed in c['refinement_seeds']:
            run(seed,dt=c['fine_dt'],label='half');run(seed,n=c['fine_n'],dt=c['fine_dt'],label='fine')
        analyze()
