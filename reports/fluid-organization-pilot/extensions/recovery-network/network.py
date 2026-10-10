"""Specified positive-orthant GLV cycle; frozen offline identification controls.

x_i' = x_i (1 - sum_j C_ij x_j), C_ii=1,
C_i,i+1=1.6*c, C_i,i+2=0.6*c (indices modulo 3).
Log-state integration avoids artificial clipping near a saddle. No noise added.
"""
import argparse,csv,hashlib,json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import nnls
from scipy.linalg import expm
from common import ROOT,P,save,interval,prediction_comparison

def matrix(coupling=1.,a=None,b=None,variant='cycle'):
    p=P['network'];a=p['a'] if a is None else a;b=p['b'] if b is None else b
    C=np.eye(3)
    for i in range(3):C[i,(i+1)%3]=a*coupling;C[i,(i+2)%3]=b*coupling
    if variant=='edge_removed':C[0,1]=0.
    if variant=='row_reversed':C[0,1],C[0,2]=C[0,2],C[0,1]
    return C

def rhs_log(t,y,C):return 1.-C@np.exp(y)

def integrate(y0,C,end,dt=None,tight=False):
    p=P['network'];dt=dt or p['sample_dt'];times=np.linspace(0,end,round(end/dt)+1)
    sol=solve_ivp(rhs_log,(0,end),y0,args=(C,),method='DOP853',t_eval=times,
                  max_step=p['max_step']/(2 if tight else 1),
                  rtol=p['rtol']/(10 if tight else 1),atol=p['atol']/(10 if tight else 1))
    if not sol.success or not np.isfinite(sol.y).all():raise RuntimeError(sol.message)
    return times,sol.y.T

def initial_state(seed):
    rng=np.random.default_rng(seed);x=rng.uniform(.05,.8,3)
    return np.log(x)

def run(seed,coupling=1.,variant='cycle',tight=False):
    p=P['network'];C=matrix(coupling,variant=variant)
    _,pre=integrate(initial_state(seed),C,p['burn'],tight=tight);y0=pre[-1]
    kick=np.zeros(3);kick[(np.argmax(y0)+1)%3]=p['kick_log_amplitude']
    t,base=integrate(y0,C,p['end'],tight=tight)
    _,pert=integrate(y0+kick,C,p['end'],tight=tight)
    name=f'{variant}_c{coupling:g}_seed{seed}'+('_tight' if tight else '')
    path=ROOT/'data'/'network'/(name+'.npz');path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,time=t,base=base,perturbed=pert,initial_log=y0,kick=kick,C=C,prehistory=pre)
    return path

def switches(t,y):
    labels=np.argmax(y,axis=1);ix=np.flatnonzero(np.diff(labels)!=0)+1
    return labels,ix

def saddle_visits(y,radius=.2):
    x=np.exp(y);distance=np.linalg.norm(x[:,None,:]-np.eye(3)[None,:,:],axis=2)
    close=distance.min(1)<radius;labels=distance.argmin(1)[close]
    return labels[np.r_[True,np.diff(labels)!=0]] if len(labels) else np.array([],int)

def waiting(t,y,index,horizon):
    label=np.argmax(y[index]);later=np.flatnonzero(np.argmax(y[index+1:],axis=1)!=label)
    wait=float(t[index+1+later[0]]-t[index]) if len(later) else float('inf')
    return min(wait,horizon),wait>horizon

def event_rows(seed,t,y,arm,coupling,variant):
    labels,ix=switches(t,y);starts=np.r_[0,ix];ends=np.r_[ix,len(t)-1]
    return [{'seed':seed,'coupling':coupling,'variant':variant,'arm':arm,'state':int(labels[i]),
             'start':float(t[i]),'duration':float(t[j]-t[i]),'left_censored':i==0,
             'right_censored':k==len(starts)-1,'next_state':int(labels[j]) if k<len(starts)-1 else None}
            for k,(i,j) in enumerate(zip(starts,ends))]

def fit_coefficients(training,shuffle=False):
    """Finite differences from training observations only; positive competition."""
    X=[];Y=[]
    for z in training:
        y=z['base'];x=np.exp(y);der=np.gradient(y,P['network']['sample_dt'],axis=0,edge_order=2)
        # Exclude endpoints, whose one-sided derivative has a different error.
        X.append(x[1:-1]);Y.append(1-der[1:-1])
    X=np.concatenate(X);Y=np.concatenate(Y)
    if shuffle:Y=Y[np.random.default_rng(332).permutation(len(Y))]
    return np.array([nnls(X,Y[:,i])[0] for i in range(3)])

def fit_linear(training):
    """Nine-parameter linear vector field: matched parameter-count comparator."""
    X=[];Y=[]
    for z in training:
        x=np.exp(z['base']);X.append(x[1:-1]);Y.append(np.gradient(x,P['network']['sample_dt'],axis=0,edge_order=2)[1:-1])
    return np.linalg.lstsq(np.concatenate(X),np.concatenate(Y),rcond=None)[0].T

def km_survival(events,state,times):
    events=[e for e in events if e['state']==state and not e['left_censored']]
    if not events:return np.ones_like(times)
    d=np.array([e['duration'] for e in events]);observed=np.array([not e['right_censored'] for e in events]);surv=np.ones_like(times,float);s=1.
    for q in np.unique(d[observed]):
        s*=1-np.sum((d==q)&observed)/np.sum(d>=q)
        surv[times>=q]=s
    return surv

def analyze():
    p=P['network'];dest=ROOT/'data'/'network';out=ROOT/'results';out.mkdir(exist_ok=True)
    allseeds=sum([p[k+'_seeds'] for k in ('train','validation','test')],[])
    sets={k:[np.load(dest/f'cycle_c1_seed{s}.npz') for s in p[k+'_seeds']] for k in ('train','validation','test')}
    learned=fit_coefficients(sets['train']);shuffled=fit_coefficients(sets['train'],True);linear=fit_linear(sets['train'])
    save(out/'learned_network_parameters.json',{'true':matrix().tolist(),'learned':learned.tolist(),'shuffled_teacher':shuffled.tolist(),
         'fixed_prior':matrix(a=1.3,b=.8).tolist(),'linear':linear.tolist(),'max_abs_coefficient_error':float(abs(learned-matrix()).max()),
         'method':'offline nonnegative least squares on finite-difference log derivatives; frozen after training; not the cited online teacher algorithm'})
    events=[];parts={};model_errors={k:[] for k in ('oracle','fixed_prior','learned','shuffled_teacher','linear')};examples={}
    response_errors={k:[] for k in model_errors};ode_wait={'oracle_dynamics':[],'learned_dynamics':[]}
    for split,zz in sets.items():
        d={'current':[],'history':[],'y':[],'groups':[],'slot':[]};records=[]
        for seed,z in zip(p[split+'_seeds'],zz):
            t=z['time']
            for arm,key in enumerate(('base','perturbed')):
                y=z[key];ev=event_rows(seed,t,y,key,1.,'cycle');events.extend(ev)
                labels,ix=switches(t,y)
                for slot,now in enumerate(p['history_sample_times']):
                    j=round(now/p['sample_dt']);before=ix[ix<=j];start=before[-1] if len(before) else 0
                    age=t[j]-t[start];prev=t[before[-1]]-t[before[-2]] if len(before)>=2 else 0.
                    h=np.concatenate([(y[j]-y[j-round(lag/p['sample_dt'])])/lag for lag in p['history_lags']])
                    wait,censored=waiting(t,y,j,p['forecast_horizon'])
                    d['current'].append(np.r_[np.exp(y[j]),y[j],rhs_log(0,y[j],matrix()),arm])
                    d['history'].append(np.r_[h,age,prev]);d['y'].append([wait]);d['groups'].append(seed);d['slot'].append(arm*len(p['history_sample_times'])+slot)
                    records.append({'state':int(labels[j]),'age':float(age),'wait':wait,'censored':censored,'seed':seed})
                    if split=='test':
                        for name,C in [('oracle_dynamics',matrix()),('learned_dynamics',learned)]:
                            tf,yf=integrate(y[j],C,p['forecast_horizon'])
                            ode_wait[name].append(waiting(tf,yf,0,p['forecast_horizon'])[0])
            if split=='test':
                errors={k:[] for k in model_errors};forecasts={k:[] for k in model_errors}
                for arm,key in enumerate(('base','perturbed')):
                    truth=np.exp(z[key][:round(p['forecast_horizon']/p['sample_dt'])+1]);y0=z[key][0]
                    for model,C in [('oracle',matrix()),('fixed_prior',matrix(a=1.3,b=.8)),('learned',learned),('shuffled_teacher',shuffled)]:
                        tt,yy=integrate(y0,C,p['forecast_horizon']);pred=np.exp(yy)
                        errors[model].append(float(np.sqrt(np.mean((pred-truth)**2))))
                        forecasts[model].append(pred)
                        if seed==p['test_seeds'][0] and arm==1:examples[model]=pred
                    pred=np.stack([expm(linear*ttt)@np.exp(y0) for ttt in tt])
                    errors['linear'].append(float(np.sqrt(np.mean((pred-truth)**2))))
                    forecasts['linear'].append(pred)
                    if seed==p['test_seeds'][0] and arm==1:examples.update({'linear':pred,'truth':truth,'time':tt})
                for model in model_errors:model_errors[model].append(float(np.mean(errors[model])))
                truth_delta=np.exp(z['perturbed'][:len(tt)])-np.exp(z['base'][:len(tt)])
                for model,ff in forecasts.items():response_errors[model].append(float(np.sqrt(np.mean(((ff[1]-ff[0])-truth_delta)**2))))
        parts[split]={k:np.array(v) for k,v in d.items()};parts[split]['records']=records
    # All fitting finished above uses only training observations. Validation tunes ridge only.
    pred,preds=prediction_comparison(parts,['restricted_wait_to_next_dominant_state'])
    train_events=[e for e in events if e['seed'] in p['train_seeds']];test=parts['test'];markov=[];semi=[]
    for r in test['records']:
        ev=[e for e in train_events if e['state']==r['state']]
        rate=sum(not e['right_censored'] for e in ev)/sum(e['duration'] for e in ev)
        markov.append(float(-np.expm1(-rate*p['forecast_horizon'])/rate) if rate else p['forecast_horizon'])
        grid=np.linspace(0,p['forecast_horizon'],301);sv=km_survival(train_events,r['state'],grid+r['age']);s0=sv[0]
        semi.append(float(np.trapezoid(sv/s0,grid)) if s0>1e-12 else markov[-1])
    for label,values in [('label_only_markov',markov),('age_conditioned_semi_markov',semi),*ode_wait.items()]:
        err=np.array([np.sqrt(np.mean((np.array(values)[test['groups']==g]-test['y'][test['groups']==g,0])**2)) for g in np.unique(test['groups'])])
        pred['restricted_wait_to_next_dominant_state'][label]={**interval(err),'per_run_rmse':err.tolist()};preds[label]=values
    pred['censored_test_observations']=sum(r['censored'] for r in test['records'])
    pred['n_test_observations']=len(test['records'])
    # Next destination learned from training transitions; this three-cycle is deliberately simple.
    transition=np.zeros((3,3),int)
    for e in train_events:
        if e['next_state'] is not None:transition[e['state'],e['next_state']]+=1
    test_events=[e for e in events if e['seed'] in p['test_seeds'] and e['next_state'] is not None]
    pred['next_destination']={'training_counts':transition.tolist(),'label_only_accuracy':float(np.mean([transition[e['state']].argmax()==e['next_state'] for e in test_events])),
                              'n_observed_test_transitions':len(test_events),'note':'Single outgoing route; high accuracy is a topology control, not evidence for memory.'}
    save(out/'network_history_prediction.json',pred)
    np.savez_compressed(out/'network_wait_predictions.npz',truth=test['y'],groups=test['groups'],**preds)
    np.savez_compressed(out/'network_forecast_example.npz',**examples)
    forecast={model:{**interval(e),'per_run_rmse':e} for model,e in model_errors.items()}
    forecast['learned_minus_fixed']=interval(np.array(model_errors['learned'])-model_errors['fixed_prior'])
    forecast['learned_minus_linear']=interval(np.array(model_errors['learned'])-model_errors['linear'])
    forecast['paired_disturbance_response']={k:{**interval(e),'per_run_rmse':e} for k,e in response_errors.items()}
    save(out/'network_forecast.json',forecast)
    sweep=[];refine=[]
    for path in sorted(dest.glob('*.npz')):
        if '_tight' in path.stem:continue
        z=np.load(path);variant=path.stem.split('_c')[0];coupling=float(path.stem.split('_c')[1].split('_')[0]);seed=int(path.stem.split('seed')[1])
        for key in ('base','perturbed'):
            labels,ix=switches(z['time'],z[key]);ev=event_rows(seed,z['time'],z[key],key,coupling,variant)
            visits=saddle_visits(z[key])
            completed=[e['duration'] for e in ev if not e['left_censored'] and not e['right_censored']]
            sweep.append({'seed':seed,'coupling':coupling,'variant':variant,'arm':key,'transitions':len(ix),
                          'saddle_neighborhood_visits':len(visits),'transitions_between_saddle_visits':max(0,len(visits)-1),
                          'mean_complete_dominance_duration':float(np.mean(completed)) if completed else None,
                          'max_complete_dominance_duration':float(max(completed)) if completed else None,
                          'last_interval_right_censored':True,'final_log_min':float(z[key][-1].min())})
        refpath=path.with_name(path.stem+'_tight.npz')
        if refpath.exists():
            ref=np.load(refpath);err=max(float(np.sqrt(np.mean((np.exp(z[k])-np.exp(ref[k]))**2))) for k in ('base','perturbed'))
            diffs=[];counts_equal=True
            for key in ('base','perturbed'):
                _,i=switches(z['time'],z[key]);_,j=switches(ref['time'],ref[key]);counts_equal &= len(i)==len(j)
                if len(i)==len(j) and len(i):diffs.append(float(np.max(abs(z['time'][i]-ref['time'][j]))))
            terr=max(diffs,default=0.)
            refine.append({'file':path.name,'state_rms_difference':err,'switch_time_max_difference':terr,'same_event_count':counts_equal,
                           'pass':bool(counts_equal and err<p['gates']['state_rms_refinement'] and terr<=p['gates']['transition_time_absolute'])})
    for name,rr in [('network_events.csv',events),('network_regimes.csv',sweep)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
    save(out/'network_validation.json',{'refinement':refine,'saddle_eigenvalues':[-1.,1.-p['a'],1.-p['b']],
              'local_cycle_contraction_ratio':((p['a']-1)/(1-p['b']))**3,
              'observed_trajectory_files':len(sweep)//2,'invariant_coordinate_planes':True,
              'noise':0,'state_definition':'Largest component; dominance intervals are not epsilon-neighborhood dwell times.'})
    print('network analysis complete',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--analyze',action='store_true');args=ap.parse_args()
    if not args.analyze:
        p=P['network']
        for seed in sum([p[k+'_seeds'] for k in ('train','validation','test')],[]):run(seed)
        for seed in p['test_seeds'][:8]:
            for c in p['coupling_values']:
                if c!=1.:run(seed,c)
            for variant in ('edge_removed','row_reversed'):run(seed,variant=variant)
        for seed in p['test_seeds'][:4]:
            for c in (0.8,1.,1.2):run(seed,c,tight=True)
        print('network trajectories complete',flush=True)
    analyze()
