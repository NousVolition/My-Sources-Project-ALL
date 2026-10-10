"""Replica-level exploratory analysis; no independent-frame significance claims."""
from pathlib import Path
import argparse, csv, json
import numpy as np
from scipy.stats import t as student_t
from scipy.cluster.vq import kmeans2, vq
from scipy.optimize import lsq_linear
from water_battery import cases, save

HERE=Path(__file__).resolve().parent;D=HERE/'data'


def summary(values):
    a=np.asarray(values,float);a=a[np.isfinite(a)];n=len(a)
    if not n:return dict(n=0,mean=None,SD=None,nominal_95pct_t_interval=None,values=[])
    sd=float(a.std(ddof=1)) if n>1 else None
    half=float(student_t.ppf(.975,n-1)*sd/np.sqrt(n)) if n>1 else None
    return dict(n=n,mean=float(a.mean()),SD=sd,nominal_95pct_t_interval=[float(a.mean()-half),float(a.mean()+half)] if half is not None else None,values=a.tolist())


def segments(mask):
    edges=np.diff(np.r_[False,mask,False].astype(int));return list(zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)))


def phase_metrics(time,polarization,frequency,threshold=.15):
    z=polarization[:,0]+1j*polarization[:,1];amp=abs(z);valid=amp>=threshold
    groups=segments(valid);start,end=max(groups,key=lambda x:x[1]-x[0]) if groups else (0,0)
    valid_fraction=float(valid.mean());fraction=(end-start)/len(time)
    dphase=np.angle(z)-2*np.pi*frequency*time
    concentration=float(abs(np.mean(np.exp(1j*dphase[valid])))) if valid.any() else None
    slope=None;excursion=None;winding=None
    if end-start>=3:
        q=np.unwrap(np.angle(z[start:end]))-2*np.pi*frequency*time[start:end]
        slope=float(np.polyfit(time[start:end],q,1)[0]/(2*np.pi*frequency))
        excursion=float(np.ptp(q));winding=float((q[-1]-q[0])/(2*np.pi))
    tracking=bool(valid_fraction>=.95 and fraction>=.8 and concentration is not None and concentration>=.9 and slope is not None and abs(slope)<=.05 and excursion<np.pi)
    response=np.mean(z*np.exp(-2j*np.pi*frequency*time))
    return dict(amplitude_threshold=threshold,valid_fraction=valid_fraction,longest_valid_fraction=fraction,phase_concentration=concentration,
        relative_phase_slope=slope,phase_excursion=excursion,net_phase_winding=winding,finite_window_tracking=tracking,
        complex_response_real=float(response.real),complex_response_imag=float(response.imag),response_magnitude=float(abs(response)),mean_lag_radians=float(-np.angle(response)))


def sequence_metrics(t,p):
    q=p@np.exp(2j*np.pi*np.arange(3)/3);omega=2*np.pi/6
    coherence=float(abs(np.mean(q*np.exp(-1j*omega*t)))/max(np.sqrt(np.mean(abs(q)**2)),1e-15))
    strength=np.linalg.norm(p,axis=1);dominant=np.argmax(p,axis=1)
    valid=(strength>=.15)&(p.max(axis=1)>.1)
    # Sample at the drive's two-ps dwell spacing rather than count frame jitter as transitions.
    stride=max(1,int(round(2/np.median(np.diff(t)))));states=dominant[::stride];ok=valid[::stride]
    forward=sum(bool(ok[i] and ok[i+1] and (states[i+1]-states[i])%3==1) for i in range(len(states)-1))
    reverse=sum(bool(ok[i] and ok[i+1] and (states[i+1]-states[i])%3==2) for i in range(len(states)-1))
    return dict(six_ps_coherence=coherence,mean_polarization_norm=float(strength.mean()),valid_fraction=float(valid.mean()),forward_sampled_transitions=forward,reverse_sampled_transitions=reverse)


def recurrence(loaded):
    runs=[loaded.get(('water',r,'baseline')) for r in range(3)]
    if any(d is None for d in runs):return None
    # Train only on replica 0. Features are spatial cell orientations, not named molecule identities.
    features=[]
    for data in runs:
        mask=data['time_ps']>=20
        features.append(data['cell_polarization'][mask].reshape(-1,24))
    center=features[0].mean(axis=0);scale=features[0].std(axis=0);scale=np.maximum(scale,.03)
    X=[(f-center)/scale for f in features];_,s,_=np.linalg.svd(X[0]-X[0].mean(axis=0),full_matrices=False)
    pca3=float(np.sum(s[:3]**2)/np.sum(s*s));cluster=[]
    for k in [2,3,4,6]:
        centres,label=kmeans2(X[0],k,minit='++',seed=630+k,iter=100)
        labels=[label]+[vq(x,centres)[0] for x in X[1:]];lag=10
        train=np.ones((k,k))*.5
        np.add.at(train,(labels[0][:-lag],labels[0][lag:]),1)
        transition=train/train.sum(axis=1,keepdims=True)
        stationary=np.bincount(labels[0],minlength=k)+.5;stationary=stationary/stationary.sum()
        rows=[]
        for rep,lab in enumerate(labels):
            a,b=lab[:-lag],lab[lag:];pred=transition[a]
            loss=float(-np.log(pred[np.arange(len(a)),b]).mean());baseline=float(-np.log(stationary[b]).mean())
            changes=np.flatnonzero(np.diff(lab)!=0)+1;dwell=np.diff(np.r_[0,changes,len(lab)])*.1
            rows.append(dict(replica=rep,role='training' if rep==0 else 'held_out',state_counts=np.bincount(lab,minlength=k),median_dwell_ps=float(np.median(dwell)),max_dwell_ps=float(np.max(dwell)),
                 markov_log_loss=loss,stationary_log_loss=baseline,markov_accuracy=float(np.mean(pred.argmax(axis=1)==b)),persistence_accuracy=float(np.mean(a==b))))
        cluster.append(dict(k=k,transition=transition,rows=rows))
    # Predict collective orientation at 1 ps and 5 ps using a fitted linear map.
    pruns=[d['polarization'][d['time_ps']>=20] for d in runs];mean=pruns[0].mean(axis=0);forecast=[]
    for lag in [10,50]:
        train=np.column_stack([pruns[0][:-lag],np.ones(len(pruns[0])-lag)])
        coef=np.linalg.lstsq(train,pruns[0][lag:],rcond=None)[0]
        for rep in [1,2]:
            a,b=pruns[rep][:-lag],pruns[rep][lag:];pred=np.column_stack([a,np.ones(len(a))])@coef
            forecast.append(dict(replica=rep,horizon_ps=lag*.1,linear_prediction_MSE=float(np.mean((pred-b)**2)),persistence_MSE=float(np.mean((a-b)**2)),constant_training_mean_MSE=float(np.mean((mean-b)**2))))
    # Deliberately try the proposed GLV form on a declared observable mapping.
    # A_i=P_i^2; three Cartesian components are chosen coordinates, not discovered chunks.
    activities=[p*p for p in pruns];A=activities[0];derivative=(A[2:]-A[:-2])/.2;middle=A[1:-1];coeff=[]
    for i in range(3):
        design=np.column_stack([middle[:,i],-middle[:,i,None]*middle])
        coeff.append(np.linalg.solve(design.T@design+1e-6*np.eye(4),design.T@derivative[:,i]))
    coeff=np.array(coeff);glv=[]
    for rep in [1,2]:
        start=activities[rep][:-10].copy();target=activities[rep][10:];pred=start.copy();failed=np.zeros(len(pred),bool)
        def f(q):return q*(coeff[:,0]-q@coeff[:,1:].T)
        with np.errstate(over='ignore',invalid='ignore'):
            for _ in range(50):
                a=f(pred);b=f(pred+.01*a);c=f(pred+.01*b);d=f(pred+.02*c);pred+=.02*(a+2*b+2*c+d)/6
                bad=(~np.isfinite(pred).all(axis=1))|(abs(pred).max(axis=1)>100)|(pred.min(axis=1)<-1e-8)
                failed|=bad;pred[bad]=0. # Only to finish the other forecasts; failed rows are never scored as success.
        glv.append(dict(replica=rep,failed_forecasts=int(failed.sum()),total_forecasts=len(pred),GLV_MSE=None if failed.any() else float(np.mean((pred-target)**2)),persistence_MSE=float(np.mean((start-target)**2)),constant_training_mean_MSE=float(np.mean((A.mean(axis=0)-target)**2))))
    # A second, explicitly constrained attempt. Its constant source approximates
    # a positive mean input, not a resolved stochastic noise process.
    constrained=[]
    for i in range(3):
        design=np.column_stack([np.ones(len(middle)),middle[:,i],-middle[:,i,None]*middle])
        fit=lsq_linear(np.vstack([design,1e-3*np.eye(5)]),np.r_[derivative[:,i],np.zeros(5)],bounds=([0,-np.inf,0,0,0],[np.inf]*5),tol=1e-12,max_iter=500)
        assert fit.success;constrained.append(fit.x)
    constrained=np.array(constrained);constrained_forecast=[]
    for rep in [1,2]:
        start=activities[rep][:-10].copy();target=activities[rep][10:];pred=start.copy()
        def drift(q):return constrained[:,0]+q*(constrained[:,1]-q@constrained[:,2:].T)
        for _ in range(50):
            a=drift(pred);b=drift(pred+.01*a);c=drift(pred+.01*b);d=drift(pred+.02*c);pred+=.02*(a+2*b+2*c+d)/6
        good=np.isfinite(pred).all() and pred.min()>=-1e-8
        constrained_forecast.append(dict(replica=rep,valid_prediction=bool(good),MSE=float(np.mean((pred-target)**2)) if good else None,persistence_MSE=float(np.mean((start-target)**2))))
    return dict(training_replica=0,held_out_replicas=[1,2],cell_features=24,first_three_PC_variance_fraction=pca3,clusters=cluster,linear_forecasts=forecast,
        GLV_observable_mapping='A_i=P_i^2 for the three global Cartesian polarization components; arbitrary axis-dependent model attempt, not physical species or inferred saddles',GLV_coefficients=coeff,GLV_forecasts=glv,
        unrestricted_GLV_has_negative_competition_coefficients=bool((coeff[:,1:]<0).any()),constrained_activity_coefficients=constrained,constrained_activity_forecasts=constrained_forecast,
        model_comparison_note='Trajectories 1 and 2 are excluded from coefficient fitting. This is exploratory model comparison, not an untouched final validation set: the nonnegative-competition/source variant was added after inspecting the unrestricted fit.',
        limitation='Clustering partitions continuous data by construction. Recurrence, persistence or a useful forecast does not identify invariant saddle states or heteroclinic connecting orbits.')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--partial',action='store_true');args=parser.parse_args()
    protocol=json.loads((D/'protocol.json').read_text());planned=cases();loaded={};rows=[];missing=[]
    for kind in ['water','rotors']:
        for rep in range(3):
            for case in planned:
                path=D/f'{kind}_r{rep}_{case["name"]}.npz'
                if not path.exists():missing.append(path.name);continue
                with np.load(path) as archive:data={k:archive[k] for k in archive.files}
                conf=json.loads(str(data['configuration_json']));assert conf['case']==case
                for k,a in data.items():
                    if a.dtype.kind in 'fci':assert np.isfinite(a).all(),(path,k)
                loaded[kind,rep,case['name']]=data
                t=data['time_ps'];late=t>=t[-1]*2/3;p=data['polarization'];e=data['field_V_nm'];norm=np.linalg.norm(e,axis=1)
                unit=e/np.maximum(norm[:,None],1e-15);projection=np.sum(p*unit,axis=1)
                row=dict(kind=kind,replica=rep,case=case['name'],duration_ps=float(t[-1]),late_Px=float(p[late,0].mean()),late_norm=float(np.linalg.norm(p[late],axis=1).mean()),late_field_projection=float(projection[late].mean()),
                    late_temperature_K=float(data['temperature_K'][late].mean()),late_hbond_neighbors=float(data['hbond_neighbors'][late].mean()),late_connected_orientation=float(data['connected_neighbor_orientation'][late].mean()),
                    late_nematic_order=float(data['nematic_order'][late].mean()),maximum_temperature_K=float(data['temperature_K'].max()))
                row['late_spatial_excess']=float(np.mean(data['connected_neighbor_orientation'][late]+(1-np.sum(p[late]*p[late],axis=1))/215))
                if case['mode']=='rotate':
                    phase_window=t>=max(t[-1]/4,1/case['frequency'])
                    row['phase_window_ps']=[float(t[phase_window][0]),float(t[-1])]
                    row['phase']=phase_metrics(t[phase_window],p[phase_window],case['frequency'])
                    row['phase_threshold_sensitivity']=[phase_metrics(t[phase_window],p[phase_window],case['frequency'],threshold) for threshold in [.1,.2]]
                if case['mode']=='sequence':
                    off=(t>=case['training']+6)&(t<=case['duration']);row['teacher_off']=sequence_metrics(t[off]-case['training'],p[off])
                    on=(t>=6)&(t<case['training']);row['teacher_on']=sequence_metrics(t[on],p[on])
                if case['mode']=='ramp':row['absolute_loop_area']=float(abs(np.trapezoid(p[:,0],e[:,0])))
                if case['mode']=='pulse':
                    row['mean_end_of_pulse_alignment']=float(p[np.isclose(t%8,1),0].mean())
                    row['mean_end_of_recovery_alignment']=float(p[np.isclose(t%8,7.9),0].mean())
                if case['mode']=='feedback':
                    q=p[late,0]+1j*p[late,1];frequency=np.fft.fftfreq(len(q),.1);power=abs(np.fft.fft(q-q.mean()))**2
                    j=int(power.argmax());row['feedback_peak_frequency_ps_inverse']=float(frequency[j]);row['feedback_peak_power_fraction']=float(power[j]/max(power.sum(),1e-30))
                    # Each recorded E[k+1] is the new controller value applied at t[k].
                    x=data['snapshot_positions_nm'][0];oh=x[:,1:]-x[:,0:1];side=protocol['box_side_nm'];oh-=side*np.rint(oh/side)
                    mu=float(np.linalg.norm(.417*oh.sum(axis=1),axis=1).mean())
                    row['controller_work_kJ_mol_per_molecule']=float(-96.48533212331002*mu*np.sum(p[:-1]*np.diff(e,axis=0)))
                rows.append(row)
    if missing and not args.partial:raise RuntimeError(f'{len(missing)} required runs missing; use --partial only for progress inspection')
    aggregates=[]
    keys=['late_Px','late_norm','late_field_projection','late_temperature_K','late_hbond_neighbors','late_connected_orientation','late_nematic_order','late_spatial_excess']
    for kind in ['water','rotors']:
        for case in planned:
            subset=[r for r in rows if r['kind']==kind and r['case']==case['name']]
            if subset:aggregates.append(dict(kind=kind,case=case['name'],**{key:summary([r[key] for r in subset]) for key in keys},tracking_count=sum(r.get('phase',{}).get('finite_window_tracking',False) for r in subset)))
    lookup={(r['kind'],r['replica'],r['case']):r for r in rows};contrasts={}
    for kind in ['water','rotors']:
        for name,a,b,metric in [('history_positive_minus_negative','history_positive','history_negative','late_Px'),('slow_minus_fast_loop','ramp_slow','ramp_fast','absolute_loop_area'),('feedback_align_minus_baseline','feedback_align','baseline','late_norm')]:
            vals=[]
            for rep in range(3):
                aa=lookup.get((kind,rep,a));bb=lookup.get((kind,rep,b))
                if aa and bb:vals.append(aa[metric]-bb[metric])
            if vals:contrasts[kind+'_'+name]=summary(vals)
    for case in planned:
        vals=[]
        for rep in range(3):
            a=lookup.get(('water',rep,case['name']));b=lookup.get(('rotors',rep,case['name']))
            if a and b:vals.append(a['late_field_projection']-b['late_field_projection'])
        if vals:contrasts['interaction_projection_'+case['name']]=summary(vals)
    sequence=[]
    for kind in ['water','rotors']:
        for rep in range(3):
            for name in ['sequence_xyz','sequence_scrambled','baseline']:
                d=loaded.get((kind,rep,name))
                if d is None:continue
                if name=='baseline':mask=d['time_ps']>=d['time_ps'][-1]-24;metric=sequence_metrics(d['time_ps'][mask],d['polarization'][mask])
                else:metric=lookup[kind,rep,name]['teacher_off']
                sequence.append(dict(kind=kind,replica=rep,case=name,**metric))
    refinements=[]
    for kind in ['water','rotors']:
        for rep in range(3):
            base=lookup.get((kind,rep,'rotate_a0.5_f1'))
            if base is None:continue
            for suffix in ['half_step','half_update']:
                path=D/f'{kind}_r{rep}_rotate_a0.5_f1_{suffix}.npz'
                if path.exists():
                    with np.load(path) as d:
                        t=d['time_ps'];mask=t>=max(t[-1]/4,1.);pm=phase_metrics(t[mask],d['polarization'][mask],1.)
                        refinements.append(dict(kind=kind,replica=rep,change=suffix,phase=pm,base_tracking=base['phase']['finite_window_tracking'],response_magnitude_difference=pm['response_magnitude']-base['phase']['response_magnitude'],temperature_K=float(d['temperature_K'][t>=t[-1]*2/3].mean())))
    analysis_checks=dict(phase_locked_fixture=phase_metrics(np.arange(0,10,.1),np.column_stack([.4*np.cos(2*np.pi*np.arange(0,10,.1)-.3),.4*np.sin(2*np.pi*np.arange(0,10,.1)-.3),np.zeros(100)]),1.)['finite_window_tracking'],
        low_amplitude_fixture_rejected=not phase_metrics(np.arange(0,10,.1),np.column_stack([.01*np.cos(2*np.pi*np.arange(0,10,.1)),.01*np.sin(2*np.pi*np.arange(0,10,.1)),np.zeros(100)]),1.)['finite_window_tracking'])
    assert all(analysis_checks.values())
    result=dict(expected_main_runs=156,completed_main_runs=len(rows),complete=not missing,missing=missing,rows=rows,aggregates=aggregates,contrasts=contrasts,
        sequence_recall=sequence,recurrence=recurrence(loaded),refinements=refinements,analysis_checks=analysis_checks,
        decision_scope='Exploratory finite-time screening. Nominal t intervals are not adjusted for the many comparisons. No water bifurcation, autonomous oscillation, heteroclinic network or learning claim follows from one positive metric.')
    save(D/('partial_results.json' if args.partial else 'results.json'),result)
    if not args.partial:
        with (D/'run_summary.csv').open('w',newline='',encoding='utf-8') as stream:
            columns=['kind','replica','case','duration_ps']+keys;writer=csv.DictWriter(stream,columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(completed=len(rows),expected=156,refinements=len(refinements),analysis_checks=analysis_checks,contrasts=contrasts)),flush=True)


if __name__=='__main__':main()
