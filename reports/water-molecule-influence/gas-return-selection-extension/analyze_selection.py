"""Past-only selection and analogue forecasting, with start-separated tuning."""
from pathlib import Path
import json, hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
SRC=HERE.parent/'gas-damping-extension/data'
PROTOCOL=json.loads((DATA/'protocol.json').read_text())
LEVELS=[.2,1,10]
H=np.array(PROTOCOL['horizons_ps'])
GRID=PROTOCOL['velocity_time_grid_ps']
METHODS=['position','direction','velocity','hold','linear']
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False,
        default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())+'\n')


def interp(p,t,dt):
    t=np.asarray(t);q=t/dt
    if np.any(q < -1e-8) or np.any(q > len(p)-1+1e-8):raise ValueError('Extrapolation')
    q=np.clip(q,0,len(p)-1);j=np.minimum(np.floor(q).astype(int),len(p)-2);f=q-j
    return p[j]*(1-f[...,None])+p[j+1]*f[...,None]


def vel(p,t,dt):return (interp(p,t,dt)-interp(p,np.asarray(t)-dt,dt))/dt


def candidate_batch(p,ids,dt):
    a=ids*dt;offset=np.arange(round(.1/dt)+1,round(1/dt)+1)
    j=ids[:,None]-offset
    base=p[j];segment=p[j+1]-base
    l2=np.sum(segment**2,axis=-1)
    frac=np.clip(np.sum((p[ids,None]-base)*segment,axis=-1)/np.maximum(l2,1e-30),0,1)
    old=base+frac[...,None]*segment;b=(j+frac)*dt
    ov=vel(p,b,dt);cv=vel(p,a,dt)
    distance2=np.sum((p[ids,None]-old)**2,axis=-1)
    speed_old=np.linalg.norm(ov,axis=-1);speed_current=np.linalg.norm(cv,axis=-1)
    valid=(speed_old>1e-12)&(speed_current[:,None]>1e-12)
    cos=np.divide(np.sum(cv[:,None]*ov,axis=-1),speed_current[:,None]*speed_old,
        out=np.full_like(distance2,-1),where=valid)
    return dict(a=a,b=b,old=old,ov=ov,cv=cv,d2=distance2,
        direction2=2-2*np.clip(cos,-1,1),velocity2=np.sum((cv[:,None]-ov)**2,axis=-1),
        valid_speed=valid,eligible=distance2<=.02**2)


def evaluate(p,dt,settings,ids=None):
    if ids is None:ids=np.rint(np.arange(1000,1990)*.01/dt).astype(int)
    out={key:dict(time=[],old_time=[],initial_distance=[],initial_angle=[],initial_speed_ratio=[],initial_velocity_defined=[],
        forecast=[],error2=[],paired_distance=[],paired_angle=[]) for key in settings}
    for batch in np.array_split(ids,max(1,int(np.ceil(len(ids)/100)))):
        c=candidate_batch(p,batch,dt);keep=c['eligible'].any(axis=1)
        if not keep.any():continue
        a=c['a'][keep];current=interp(p,a,dt);cv=c['cv'][keep]
        truth=interp(p,a[:,None]+H,dt);future_v=vel(p,a[:,None]+H,dt)
        for key,(method,weight) in settings.items():
            penalty=0 if method=='position' else weight**2*c[method+'2']
            score=np.where(c['eligible'],c['d2']+penalty,np.inf)
            picks=np.argmin(score,axis=1)[keep];rr=np.flatnonzero(keep)
            b=c['b'][rr,picks];old=c['old'][rr,picks];ov=c['ov'][rr,picks]
            assert np.all(b+.1<=a+1e-10)
            past_continuation=interp(p,b[:,None]+H,dt)
            prediction=current[:,None]+past_continuation-old[:,None]
            old_future_v=vel(p,b[:,None]+H,dt)
            denom=np.linalg.norm(future_v,axis=-1)*np.linalg.norm(old_future_v,axis=-1)
            cosine=np.divide(np.sum(future_v*old_future_v,axis=-1),denom,
                out=np.full_like(denom,-1),where=denom>1e-20)
            o=out[key]
            values=dict(time=a,old_time=b,initial_distance=np.sqrt(c['d2'][rr,picks]),
                initial_angle=np.rad2deg(np.arccos(np.clip(1-c['direction2'][rr,picks]/2,-1,1))),
                initial_speed_ratio=np.linalg.norm(cv,axis=-1)/np.maximum(np.linalg.norm(ov,axis=-1),1e-12),
                initial_velocity_defined=c['valid_speed'][rr,picks],
                forecast=prediction,error2=np.sum((prediction-truth)**2,axis=-1),
                paired_distance=np.linalg.norm(past_continuation-truth,axis=-1),
                paired_angle=np.rad2deg(np.arccos(np.clip(cosine,-1,1))))
            for name,value in values.items():o[name].append(value)
    for key in out:
        for name in out[key]:out[key][name]=np.concatenate(out[key][name])
    a=out[next(iter(out))]['time'];current=interp(p,a,dt);cv=vel(p,a,dt);truth=interp(p,a[:,None]+H,dt)
    for name,prediction in [('hold',np.repeat(current[:,None,:],len(H),axis=1)),('linear',current[:,None]+H[None,:,None]*cv[:,None])]:
        out[name]=dict(time=a,forecast=prediction,error2=np.sum((prediction-truth)**2,axis=-1))
    return out


def summarize(out,mask=None):
    summaries={}
    for key,values in out.items():
        keep=np.ones(len(values['time']),dtype=bool) if mask is None else mask
        summary=dict(n=int(keep.sum()),mse=values['error2'][keep].mean(axis=0))
        if 'old_time' in values:
            same=(values['paired_distance']<=.02)&(values['paired_angle']<=45)
            summary.update(near_same_fraction=same[keep].mean(axis=0),
                undefined_initial_velocities=int((~values['initial_velocity_defined'][keep]).sum()),
                median_initial_distance=float(np.median(values['initial_distance'][keep])),
                median_initial_angle=float(np.median(values['initial_angle'][keep])),
                median_initial_speed_ratio=float(np.median(values['initial_speed_ratio'][keep])))
        summaries[key]=summary
    return summaries


def checks():
    dt=.01;t=np.arange(0,5+dt/2,dt);p=np.column_stack([.002*t,.003*t])
    ids=np.array([200,201,202]);settings={'position':('position',0),'velocity':('velocity',.01),'direction':('direction',.01)}
    result=evaluate(p,dt,settings,ids)
    assert np.max(result['position']['error2'])<1e-25
    assert np.max(result['linear']['error2'])<1e-25
    assert np.min(result['hold']['error2'])>0
    assert np.allclose(vel(p,2.315,dt),[.002,.003])
    # Changing any data after the current anchor must leave the chosen old
    # passage and every forecast unchanged; only evaluated error may change.
    rng=np.random.default_rng(73);p=np.cumsum(rng.normal(0,.001,size=(501,2)),axis=0)
    old=evaluate(p,dt,settings,np.array([200]));q=p.copy();q[201:]+=7
    changed=evaluate(q,dt,settings,np.array([200]))
    for key in old:
        assert np.array_equal(old[key]['forecast'],changed[key]['forecast'])
        if 'old_time' in old[key]:assert np.array_equal(old[key]['old_time'],changed[key]['old_time'])
        assert not np.allclose(old[key]['error2'],changed[key]['error2'])
    return dict(passed=True,constant_velocity_analogue=True,linear_baseline=True,
        backward_fractional_derivative=True,future_mutation_leaves_selection_and_forecast_unchanged=True)


def main():
    inputs={};paths={}
    for rep in range(3):
        for g in LEVELS:
            path=SRC/f'vapor_r{rep}_gamma{g:g}.npz';inputs[path.name]=sha(path)
            with np.load(path) as z:paths[rep,g]=z['state'][:,:2]
    candidates={'position':('position',0)}
    for method in ['direction','velocity']:
        for w in GRID:candidates[f'{method}_{w:g}']=(method,w)
    train={g:evaluate(paths[0,g],.01,candidates) for g in LEVELS}
    grid=[];chosen={}
    for method in ['direction','velocity']:
        candidates_scores=[]
        for w in GRID:
            errors=[float(train[g][f'{method}_{w:g}']['error2'][:,0].mean()) for g in LEVELS]
            row=dict(method=method,weight=w,per_damping_mse=errors,mean_mse=float(np.mean(errors)))
            candidates_scores.append(row);grid.append(row)
        chosen[method]=min(candidates_scores,key=lambda s:(s['mean_mse'],s['weight']))['weight']
    lock=dict(training_start=1,training_horizon_ps=.02,weights=chosen,grid=grid,
        script_sha256=sha(__file__),protocol_sha256=sha(DATA/'protocol.json'))
    save(DATA/'locked_training.json',lock)
    settings={'position':('position',0),**{m:(m,w) for m,w in chosen.items()}}
    rows=[];coarse={}
    for rep in range(3):
        for g in LEVELS:
            result=evaluate(paths[rep,g],.01,settings);coarse[rep,g]=result
            row=dict(start=rep+1,friction_per_ps=g,role='training' if rep==0 else 'test',
                eligible=len(result['position']['time']),anchors=990,methods=summarize(result))
            rows.append(row)
            arrays={f'{method}__{key}':value for method,values in result.items() for key,value in values.items()}
            np.savez_compressed(DATA/f'predictions_r{rep}_gamma{g:g}.npz',**arrays)
    fine=[]
    for g in LEVELS:
        path=SRC/f'vapor_r1_gamma{g:g}_fine.npz';inputs[path.name]=sha(path)
        with np.load(path) as z:p=z['state'][:,:2]
        assert np.array_equal(p[::10],paths[1,g])
        result=evaluate(p,.001,settings)
        ct=np.rint(coarse[1,g]['position']['time']*100).astype(int)
        ft=np.rint(result['position']['time']*100).astype(int)
        common=np.intersect1d(ct,ft)
        fine.append(dict(start=2,friction_per_ps=g,native=summarize(result),
            common_anchors=len(common),fine_common=summarize(result,np.isin(ft,common)),
            coarse_common=summarize(coarse[1,g],np.isin(ct,common))))
        np.savez_compressed(DATA/f'predictions_r1_gamma{g:g}_fine.npz',
            **{f'{method}__{key}':value for method,values in result.items() for key,value in values.items()})
    results=dict(rows=rows,fine_recording_check=fine,input_sha256=inputs,
        training=lock,validation=checks())
    save(DATA/'results.json',results)
    print(json.dumps(dict(weights=chosen,test=[dict(start=r['start'],g=r['friction_per_ps'],n=r['eligible'],
        mse={k:v['mse'].tolist() for k,v in r['methods'].items()}) for r in rows if r['role']=='test']),indent=2))


if __name__=='__main__':main()
