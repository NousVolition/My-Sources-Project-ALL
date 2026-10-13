"""Prospective viscosity/grid pilots, with no training or test observations."""
import argparse,time
import numpy as np
from solver import Flow,initial,observe,resize
from simulate import response
from common import ROOT,protocol,save,digest,sources,now

def baseline(case):
    p=protocol();c=p['pilot_cases'][case];f=Flow(c['n'],p['nu'],p['power'])
    h=initial(f,c['seed'],p['initial_speed'])[None];e0=float(f.inner(h,h)[0]/2)
    budget=np.zeros((1,2));cfl=0.;rows=[];spectra=[];times=[];start=time.time()
    stem=ROOT/'data/pilot'/f"pilot_{case}"
    stem.parent.mkdir(parents=True,exist_ok=True)
    dt=c['dt'];steps=round(p['pilot_end']/dt);stride=round(p['baseline_save_dt']/dt)
    src={**sources(),'pilot.py':digest(ROOT/'pilot.py')};started=now()
    for step in range(steps+1):
        t=step*dt
        if step%stride==0:
            _,s,sp=observe(f,h[0]);s.update(t=t,cfl=cfl,balance=float((s['energy']-e0-budget[0,0]+budget[0,1])/e0))
            rows.append(s);spectra.append(sp);times.append(t)
            if step%round(4/dt)==0:
                save(f'data/pilot/pilot_{case}_progress.json',dict(case=case,config=c,rows=rows,source_sha256=src,elapsed_seconds=time.time()-start))
                print(dict(case=case,t=t,re_lambda=round(s['re_lambda'],3),keta=round(s['keta'],4),cfl=round(cfl,4),elapsed=round(time.time()-start)),flush=True)
        if step==steps:break
        h,b,_,cf=f.step(h,dt);budget+=b;cfl=max(cfl,cf)
        if not np.isfinite(h).all() or cf>p['gates']['cfl_max']:raise RuntimeError('Pilot CFL or finite-value rejection')
    np.savez_compressed(str(stem)+'.npz',final=h[0],spectra=np.asarray(spectra),time=np.asarray(times))
    cut=f.cut;g=p['gates']
    checks=dict(cfl=cfl<=g['cfl_max'],diffusion=3*p['nu']*cut**2*dt<=g['diffusion_max'],resolution=min(r['keta'] for r in rows)>=g['keta_min'],tail=max(r['tail_energy'] for r in rows)<=g['tail_energy_max'],divergence=max(r['divergence'] for r in rows)<=g['divergence_max'],total_budget=max(abs(r['balance']) for r in rows)<=g['budget_relative_max'])
    save(f'data/pilot/pilot_{case}.json',dict(case=case,config=c,baseline=rows,checks=checks,started_utc=started,completed_utc=now(),source_sha256=src,sha256=digest(str(stem)+'.npz'),elapsed_seconds=time.time()-start))
    print(dict(case=case,checks=checks),flush=True)
    if c['n']==p['n'] and c['seed']==p['response_pilot']['seed']:
        for label in ['base','half','fine']:pilot_response(label)

def pilot_response(label):
    p=protocol();case=p['response_pilot']['case'];source=ROOT/'data/pilot'/f'pilot_{case}.npz'
    c=p.copy();c['n']=p['fine_n'] if label=='fine' else p['n'];c['response_dt']=p['response_dt']/(1 if label=='base' else 2)
    with np.load(source) as z:h=resize(z['final'],c['n'])
    f=Flow(c['n'],c['nu'],c['power']);start=time.time();started=now();z=response(f,f.project(h),c)
    out=ROOT/'data/pilot'/f'response_{label}.npz';np.savez_compressed(out,**z)
    q=z['response_checks'];g=p['gates']
    checks=dict(cfl=float(z['response_cfl'])<=g['cfl_max'],resolution=float(q[:,0].min())>=g['keta_min'],tail=float(q[:,1].max())<=g['tail_energy_max'],total_budget=float(q[:,2].max())<=g['budget_relative_max'],difference_budget=float(q[:,3].max())<=g['budget_relative_max'])
    save(f'data/pilot/response_{label}.json',dict(config={k:c[k] for k in ['n','nu','power','response_dt','horizon','strengths','impulse_width']},input_path=source.relative_to(ROOT).as_posix(),input_sha256=digest(source),source_sha256={**sources(),'pilot.py':digest(ROOT/'pilot.py')},started_utc=started,completed_utc=now(),elapsed_seconds=time.time()-start,sha256=digest(out),checks=checks,minimum_keta=float(q[:,0].min()),maximum_cfl=float(z['response_cfl'])))
    print(label,checks,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('case',type=int);baseline(a.parse_args().case)
