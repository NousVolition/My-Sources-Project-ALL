"""Independent stationary baselines and paired interventions; no labels used in setup."""
import argparse,json,time
from pathlib import Path
import numpy as np
from solver import Flow,initial,impulse,resize,observe
from common import ROOT,digest,save,read,protocol,sources,now

def basic(f,h):
 density=np.sum(abs(h)**2,axis=-4)*f.weights/f.n**6/2
 E=density.sum(axis=(-3,-2,-1));eps=2*f.nu*np.sum(density*f.k2,axis=(-3,-2,-1))
 return E,eps,f.cut*(f.nu**3/eps)**.25,np.sum(density*f.high,axis=(-3,-2,-1))/E

def response(f,h,config):
 strength=np.array(config['strengths']);speed=np.sqrt(f.inner(h,h));hs=np.stack([h]+[h+impulse(f,a*speed,config['impulse_width']) for a in strength]);initial_h=hs.copy()
 E0=.5*f.inner(hs,hs);delta=hs[1:]-hs[0];de0=.5*f.inner(delta,delta)
 _,_,rates,_=f.parts(hs)
 budget=np.zeros((len(hs),2));pdq=np.zeros((len(hs)-1,3));curves=[];energies=[];budgets=[];diffbudgets=[];times=[];diagnostics=[];cfl=0
 dt=config['response_dt'];steps=round(config['horizon']/dt);stride=round(config['save_dt']/dt)
 for step in range(steps+1):
  if step%stride==0:
   E,eps,keta,tail=basic(f,hs);d=hs[1:]-hs[0];de=.5*f.inner(d,d)
   curves.append(de/de0);energies.append(E);budgets.append(budget.copy());diffbudgets.append(pdq.copy()/de0[:,None]);times.append(step*dt)
   diagnostics.append([float(keta.min()),float(tail.max()),float(np.max(abs(E-E0-budget[:,0]+budget[:,1])/E0)),float(np.max(abs(de-de0-pdq[:,0]+pdq[:,1]-pdq[:,2])/de0))])
  if step==steps:break
  hs,b,d,c=f.step(hs,dt);budget+=b;pdq+=d;cfl=max(cfl,c)
  if not np.isfinite(hs).all() or c>1.2:raise RuntimeError('Response failed stability guard')
 return dict(response_time=np.array(times),energy_ratio=np.array(curves),total_energy=np.array(energies),total_budget=np.array(budgets),difference_budget=np.array(diffbudgets),response_checks=np.array(diagnostics),initial_rates=rates/de0[:,None],initial_difference_energy=de0,response_cfl=np.array(cfl),initial_state_sha256=np.array(__import__('hashlib').sha256(initial_h.tobytes()).hexdigest()),final_state_sha256=np.array(__import__('hashlib').sha256(hs.tobytes()).hexdigest()))

def run(seed,label='base',n=None,dt=None,screen=False):
 p=protocol();c={k:p[k] for k in ['n','nu','power','spin_dt','response_dt','events','horizon','save_dt','baseline_save_dt','strengths','impulse_width','initial_speed','history_lags']}
 c['n']=n or c['n'];c['spin_dt']=dt or c['spin_dt'];c['response_dt']=min(c['response_dt'],c['spin_dt'])
 if label=='half':c['response_dt']=p['response_dt']/2
 folder=ROOT/'data'/('pilot' if screen else 'runs');folder.mkdir(parents=True,exist_ok=True);stem=f'{label}_{seed}';npz=folder/(stem+'.npz');meta=folder/(stem+'.json');src=sources()
 if npz.exists():
  old=json.loads(meta.read_text(encoding='utf-8'));assert old['config']==c and old['source_sha256']==src and old['sha256']==digest(npz);print('Cached '+stem,flush=True);return
 if not screen and seed in p['splits']['test']:
  from forecast import verify_freeze
  verify_freeze()
 started=now();start=time.time();f=Flow(c['n'],c['nu'],c['power']);h=initial(f,seed,c['initial_speed'])[None];e0=float(f.inner(h,h)[0]/2);budget=np.zeros((1,2));rows=[];feature_rows=[];spectra=[];times=[];states=[];state_times=[];event_results=[];history=[];cfl=0
 stepdt=c['spin_dt'];end=max(c['events']);steps=round(end/stepdt);stride=round(c['baseline_save_dt']/stepdt)
 for step in range(steps+1):
  t=step*stepdt
  if step%stride==0:
   obs,stats,sp=observe(f,h[0]);stats.update(t=t,cfl=cfl,balance=float((stats['energy']-e0-budget[0,0]+budget[0,1])/e0));times.append(t);feature_rows.append(obs);spectra.append(sp);rows.append(stats)
   if seed in p['refinement_seeds'] and any(abs(t-(event-max(c['history_lags'])))<1e-8 for event in p['refinement_events']):states.append(h[0].copy());state_times.append(t)
   if any(abs(t-event)<1e-8 for event in c['events']):
    index=[int(np.argmin(abs(np.array(times)-(t-lag)))) for lag in [0]+c['history_lags']]
    history.append(np.array(feature_rows)[index]);event_results.append(response(f,h[0],c))
   if abs(t/4-round(t/4))<1e-9:
    save(str(folder.relative_to(ROOT)/(stem+'_progress.json')),{'seed':seed,'t':t,'elapsed_seconds':time.time()-start,'latest':stats})
    print(json.dumps({'seed':seed,'label':label,'time':t,'elapsed':round(time.time()-start),'Re_lambda':round(stats['re_lambda'],3)}),flush=True)
  if step==steps:break
  h,b,_,cf=f.step(h,stepdt);budget+=b;cfl=max(cfl,cf)
  if not np.isfinite(h).all() or cf>1.2:raise RuntimeError('Baseline failed stability guard')
 arrays={k:np.array([v[k] for v in event_results]) for k in event_results[0]};arrays.update(baseline_time=np.array(times),baseline_features=np.array(feature_rows),spectra=np.array(spectra),history=np.array(history),strengths=np.array(c['strengths']),events=np.array(c['events']))
 # Keep full restart fields only on predeclared refinement seeds; no test labels used.
 if seed in p['refinement_seeds'] or screen:
  checkpoint=folder/(stem+'_checkpoints.npz');np.savez_compressed(checkpoint,states=np.array(states),times=np.array(state_times));checkpoint_info={'path':checkpoint.name,'sha256':digest(checkpoint)}
 else:checkpoint_info=None
 np.savez_compressed(npz,**arrays)
 save(str(meta.relative_to(ROOT)),{'seed':seed,'label':label,'config':c,'source_sha256':src,'started_utc':started,'completed_utc':now(),'elapsed_seconds':time.time()-start,'sha256':digest(npz),'baseline':rows,'checkpoints':checkpoint_info})
 print('Completed '+stem,flush=True)

def refine(seed,label):
 p=protocol();base=ROOT/'data/runs'/f'base_{seed}';original=read(str(base.relative_to(ROOT))+'.json');assert original['sha256']==digest(str(base)+'.npz')
 c=original['config'].copy();c['events']=p['refinement_events'];c['n']=p['fine_n'] if label=='fine' else p['n'];c['spin_dt']=p['spin_dt']/2;c['response_dt']=p['response_dt']/2
 f=Flow(c['n'],c['nu'],c['power']);src=sources();out=ROOT/'data/runs'/f'{label}_{seed}.npz';checkpoint=ROOT/'data/runs'/original['checkpoints']['path'];assert digest(checkpoint)==original['checkpoints']['sha256']
 if out.exists():
  q=read(str(out.with_suffix('.json').relative_to(ROOT)));assert q['source_sha256']==src and q['sha256']==digest(out);print('Cached '+out.name,flush=True);return
 start=time.time();results=[];histories=[];prechecks=[]
 with np.load(checkpoint) as z:states=z['states'];state_times=z['times']
 for j,(state,t0) in enumerate(zip(states,state_times)):
  h=f.project(resize(state,c['n']))[None];features=[];tt=[];stats=[];cfl=0;E0=float(f.inner(h,h)[0]/2);budget=np.zeros((1,2))
  for step in range(round((c['events'][j]-t0)/c['spin_dt'])+1):
   t=t0+step*c['spin_dt']
   if step%round(c['baseline_save_dt']/c['spin_dt'])==0:
    obs,s,_=observe(f,h[0]);s.update(t=t,cfl=cfl,balance=float((s['energy']-E0-budget[0,0]+budget[0,1])/E0));features.append(obs);tt.append(t);stats.append(s)
   if abs(t-c['events'][j])<1e-8:break
   h,b,_,cf=f.step(h,c['spin_dt']);budget+=b;cfl=max(cfl,cf)
  idx=[int(np.argmin(abs(np.array(tt)-(c['events'][j]-lag)))) for lag in [0]+c['history_lags']]
  histories.append(np.array(features)[idx]);results.append(response(f,h[0],c));prechecks.extend(stats)
  print(f'{label} {seed} event {j+1} complete',flush=True)
 arrays={k:np.array([r[k] for r in results]) for k in results[0]};arrays.update(history=np.array(histories),strengths=np.array(c['strengths']),events=np.array(c['events']))
 np.savez_compressed(out,**arrays);save(str(out.with_suffix('.json').relative_to(ROOT)),{'seed':seed,'label':label,'config':c,'source_sha256':src,'sha256':digest(out),'baseline':prechecks,'elapsed_seconds':time.time()-start,'checkpoint_sha256':digest(checkpoint),'comparison':'Matched restart from base grid at each event minus longest history lag; not an independent realization or full-spin convergence test'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('seed',type=int);ap.add_argument('--label',default='base');ap.add_argument('--screen',action='store_true');a=ap.parse_args()
 if a.label=='base':run(a.seed,screen=a.screen)
 else:refine(a.seed,a.label)
