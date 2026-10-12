import json
import numpy as np
from common import ROOT,read,save,protocol,digest

def numerical(path,meta,p):
 with np.load(path) as z:
  rc=z['response_checks'];cfl=float(z['response_cfl'].max());finite=all(np.isfinite(z[k]).all() for k in z.files if z[k].dtype.kind in 'fci')
 rows=meta['baseline'];g=p['gates'];c=meta['config'];cut=int(np.ceil(c['n']/3)-1)*2*np.pi/6
 return dict(finite=bool(finite),cfl=max(cfl,max(r['cfl'] for r in rows))<=g['cfl_max'],diffusion=3*c['nu']*cut**2*max(c['spin_dt'],c['response_dt'])<=g['diffusion_max'],
             total_budget=max(float(rc[...,2].max()),max(abs(r['balance']) for r in rows))<=g['budget_relative_max'],difference_budget=float(rc[...,3].max())<=g['budget_relative_max'],
             resolution=min(float(rc[...,0].min()),min(r['keta'] for r in rows))>=g['keta_min'],tail=max(float(rc[...,1].max()),max(r['tail_energy'] for r in rows))<=g['tail_energy_max'],
             divergence=max(r['divergence'] for r in rows)<=g['divergence_max'])

def screening():
 p=protocol();records={};hashes={};arrays={}
 for label in ['base','half','fine']:
  path=ROOT/'data/screen'/f'response_{label}_300000.npz';q=read(str(path.with_suffix('.json').relative_to(ROOT)));assert digest(path)==q['sha256']
  for key in ['nu','power','horizon','strengths','impulse_width']:assert q['config'][key]==p[key]
  assert q['config']['n']==(p['fine_n'] if label=='fine' else p['n'])
  assert q['config']['response_dt']==p['response_dt']/(1 if label=='base' else 2)
  assert q['input_path']=='data/screen/n64_nu0.025_dt0.025_seed300000.npz'
  assert digest(ROOT/q['input_path'])==q['input_sha256']
  for source,sha in q['source_sha256'].items():assert digest(ROOT/source)==sha,source
  hashes[q['input_path']]=q['input_sha256']
  records[label]=q['checks'];hashes[path.relative_to(ROOT).as_posix()]=digest(path);hashes[path.with_suffix('.json').relative_to(ROOT).as_posix()]=digest(path.with_suffix('.json'))
  with np.load(path) as z:arrays[label]=z['energy_ratio']
 step=float(np.max(abs(arrays['half']-arrays['base']))/np.max(abs(arrays['base'])));grid=float(np.max(abs(arrays['fine']-arrays['half']))/np.max(abs(arrays['half'])))
 return dict(checks=records,step_curve_change=step,grid_curve_change=grid,sha256=hashes,passed=bool(all(all(q.values()) for q in records.values()) and step<=p['gates']['step_curve_change'] and grid<=p['gates']['grid_curve_change']))

def verify_all(write=True):
 p=protocol();checks={};refinements=[]
 for split in ['train','validation','test']:
  for seed in p['splits'][split]:
   file=ROOT/'data/runs'/f'base_{seed}.npz';m=read(f'data/runs/base_{seed}.json');assert digest(file)==m['sha256'];checks[file.name]=numerical(file,m,p)
 for seed in p['refinement_seeds']:
  zs={}
  for label in ['base','half','fine']:
   file=ROOT/'data/runs'/f'{label}_{seed}.npz';m=read(f'data/runs/{label}_{seed}.json');assert digest(file)==m['sha256'];checks[file.name]=numerical(file,m,p)
   with np.load(file) as z:zs[label]={k:z[k] for k in z.files}
  # Half and fine have identical smaller timesteps: the mesh contrast is isolated.
  for label,reference in [('half','base'),('fine','half')]:
   a=zs[label]['energy_ratio'];b=zs[reference]['energy_ratio'];relative=float(np.max(abs(a-b))/max(float(np.max(abs(b))),1e-12))
   from features import target,inputs
   ya=target(zs[label]);yb=target(zs[reference]);inputs_a=inputs(zs[label]);inputs_b=inputs(zs[reference])
   refinements.append(dict(seed=seed,label=label,reference=reference,max_curve_relative_change=relative,max_target_relative_change=float(np.max(abs(ya-yb)/abs(yb))),passed=relative<=p['gates']['step_curve_change' if label=='half' else 'grid_curve_change']))
 screen=screening();result=dict(numerical=checks,refinements=refinements,screening=screen,all_passed=bool(screen['passed'] and all(all(q.values()) for q in checks.values()) and all(q['passed'] for q in refinements)))
 if write:save('results/numerical_checks.json',result)
 print(json.dumps({'all_passed':result['all_passed'],'runs_checked':len(checks),'refinement_pairs':len(refinements)}));return result
if __name__=='__main__':verify_all()
