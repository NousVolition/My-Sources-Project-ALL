"""Summarize retained screening cases separately from the independent prediction ensemble."""
import csv,json
import numpy as np
from common import ROOT,read,save,digest

def main():
 files=['solver.py','screen_original_solver.py','screen_optimized_solver.py','numpy_solver_reference.py','numexpr_solver_reference.py']
 known={digest(ROOT/f):f for f in files};rows=[];provenance=[]
 failed={'n80_nu0.02_dt0.025_seed300003':'Stopped at the initial-step CFL guard (>1); not a completed integration'}
 for path in sorted((ROOT/'data/screen').glob('n*.json')):
  q=json.loads(path.read_text(encoding='utf-8'));c=q['config'];tag=f"n{c['n']}_nu{c['nu']}_dt{c['dt']}_seed{c['seed']}";recording=ROOT/'data/screen'/(tag+'.npz');complete=recording.exists() and q['rows'][-1]['t']>=c['end']-1e-9
  status='completed' if complete else failed.get(tag,'unfinished; inspect retained log before publication')
  window=[r for r in q['rows'] if 20<=r['t']<45];source=q.get('source_sha256');module=known.get(source)
  note='Source hash captured at launch.'
  if path.name=='n48_nu0.025_dt0.json':
   module='screen_original_solver.py';note='First driver hashed the on-disk solver at completion, after optimization, rather than the module loaded at launch. The executed pre-optimization source was preserved before editing; its actual hash is supplied separately here. Its JSON filename also reflects a with_suffix naming bug. Numerical recordings are retained unchanged.'
  provenance.append(dict(metadata=path.relative_to(ROOT).as_posix(),recording=recording.relative_to(ROOT).as_posix() if complete else None,recorded_metadata_source_sha256=source,executed_module=module,executed_module_sha256=digest(ROOT/module) if module else None,kernels_sha256=digest(ROOT/'kernels.py') if module=='solver.py' else None,note=note,status=status))
  rows.append(dict(n=c['n'],nu=c['nu'],step=c['dt'],seed=c['seed'],requested_end=c['end'],last_saved_time=q['rows'][-1]['t'],status=status,mean_re_after20_available=float(np.mean([r['re_lambda'] for r in window])) if window else '',mean_energy_after20_available=float(np.mean([r['energy'] for r in window])) if window else '',mean_dissipation_after20_available=float(np.mean([r['dissipation'] for r in window])) if window else '',minimum_sampled_keta=min(r['keta'] for r in q['rows']),max_cfl=max(r['cfl_max'] for r in q['rows']),max_absolute_energy_balance=max(abs(r['balance']) for r in q['rows'])))
 dest=ROOT/'results/screening_regimes.csv'
 with dest.open('w',newline='',encoding='utf-8') as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 save('results/screening_provenance.json',provenance)
 path=ROOT/'data/screen/n80_nu0.02_dt0.025_seed300003.json';q=read(str(path.relative_to(ROOT)));initial=q['rows'][0];bound=initial['peak_speed']*.025/(6/80)
 save('data/screen/events.json',dict(events=[dict(case=path.name,event='Initial-step CFL guard rejection',observed_initial_speed=initial['peak_speed'],initial_stage_cfl_lower_bound=bound,limit=1.,completed=False,note='The tool execution returned a RuntimeError at the guard. This is a numerical screening rejection, not evidence of a physical instability.')],coarse_grid_note='The 48-cubed pilot violates the later declared cutoff-times-Kolmogorov requirement of 1.5 and is excluded from confirmatory predictions.'))
 # Long trajectory comparisons are descriptive single-seed checks, distinct from paired restarts.
 comparisons=[]
 for nu,seed,a,b,da,db in [(.03,300001,56,64,.025,.025),(.025,300000,64,80,.025,.0125),(.02,300003,64,80,.025,.0125)]:
  tags=[f'n{n}_nu{nu}_dt{dt}_seed{seed}' for n,dt in [(a,da),(b,db)]]
  paths=[ROOT/'data/screen'/(tag+'.json') for tag in tags]
  if not all((ROOT/'data/screen'/(tag+'.npz')).exists() for tag in tags):continue
  qs=[read(str(path.relative_to(ROOT))) for path in paths];means=[{key:float(np.mean([r[key] for r in q['rows'] if 20<=r['t']<45])) for key in ['energy','dissipation','re_lambda','flux_k4']} for q in qs]
  comparisons.append(dict(nu=nu,seed=seed,coarse_n=a,fine_n=b,coarse_step=da,fine_step=db,coarse_means=means[0],fine_means=means[1],relative_changes={key:means[1][key]/means[0][key]-1 for key in means[0]},interpretation='Descriptive time-window means from one matched seed; not independent ensemble uncertainty. Space and time are jointly changed when steps differ.'))
 save('results/long_refinement_comparisons.json',comparisons)
 print(json.dumps(dict(cases=len(rows),completed=sum(r['status']=='completed' for r in rows),unfinished=[r for r in rows if r['status'].startswith('unfinished')],long_comparisons=len(comparisons))))
if __name__=='__main__':main()
