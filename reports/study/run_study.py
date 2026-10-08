"""Resumable, ordered execution of the user-requested three-case matrix."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
import datetime as dtlib
import gc
import json
import math
import msvcrt
import shutil
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from protocol import ROOT,GAMMA,GRIDS,CASES,METHODS,CFL_LIMIT,jobs,protocol,source_hashes,atomic_json
from solver import Solver,seed_markers,validate
from initial_design import inspect

def now():
    return dtlib.datetime.now(dtlib.timezone.utc).isoformat()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def save_checkpoint(path,h,markers,accum,initial,step,t,job,hashes):
    temp=path.with_suffix('.tmp')
    with temp.open('wb') as f:
        np.savez(f,h=h,markers=markers,accum=accum,initial=initial,step=step,t=t,
                 job_json=json.dumps(job),source_hashes=json.dumps(hashes))
    temp.replace(path)

def measure(s,h,t,markers,accum,initial):
    row=s.observe(h,t,markers,accum,initial)
    row['max_speed']=float(np.sqrt(np.sum(s.real(h)**2,axis=0)).max())
    row['scaled_native_divergence']=row['divergence_native']*s.dx/max(row['max_speed'],1e-300)
    row['energy_balance_relative']=row['energy_balance_residual']/initial[0]
    row['enstrophy_balance_relative']=row['native_enstrophy_balance_residual']/initial[1]
    row['resolution_flags']=[]
    if row['global_width']['minimum_chord_cells'] is None or row['global_width']['minimum_chord_cells']<6:
        row['resolution_flags'].append('global peak has fewer than six cells across a measured half-peak chord')
    if row['high_band_enstrophy_fraction']>.01:
        row['resolution_flags'].append('more than one percent of physical enstrophy is in the highest retained band')
    # JSON serialization rejects NaN/Infinity in every diagnostic, including nested ones.
    json.dumps(row,allow_nan=False)
    return row

def verify_row(row):
    if row['scaled_native_divergence']>1e-10:
        raise RuntimeError('Native discrete divergence gate failed')
    if abs(row['energy_balance_relative'])>.005:
        raise RuntimeError('Energy balance residual exceeded 0.5 percent of initial energy')

def preflight():
    expected=source_hashes()
    path=ROOT/'preflight.json'
    if path.exists():
        prev=read(path)
        if prev.get('status')=='passed' and prev.get('source_hashes')==expected:
            print('Reusing passed preflight.',flush=True);return prev
    result={'status':'running','started':now(),'source_hashes':expected,
            'solver_validation':validate(),'initial_rows':[]}
    # This is a material design control: global metrics must not compare two
    # fields that merely translate the same axially uniform central tube.
    probe=Solver(48)
    ha,_,_=probe.initial('aligned',GAMMA);hb,_,_=probe.initial('compressive',GAMMA)
    translated_error=float(np.max(abs(probe.real(ha)-np.roll(probe.real(hb),24,axis=3))))
    rates=[]
    for field in (ha,hb):
        rhs,_,_=probe.rhs(field)
        rates.append(float(probe.real(probe.curl(rhs))[(2,24,24,24)]))
    assert translated_error>1 and rates[0]>0 and rates[1]<0
    result['distinct_initial_interactions']={'half_box_translation_max_velocity_difference':translated_error,
        'aligned_origin_domega_z_dt':rates[0],'compressive_origin_domega_z_dt':rates[1]}
    del probe,ha,hb;gc.collect()
    atomic_json(ROOT/'protocol.json',protocol());atomic_json(path,result)
    for case in CASES:
        for n in GRIDS:
            geometry=inspect(n,case,GAMMA)
            assert max(abs(x) for x in geometry['mean_velocity'])<1e-10
            assert geometry['core_radius_cells']>=4
            if case=='aligned':
                assert geometry['parallel_stretch']>3.99 and geometry['alignment_to_largest_eigenvector']>.999
            if case=='compressive':assert geometry['parallel_stretch']< -3.99
            if case=='exodus':assert min(geometry['outer_mean_outward_speeds'])>0
            for method in METHODS:
                s=Solver(n,method);h,pos,change=s.initial(case,GAMMA)
                initial=np.array([.5*s.inner(h,h),.5*s.inner(s.curl(h),s.curl(h))])
                row=measure(s,h,0.,seed_markers(pos),np.zeros(7),initial)
                verify_row(row)
                assert row['core_width']['minimum_chord_cells']>=8
                if case=='aligned':assert row['origin_axial_strain']>3.98
                if case=='compressive':assert row['origin_axial_strain']< -3.98
                result['initial_rows'].append({'case':case,'n':n,'method':method,
                    'initial_projection_relative_l2_change':change,'geometry':geometry,'diagnostics':row})
                atomic_json(path,result)
                print(f'preflight {case} {method} {n}: W={row["Wmax"]:.9g}, width={row["core_width"]["minimum_chord_cells"]:.2f} cells',flush=True)
                del s,h,row;gc.collect()
    result['status']='passed';result['finished']=now();atomic_json(path,result)
    return result

def run_job(job):
    folder=ROOT/'runs'/job['id'];folder.mkdir(parents=True,exist_ok=True)
    result_path=folder/'result.json';checkpoint=folder/'checkpoint.npz'
    hashes=source_hashes()
    s=Solver(job['n'],job['method'],job['nu'],workers=4)
    started=time.perf_counter();record=None
    try:
        if result_path.exists():
            record=read(result_path)
            if record['source_hashes']!=hashes or record['job']!=job:
                raise RuntimeError('Source or settings changed; refusing to overwrite previous computation')
            if record['status']=='complete':return {'id':job['id'],'status':'complete','reused':True}
            if record['status']=='failed':raise RuntimeError('Saved numerical failure requires investigation before resuming')
        if checkpoint.exists():
            with np.load(checkpoint,allow_pickle=False) as a:
                assert json.loads(str(a['source_hashes']))==hashes
                assert json.loads(str(a['job_json']))==job
                h=a['h'];markers=a['markers'];accum=a['accum'];initial=a['initial']
                step=int(a['step']);t=float(a['t'])
            if record is None:raise RuntimeError('Checkpoint has no matching observation record')
            # An interrupted write may leave a newer field than the JSON table.
            rows=[r for r in record['rows'] if r['t']<=t+1e-12]
            if not rows or abs(rows[-1]['t']-t)>1e-12:
                rows.append(measure(s,h,t,markers,accum,initial))
            record['rows']=rows
        else:
            h,pos,change=s.initial(job['case'],job['gamma']);markers=seed_markers(pos)
            accum=np.zeros(7);step=0;t=0.
            initial=np.array([.5*s.inner(h,h),.5*s.inner(s.curl(h),s.curl(h))])
            first=measure(s,h,t,markers,accum,initial);verify_row(first)
            record={'status':'running','job':job,'source_hashes':hashes,'started':now(),
                'initial_projection_relative_l2_change':change,'rows':[first],
                'maximum_stage_cfl':0.,'largest_stage_W':first['Wmax'],'elapsed_seconds':0.}
            save_checkpoint(checkpoint,h,markers,accum,initial,step,t,job,hashes)
            shutil.copyfile(checkpoint,folder/'field-t0.00.npz')
        record.update(status='running',pid=os.getpid(),updated=now());atomic_json(result_path,record)
        previous_elapsed=record['elapsed_seconds']
        block=job['steps_per_output'];blocks=round(job['end']/job['output_dt']);total=block*blocks
        if step%block!=0:raise RuntimeError('Resume is not at an observation checkpoint')
        u=s.real(h)
        for stop_block in range(step//block+1,blocks+1):
            for _ in range(block):
                new,added,cfl=s.step(h,job['dt'])
                if not np.isfinite(new).all() or not np.isfinite(added).all() or not math.isfinite(cfl):
                    with (folder/'nonfinite-evidence.npz').open('wb') as f:np.savez(f,h_before=h,h_after=new,step=step)
                    raise RuntimeError('Nonfinite field or diagnostic detected')
                next_u=s.real(new)
                markers=s.move_markers(markers,u,next_u,job['dt'])
                h=new;u=next_u;accum+=added;step+=1;t=step*job['dt']
                record['maximum_stage_cfl']=max(record['maximum_stage_cfl'],cfl)
                if cfl>CFL_LIMIT:
                    save_checkpoint(folder/'cfl-failure.npz',h,markers,accum,initial,step,t,job,hashes)
                    raise RuntimeError(f'Stage CFL {cfl:.6g} exceeds {CFL_LIMIT}')
            t=stop_block*job['output_dt']
            row=measure(s,h,t,markers,accum,initial)
            row['max_stage_cfl_so_far']=record['maximum_stage_cfl']
            record['rows'].append(row)
            record.update(t=t,step=step,total_steps=total,updated=now(),elapsed_seconds=previous_elapsed+time.perf_counter()-started)
            save_checkpoint(checkpoint,h,markers,accum,initial,step,t,job,hashes)
            if stop_block%5==0:shutil.copyfile(checkpoint,folder/f'field-t{t:.2f}.npz')
            atomic_json(result_path,record)
            verify_row(row)
            print(f'{job["id"]}: t={t:.2f} W={row["Wmax"]:.6g} I={row["I"]:.6g} energy_res={row["energy_balance_relative"]:.2e}',flush=True)
        record.update(status='complete',finished=now());atomic_json(result_path,record)
        return {'id':job['id'],'status':'complete','elapsed_seconds':record['elapsed_seconds']}
    except Exception as e:
        if record is None:record={'job':job,'source_hashes':hashes,'rows':[]}
        record.update(status='failed',updated=now(),error=str(e),traceback=traceback.format_exc())
        atomic_json(result_path,record)
        print(f'FAILED {job["id"]}: {e}',flush=True)
        return {'id':job['id'],'status':'failed','error':str(e)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--preflight-only',action='store_true')
    p.add_argument('--case',choices=CASES);p.add_argument('--n',type=int);p.add_argument('--method',choices=METHODS)
    p.add_argument('--half',action='store_true');args=p.parse_args()
    # Windows kernel lock automatically releases if the process exits or crashes.
    lock=(ROOT/'study.lock').open('a+b');lock.seek(0)
    if lock.read(1)==b'':lock.write(b'0');lock.flush()
    lock.seek(0)
    try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    except OSError:raise SystemExit('Another study process already holds the lock; no duplicate started.')
    state_path=ROOT/'task-state.json'
    state={'status':'preflight','pid':os.getpid(),'started':now(),'completed':[],'active':[],'source_hashes':source_hashes()}
    atomic_json(state_path,state)
    try:
        preflight()
        if args.preflight_only:
            state.update(status='preflight_complete',updated=now());atomic_json(state_path,state);return
        selected=list(jobs(args.case))
        if args.n:selected=[j for j in selected if j['n']==args.n and j['half']==args.half]
        if args.method:selected=[j for j in selected if j['method']==args.method]
        state.update(status='running',planned_jobs=[j['id'] for j in selected]);atomic_json(state_path,state)
        with ProcessPoolExecutor(max_workers=2) as pool:
            for case in CASES:
                case_jobs=[j for j in selected if j['case']==case]
                # Adjacent jobs share grid and timestep, and differ only in spatial method.
                for start in range(0,len(case_jobs),2):
                    pair=case_jobs[start:start+2]
                    state.update(active=[j['id'] for j in pair],updated=now());atomic_json(state_path,state)
                    futures=[pool.submit(run_job,j) for j in pair]
                    outcomes=[f.result() for f in futures]
                    state['completed'].extend(outcomes)
                    if any(r['status']=='failed' for r in outcomes):
                        state.update(status='stopped_for_numerical_failure',active=[],updated=now());atomic_json(state_path,state);return
                if case_jobs:
                    state.update(last_completed_case=case,updated=now());atomic_json(state_path,state)
        state.update(status='complete',active=[],finished=now());atomic_json(state_path,state)
    except Exception:
        state.update(status='failed',error=traceback.format_exc(),updated=now());atomic_json(state_path,state);raise
    finally:
        lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()

if __name__=='__main__':main()
