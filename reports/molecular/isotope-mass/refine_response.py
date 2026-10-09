"""Exploratory follow-up prompted by failed 0.5-ps response convergence.
Preserves the originally declared primary endpoint and its failed checks.
"""
from pathlib import Path
import concurrent.futures, json
import numpy as np
from simulate import system,response,dump
from analyze import response_metrics,summarize

def run(arg):
    c,heavy,dt,j=arg;z=np.load(f'data/main_c{c:02d}_v0_A.npz');q=z['snapshot_q_nm'][1];v=z['snapshot_v_nmps'][1];s,_=system(9,heavy)
    raw,base,e,lag=response(s,q,v,dt,j)
    zz=dict(response_positions_nm=raw[None],response_baseline_nm=base[None],response_energy_kjmol=e[None],snapshot_q_nm=q[None],impulse_da_nmps=j)
    name=f'c{c:02d}_h{heavy}_dt{dt}_j{j}';np.savez_compressed(Path('refinement')/(name+'.npz'),**zz)
    met,score,mats,scores=response_metrics(zz)
    return dict(config=c,heavy=heavy,dt=dt,j=j,metrics=met,matrices=mats)

def main():
    Path('refinement').mkdir(exist_ok=True)
    # 8 configs, three placements, two finer steps and two smaller impulses.
    jobs=[(c,h,dt,j) for c in range(8) for h in [-1,0,8] for dt,j in [(.00025,.0125),(.000125,.0125),(.000125,.00625)]]
    with concurrent.futures.ProcessPoolExecutor(4) as pool:results=list(pool.map(run,jobs))
    dump(Path('refinement/results.json'),results);errors=[];effects={}
    for c in range(8):
        for h in [-1,0,8]:
            arr=[r for r in results if r['config']==c and r['heavy']==h]
            for name,a,b in [('step_025_to_0125',arr[0],arr[1]),('impulse_0125_to_00625',arr[1],arr[2])]:
                errors.append(dict(config=c,heavy=h,check=name,relative_by_lag=[float(np.linalg.norm(a['matrices'][l]-b['matrices'][l])/max(np.linalg.norm(b['matrices'][l]),1e-30)) for l in range(4)]))
    for metric in ['response_002','response_010','response_offdiag']:
        diffs=[]
        for c in range(8):
            r=[x for x in results if x['config']==c and x['dt']==.000125 and x['j']==.00625];a=next(x for x in r if x['heavy']==-1)
            diffs.append(np.mean([x['metrics'][metric]-a['metrics'][metric] for x in r if x['heavy']>=0]))
        effects[metric]=summarize(diffs)
    dump(Path('refinement/summary.json'),dict(exploratory_after_failed_check=True,assays=len(results),errors=errors,matched_state_effects=effects))
    print('Completed 72 finer-step/smaller-impulse conditional response assays.',flush=True)

if __name__=='__main__':main()
