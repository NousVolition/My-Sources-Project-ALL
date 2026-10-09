"""Conditional mass switch at identical positions AND velocities; not equilibrium.
The instantaneous change changes kinetic energy. Do not pool with equilibrated arms.
"""
from pathlib import Path
import argparse, concurrent.futures
import numpy as np
from simulate import system,response,dump
from analyze import response_metrics,summarize

def run(arg):
    c,heavy,out=arg;out=Path(out);z=np.load(f'data/main_c{c:02d}_v0_A.npz');q=z['snapshot_q_nm'][1];v=z['snapshot_v_nmps'][1]
    s,_=system(9,heavy);raw,base,ene,lags=response(s,q,v)
    zz=dict(response_positions_nm=raw[None],response_baseline_nm=base[None],response_energy_kjmol=ene[None],snapshot_q_nm=q[None],snapshot_v_nmps=v[None],lags_ps=lags,impulse_da_nmps=.05)
    np.savez_compressed(out/f'c{c:02d}_h{heavy}.npz',**zz);met,score,g,scores=response_metrics(zz)
    return dict(config=c,heavy=heavy,metrics=met,scores=score,matrices=g)

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);args=p.parse_args();out=Path('switch');out.mkdir(exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:r=list(pool.map(run,[(c,h,str(out)) for c in range(16) for h in [-1]+list(range(9))]))
    dump(out/'results.json',r);effects=[];tag=[];position=[]
    for c in range(16):
        a=next(x for x in r if x['config']==c and x['heavy']==-1);ds=[x for x in r if x['config']==c and x['heavy']>=0]
        effects.append(np.mean([x['metrics']['response_offdiag']-a['metrics']['response_offdiag'] for x in ds]))
        tag.append(np.mean([x['scores'][x['heavy']]-a['scores'][x['heavy']] for x in ds]))
        position.append(ds[0]['metrics']['response_offdiag']-ds[-1]['metrics']['response_offdiag'])
    # Equal-impulse amplitude sensitivity and same-state half-step sensitivity, all retained.
    sensitivity=[]
    for c in range(8):
        z=np.load(f'data/main_c{c:02d}_v0_A.npz');q=z['snapshot_q_nm'][1];v=z['snapshot_v_nmps'][1];s,_=system()
        ref=next(x for x in r if x['config']==c and x['heavy']==-1)['matrices'][-1]
        for name,dt,j in [('half_impulse',.001,.025),('half_step',.0005,.05)]:
            raw,base,e,l=response(s,q,v,dt,j)
            z2=dict(response_positions_nm=raw[None],response_baseline_nm=base[None],response_energy_kjmol=e[None],snapshot_q_nm=q[None],impulse_da_nmps=j)
            _,_,m,_=response_metrics(z2)
            relative=float(np.linalg.norm(m[-1]-ref)/np.linalg.norm(ref))
            np.savez_compressed(out/f'check_c{c:02d}_{name}.npz',**z2)
            sensitivity.append(dict(config=c,check=name,relative_frobenius_difference=relative))
    dump(out/'summary.json',dict(conditional_switch=True,independent_configs=16,states=160,group_response_change=summarize(effects),tagged_response_change=summarize(tag),inner_outer_difference=summarize(position),sensitivity=sensitivity))
    print('Completed 160 conditional mass-switch assays and 16 numerical response checks.',flush=True)

if __name__=='__main__':main()
