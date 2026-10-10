"""Describe paired interventions; do not refit or reinterpret prediction models."""
from pathlib import Path
import json
import sys
import csv
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from features import future_target,index_at,neighbors,wrap
from simulate import sha,save_json


def same_source(path,digest):
    import hashlib
    raw=Path(path).read_bytes(); lf=raw.replace(b'\r\n',b'\n')
    return digest in {hashlib.sha256(b).hexdigest() for b in [raw,lf,lf.replace(b'\n',b'\r\n')]}


def compare(a,b,cutoff=None):
    delta=b-a
    out={'same_label_RMSE':float(np.sqrt(np.mean(delta**2))),'same_label_MAE':float(np.mean(abs(delta))),
         'mean_reversed_minus_normal':float(delta.mean()),'normal_mean':float(a.mean()),'reversed_mean':float(b.mean()),
         'correlation':float(np.corrcoef(a,b)[0,1]) if a.std()>0 and b.std()>0 else None}
    if cutoff is not None:
        aa=a>=cutoff; bb=b>=cutoff; union=(aa|bb).sum()
        out.update({'normal_events':int(aa.sum()),'reversed_events':int(bb.sum()),'intersection':int((aa&bb).sum()),
                    'event_disagreement_fraction':float(np.mean(aa!=bb)),'event_Jaccard':float((aa&bb).sum()/union) if union else None})
    return out


def bootstrap(values):
    a=np.asarray(values); rng=np.random.default_rng(2006)
    means=a[rng.integers(0,len(a),(2000,len(a)))].mean(1)
    return {'mean':float(a.mean()),'descriptive_ci95':np.quantile(means,[.025,.975]).tolist(),'independent_runs':len(a)}


def main():
    p=json.loads((HERE/'protocol.json').read_text()); execution=json.loads((HERE/'execution.json').read_text())
    if not execution['complete'] or len(execution['jobs'])!=9: raise ValueError('Incomplete authorized pilot')
    cutoff=p['event_threshold']
    if cutoff!=json.loads((HERE.parent/'results'/'primary.json').read_text())['train_event_cutoff']:
        raise ValueError('Changed original event cutoff')
    result={'protocol':p,'runs':{},'grid_comparisons':{},'half_step_comparison':{},'summary':{},
            'source_sha256':{'analyze.py':sha(HERE/'analyze.py'),'run_pairs.py':sha(HERE/'run_pairs.py'),'protocol.json':sha(HERE/'protocol.json')}}
    data={}; rows=[]; targets={}
    for job in execution['jobs']:
        name=job['name']; path=HERE/'recorded-data'/(name+'.npz'); meta=json.loads(path.with_suffix('.json').read_text())
        if sha(path)!=meta['sha256'] or sha(path.with_name(name+'_branch-state.npz'))!=meta['checkpoint_sha256']:
            raise ValueError('Saved result hash mismatch')
        for rel,digest in meta['sources'].items():
            if not same_source(HERE.parent.parent/rel,digest): raise ValueError('Executed source changed')
        with np.load(path,allow_pickle=False) as z: z=dict(z)
        if not all(np.isfinite(v).all() for v in z.values()): raise ValueError('Nonfinite output')
        data[name]=z; rec={'seed':meta['seed'],'n':meta['n'],'dt':meta['dt'],'data_sha256':meta['sha256'],
                         'checkpoint_sha256':meta['checkpoint_sha256'],'replay_match':meta['replay_match'],'checks':meta['checks'],
                         'diagnostics':meta['diagnostics'],'horizons':{}}
        for horizon in p['horizons']:
            i=index_at(z['normal_time'],horizon)
            a=future_target(z['normal_tangent'][0],z['normal_tangent'][i],horizon)
            b=future_target(z['reversed_tangent'][0],z['reversed_tangent'][i],horizon)
            rec['horizons'][str(horizon)]=compare(a,b,cutoff if horizon==p['primary_horizon'] else None)
            targets[name+'_'+str(horizon)+'_normal']=a; targets[name+'_'+str(horizon)+'_reversed']=b
        idx=neighbors(z['normal_positions'][0],p['neighbors'])
        d0=np.linalg.norm(wrap(z['normal_positions'][0,idx]-z['normal_positions'][0,:,None]),axis=-1)
        geo=[]
        for label in ['normal','reversed']:
            x=z[label+'_positions'][-1]; distance=np.linalg.norm(wrap(x[idx]-x[:,None]),axis=-1)
            geo.append(np.log(np.maximum(distance,1e-12)/np.maximum(d0,1e-12)).mean(axis=1)/p['duration'])
        rec['finite_neighbor_log_rate']=compare(*geo)
        rec['retrace_position_RMS']=[float(np.sqrt(np.mean(np.sum(wrap(x-y)**2,axis=1))))
                                    for x,y in zip(z['reversed_positions'],z['recorded_past_positions'])]
        result['runs'][name]=rec
        rows.append({'name':name,'seed':meta['seed'],'grid':meta['n'],'dt':meta['dt'],**rec['horizons']['0.2'],
                     'retrace_position_RMS_at_0.2':rec['retrace_position_RMS'][-1]})
    for seed in p['seeds']:
        coarse=f's{seed}_n38_dt0.005'; fine=f's{seed}_n56_dt0.005'; checks={}
        for label in ['normal','reversed']:
            a=targets[coarse+'_0.2_'+label]; b=targets[fine+'_0.2_'+label]
            checks[label]=compare(a,b,cutoff)
        dc=targets[coarse+'_0.2_reversed']-targets[coarse+'_0.2_normal']
        df=targets[fine+'_0.2_reversed']-targets[fine+'_0.2_normal']
        checks['branch_contrast_grid_RMSE']=float(np.sqrt(np.mean((dc-df)**2)))
        checks['branch_contrast_grid_error_over_fine_contrast_RMS']=float(np.linalg.norm(dc-df)/np.linalg.norm(df))
        checks['starting_marker_grid_RMS']=float(np.sqrt(np.mean(np.sum(wrap(data[coarse]['normal_positions'][0]-data[fine]['normal_positions'][0])**2,axis=1))))
        result['grid_comparisons'][str(seed)]=checks
    for label in ['normal','reversed']:
        a=targets['s1024_n56_dt0.005_0.2_'+label]; b=targets['s1024_n56_dt0.0025_0.2_'+label]
        result['half_step_comparison'][label]=compare(a,b,cutoff)
    for n in p['grids']:
        selected=[result['runs'][f's{s}_n{n}_dt0.005']['horizons']['0.2'] for s in p['seeds']]
        result['summary'][str(n)]={metric:bootstrap([x[metric] for x in selected]) for metric in
                                  ['same_label_RMSE','event_disagreement_fraction','event_Jaccard','mean_reversed_minus_normal','correlation']}
    save_json(HERE/'results.json',result); np.savez_compressed(HERE/'targets.npz',**targets)
    with (HERE/'per-run.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    name='s1024_n56_dt0.005'; z=data[name]
    a=targets[name+'_0.2_normal']; b=targets[name+'_0.2_reversed']
    ax=axes[0,0]; ax.scatter(a,b,s=9,alpha=.5,color='#306fa0')
    lo=min(a.min(),b.min()); hi=max(a.max(),b.max()); ax.plot([lo,hi],[lo,hi],color='gray',lw=1)
    ax.set(xlabel='Normal branch FTLE',ylabel='Reversed-velocity branch FTLE',title='Same 512 starting labels; seed 1024, grid 56')
    ax=axes[0,1]
    for s in p['seeds']:
        ax.plot(p['grids'],[result['runs'][f's{s}_n{n}_dt0.005']['horizons']['0.2']['same_label_RMSE'] for n in p['grids']],'-o',label=str(s))
    ax.set(xlabel='Grid points per direction',ylabel='RMS difference between branch FTLEs',title='Physical branch difference at elapsed 0.2'); ax.set_xticks(p['grids']); ax.legend(title='Seed')
    ax=axes[1,0]; x=wrap(z['normal_positions'][0]); ea=a>=cutoff; eb=b>=cutoff
    for mask,color,label in [(~(ea|eb),'#bbbbbb','Neither'),(ea&eb,'#7b3294','Both'),(ea&~eb,'#2166ac','Normal only'),(eb&~ea,'#d6604d','Reversed only')]:
        ax.scatter(x[mask,0],x[mask,1],s=12,color=color,label=label,alpha=.7)
    ax.set(xlabel='Starting x',ylabel='Starting y',title='Future high-stretch labels: initial xy projection'); ax.legend(fontsize=8)
    ax=axes[1,1]
    for s in p['seeds']:
        name=f's{s}_n56_dt0.005'; ax.plot(data[name]['normal_time'],result['runs'][name]['retrace_position_RMS'],label=str(s))
    ax.set(xlabel='Elapsed time after velocity reversal',ylabel='RMS distance from recorded backward path',title='Reversing velocity does not exactly retrace history'); ax.legend(title='Seed')
    fig.suptitle('Forward physical intervention: same starting geometry, opposite velocity\n4 initial conditions; positive viscosity; two-grid pilot, not established convergence',fontsize=12)
    fig.savefig(HERE/'comparison.png',dpi=150); fig.savefig(HERE/'comparison.svg'); plt.close(fig)
    print(json.dumps({'summary':result['summary'],'grid_comparisons':result['grid_comparisons'],'half_step':result['half_step_comparison']},indent=2))


if __name__=='__main__':main()
