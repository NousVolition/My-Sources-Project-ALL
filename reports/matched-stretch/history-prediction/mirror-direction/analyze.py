"""Four-cell and anchored-neighborhood comparisons; no predictive refitting."""
from pathlib import Path
import sys
import json
import csv
import importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from run_matrix import R,cloud_target,sha,save_json,wrap,source_matches
from features import future_target,index_at
spec=importlib.util.spec_from_file_location('prior_velocity_analysis',HERE.parent/'velocity-reversal'/'analyze.py')
prior=importlib.util.module_from_spec(spec); spec.loader.exec_module(prior)


def factorial(cells):
    op,om,mp,mm=[cells[k] for k in ['original_normal','original_reversed','mirror_normal','mirror_reversed']]
    contrasts={'arrangement_average':.5*((mp-op)+(mm-om)),
               'direction_average':.5*((om-op)+(mm-mp)),
               'interaction_difference_in_differences':(mm-mp)-(om-op)}
    return {k:{'mean':float(v.mean()),'RMS':float(np.sqrt(np.mean(v*v))),'MAE':float(abs(v).mean())} for k,v in contrasts.items()}


def main():
    p=json.loads((HERE/'protocol.json').read_text()); execution=json.loads((HERE/'execution.json').read_text())
    if not execution['complete'] or len(execution['jobs'])!=9: raise ValueError('Incomplete authorized matrix')
    cutoff=p['event_threshold']; result={'protocol':p,'runs':{},'summary':{},'grid_comparisons':{},'half_step_comparison':{},
                                       'sources':{'analyze.py':sha(HERE/'analyze.py')}}
    all_targets={}; raw={}; table=[]
    for job in execution['jobs']:
        name=job['name']; path=HERE/'recorded-data'/(name+'.npz'); meta=json.loads(path.with_suffix('.json').read_text())
        if sha(path)!=meta['sha256']: raise ValueError('New data changed')
        for rel,digest in meta['sources'].items():
            if not source_matches(HERE.parent.parent/rel,digest): raise ValueError('Executed source changed')
        oldpath=HERE.parent/'velocity-reversal'/'recorded-data'/(name+'.npz')
        if sha(oldpath)!=meta['reused_archive_sha256']: raise ValueError('Reused data changed')
        with np.load(path,allow_pickle=False) as z: new=dict(z)
        with np.load(oldpath,allow_pickle=False) as z: old=dict(z)
        if not all(np.isfinite(v).all() for v in new.values()): raise ValueError('Nonfinite data')
        raw[name]=(new,old)
        rec={'seed':meta['seed'],'n':meta['n'],'dt':meta['dt'],'data_sha256':meta['sha256'],
             'checks':meta['checks'],'diagnostics':meta['diagnostics'],'horizons':{}}
        for horizon in p['horizons']:
            i=index_at(new['normal_time'],horizon); globals_={}; clouds={}; pairs={}
            for label in ['normal','reversed']:
                globals_['original_'+label]=future_target(old[label+'_tangent'][0],old[label+'_tangent'][i],horizon)
                globals_['mirror_'+label]=future_target(new[label+'_global_tangent'][0],new[label+'_global_tangent'][i],horizon)
                centers=old[label+'_positions'][i]; idx=new['neighbor_indices']
                r1=wrap(centers[idx]-centers[:,None]); rm=wrap(new[label+'_local_neighbor_positions'][i]-centers[:,None])
                r0=new['initial_offsets']; r0m=new['initial_local_mirror_offsets']
                for arrangement,initial,final in [('original',r0,r1),('mirror',r0m,rm)]:
                    key=arrangement+'_'+label
                    clouds[key]=cloud_target(initial,final,horizon)
                    pairs[key]=np.log(np.maximum(np.linalg.norm(final,axis=-1),1e-12)/np.maximum(np.linalg.norm(initial,axis=-1),1e-12)).mean(axis=1)/horizon
                if label+'_full_mirror_tangent' in new:
                    target=future_target(new[label+'_full_mirror_tangent'][0],new[label+'_full_mirror_tangent'][i],horizon)
                    rec['checks'][label+'_full_mirror_target_error_'+str(horizon)]=float(abs(target-globals_['original_'+label]).max())
            comparisons={}
            for key,a,b in [('direction_original','original_normal','original_reversed'),('direction_global_mirror','mirror_normal','mirror_reversed'),
                            ('global_mirror_normal','original_normal','mirror_normal'),('global_mirror_reversed','original_reversed','mirror_reversed')]:
                comparisons[key]=prior.compare(globals_[a],globals_[b],cutoff if horizon==.2 else None)
            local={key:prior.compare(clouds['original_'+key],clouds['mirror_'+key]) for key in ['normal','reversed']}
            rec['horizons'][str(horizon)]={'global_point_FTLE':comparisons,'anchored_cloud_mirror':local,
                                          'anchored_cloud_factorial':factorial(clouds),'anchored_pair_distance_factorial':factorial(pairs)}
            for family,values in [('global',globals_),('cloud',clouds),('pair',pairs)]:
                for key,value in values.items(): all_targets[name+'_'+str(horizon)+'_'+family+'_'+key]=value
        result['runs'][name]=rec
        table.append({'name':name,'seed':meta['seed'],'grid':meta['n'],'dt':meta['dt'],
                      **{k+'_point_RMSE':v['same_label_RMSE'] for k,v in rec['horizons']['0.2']['global_point_FTLE'].items()},
                      **{k+'_cloud_RMS':v['RMS'] for k,v in rec['horizons']['0.2']['anchored_cloud_factorial'].items()}})
    for n in p['grids']:
        selected=[result['runs'][f's{s}_n{n}_dt0.005']['horizons']['0.2'] for s in p['seeds']]
        out={'global_point_FTLE':{},'anchored_cloud_factorial':{},'anchored_pair_distance_factorial':{}}
        for key in selected[0]['global_point_FTLE']:
            out['global_point_FTLE'][key]={m:prior.bootstrap([r['global_point_FTLE'][key][m] for r in selected]) for m in
                                         ['same_label_RMSE','event_disagreement_fraction','event_Jaccard','correlation']}
        for family in ['anchored_cloud_factorial','anchored_pair_distance_factorial']:
            out[family]={k:prior.bootstrap([r[family][k]['RMS'] for r in selected]) for k in selected[0][family]}
        result['summary'][str(n)]=out
    for seed in p['seeds']:
        result['grid_comparisons'][str(seed)]={}
        for family in ['global','cloud','pair']:
            result['grid_comparisons'][str(seed)][family]={}
            for cell in ['original_normal','original_reversed','mirror_normal','mirror_reversed']:
                a=all_targets[f's{seed}_n38_dt0.005_0.2_{family}_{cell}']; b=all_targets[f's{seed}_n56_dt0.005_0.2_{family}_{cell}']
                result['grid_comparisons'][str(seed)][family][cell]=prior.compare(a,b,cutoff if family=='global' else None)
    for family in ['global','cloud','pair']:
        result['half_step_comparison'][family]={}
        for cell in ['original_normal','original_reversed','mirror_normal','mirror_reversed']:
            a=all_targets[f's1024_n56_dt0.005_0.2_{family}_{cell}']; b=all_targets[f's1024_n56_dt0.0025_0.2_{family}_{cell}']
            result['half_step_comparison'][family][cell]=prior.compare(a,b,cutoff if family=='global' else None)
    save_json(HERE/'results.json',result); np.savez_compressed(HERE/'targets.npz',**all_targets)
    with (HERE/'per-run.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0])); w.writeheader(); w.writerows(table)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained'); name='s1024_n56_dt0.005'
    new,old=raw[name]
    for ax,arrangement in zip(axes[0],['original','mirror']):
        a=all_targets[name+'_0.2_global_'+arrangement+'_normal']; b=all_targets[name+'_0.2_global_'+arrangement+'_reversed']
        x=wrap(old['normal_positions'][0] if arrangement=='original' else new['normal_global_positions'][0])
        for mask,color,label in [(~((a>=cutoff)|(b>=cutoff)),'#bbbbbb','Neither'),((a>=cutoff)&(b>=cutoff),'#7b3294','Both'),((a>=cutoff)&(b<cutoff),'#2166ac','Normal only'),((a<cutoff)&(b>=cutoff),'#d6604d','Reversed only')]:
            ax.scatter(x[mask,0],x[mask,1],s=13,alpha=.7,color=color,label=label)
        ax.set(xlabel='Initial x',ylabel='Initial y',title=arrangement.capitalize()+' starting cloud; seed 1024, grid 56'); ax.legend(fontsize=8)
    ax=axes[1,0]
    a=all_targets[name+'_0.2_cloud_original_normal']; b=all_targets[name+'_0.2_cloud_mirror_normal']
    ax.scatter(a,b,s=10,alpha=.5,color='#2166ac');lo=min(a.min(),b.min());hi=max(a.max(),b.max());ax.plot([lo,hi],[lo,hi],color='gray',lw=1)
    ax.set(xlabel='Original local-cloud deformation rate',ylabel='Mirrored local-cloud deformation rate',title='Same centers, changed finite neighbor layout')
    ax=axes[1,1]; keys=['arrangement_average','direction_average','interaction_difference_in_differences']
    labels=['Arrangement contrast','Direction contrast','Interaction (difference in differences)']
    for i,key in enumerate(keys):
        vals=[result['runs'][f's{s}_n56_dt0.005']['horizons']['0.2']['anchored_cloud_factorial'][key]['RMS'] for s in p['seeds']]
        ax.scatter(vals,[i]*len(vals),s=30,color='#2166ac');ax.scatter([np.mean(vals)],[i],marker='D',s=45,color='#d6604d')
    ax.set_yticks(range(3),labels);ax.invert_yaxis();ax.set(xlabel='RMS of predeclared contrast',title='Finite-cloud contrasts: four independent starts');ax.grid(axis='x',alpha=.2)
    fig.suptitle('Mirror arrangement × reverse physical velocity\nPoint markers are passive; initial xy projections; local contrasts have different scales',fontsize=12)
    fig.savefig(HERE/'comparison.png',dpi=150);fig.savefig(HERE/'comparison.svg');plt.close(fig)
    print(json.dumps({'summary':result['summary'],'half_step':result['half_step_comparison']},indent=2))


if __name__=='__main__':main()
