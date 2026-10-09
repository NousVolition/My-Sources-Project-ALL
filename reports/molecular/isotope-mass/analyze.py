"""Reanalyze raw trajectories. Uncertainty is over initial-configuration clusters."""
from pathlib import Path
import argparse, concurrent.futures, csv, json
import numpy as np
from scipy.stats import rankdata
from simulate import dump, R_GAS

def geometry(x, cutoff=.35, angle=150):
    o=x[:,:,0]; h=x[:,:,1:]
    oo=o[:,:,None,:]-o[:,None,:,:]; d=np.linalg.norm(oo,axis=-1)
    n=o.shape[1];d[:,np.arange(n),np.arange(n)]=np.inf
    ha=o[:,None,:,None,:]-h[:,:,None,:,:];oh=h-o[:,:,None,:]
    cosine=np.einsum('fdhc,fdahc->fdah',oh,ha)/(np.linalg.norm(oh,axis=-1)[:,:,None,:]*np.linalg.norm(ha,axis=-1))
    directed=(d<cutoff)&(cosine>=np.cos(np.deg2rad(180-angle))).any(axis=-1)
    hb=directed|directed.transpose(0,2,1)
    return d,hb,d<cutoff

def episode_data(contact, sample):
    """Episodes beginning after observation starts; initial episodes are left censored."""
    a=np.asarray(contact,bool).reshape(len(contact),-1); lengths=[];events=[];left=0
    for col in a.T:
        starts=np.flatnonzero(np.diff(col.astype(int),prepend=0)==1)
        stops=np.flatnonzero(np.diff(col.astype(int),append=0)==-1)+1
        for start,stop in zip(starts,stops):
            if start==0:left+=1;continue
            event=stop<len(col)
            lengths.append(((stop if event else stop-1)-start)*sample);events.append(event)
    return np.array(lengths),np.array(events,bool),left

def km_rmst(lengths,events,tau=2.):
    if not len(lengths):return float('nan')
    s=1.;area=0.;last=0.
    for t in np.unique(lengths[lengths<=tau]):
        area+=s*(t-last);risk=np.sum(lengths>=t);deaths=np.sum((lengths==t)&events)
        s*=1-deaths/risk;last=t
    return area+s*(tau-last)

def lifetimes(contact,sample):
    lengths,events,left=episode_data(contact,sample)
    return dict(rmst_ps=km_rmst(lengths,events),episodes=len(lengths),right_censored=int((~events).sum()),left_excluded=left,
                mean_completed_ps=float(np.mean(lengths[events])) if events.any() else float('nan'))

def corr(a,b):
    a=rankdata(a);b=rankdata(b)
    return float(np.corrcoef(a,b)[0,1]) if np.std(a)*np.std(b)>0 else float('nan')

def response_metrics(z):
    raw=z['response_positions_nm'];j=float(z['impulse_da_nmps']);n=raw.shape[-2]
    derivative=(raw[:,:,:,1]-raw[:,:,:,0])/(2*j) # snapshot,source,axis,lag,target,coord
    g=np.sqrt(np.mean(np.sum(derivative**2,axis=-1),axis=2))
    for i in range(n):g[:,i,:,i]=0
    scores=g.sum(axis=-1)/(n-1)
    final=g[:,:,-1];off=final.sum(axis=-1)
    participation=off**2/np.maximum((final**2).sum(axis=-1),1e-40)
    q=z['snapshot_q_nm'][:,::3];dist=np.linalg.norm(q[:,:,None]-q[:,None,:],axis=-1)
    far=(final*(dist>.35)).sum(axis=-1)/np.maximum(off,1e-40)
    metrics=dict(response_offdiag=float(scores[:,:,-1].mean()),response_002=float(scores[:,:,1].mean()),response_010=float(scores[:,:,2].mean()),response_participation=float(participation.mean()),response_far_share=float(far.mean()),
                 response_top_repeat=float(np.argmax(scores[0,:,-1])==np.argmax(scores[-1,:,-1])),
                 response_rank_stability=corr(scores[0,:,-1],scores[-1,:,-1]),
                 response_nve_energy_range=float(np.max(np.ptp(z['response_energy_kjmol'],axis=-1))))
    return metrics,scores[:,:,-1].mean(axis=0),g.mean(axis=0).transpose(1,0,2),scores

def one(path):
    meta=json.loads(path.read_text());z=np.load(path.with_suffix('.npz'));x=z['positions_nm'];dt=meta['sample_ps'];n=meta['n'];nf=len(x);half=nf//2
    d,hb,contact=geometry(x);degree=hb.sum(axis=-1);coord=contact.sum(axis=-1);nn=d.argmin(axis=-1)
    duration=(nf-1)*dt;switch=(nn[1:]!=nn[:-1]).sum(axis=0)/duration
    upper=np.triu_indices(n,1);life=lifetimes(hb[:,upper[0],upper[1]],dt);res=lifetimes(contact[:,upper[0],upper[1]],dt)
    lag=round(.5/dt);sums=np.concatenate([np.zeros((1,n,n),int),np.cumsum(hb,axis=0)],axis=0)
    continuously=(sums[lag+1:]-sums[:-(lag+1)])==lag+1
    survival=continuously.sum()/max(1,hb[:-lag].sum())
    meanadj=hb.mean(axis=0);ev,vec=np.linalg.eigh(meanadj);central=np.abs(vec[:,-1]);central/=central.sum()
    top=degree==degree.max(axis=1,keepdims=True);topshare=(top/top.sum(axis=1,keepdims=True)).mean(axis=0)
    ene=z['energies_kjmol'];temp=2*ene[:,1]/(R_GAS*6*n)
    row={k:meta[k] for k in ['suite','config','velocity','heavy','label','n','production_ps','equil_ps','dt_ps','gamma_ps','radius_nm']}
    row.update(stem=path.stem,hb_degree=float(degree.mean()),coordination=float(coord.mean()),degree_sd=float(degree.std(axis=1).mean()),nn_exchange=float(switch.mean()),
               hb_rmst=life['rmst_ps'],residence_rmst=res['rmst_ps'],hb_survival_05=float(survival),
               hb_completed_mean=life['mean_completed_ps'],hb_episodes=life['episodes'],hb_right=life['right_censored'],hb_left=life['left_excluded'],
               residence_episodes=res['episodes'],residence_right=res['right_censored'],residence_left=res['left_excluded'],
               rank_stability=corr(degree[:half].mean(axis=0),degree[half:].mean(axis=0)),top_repeat=int(np.argmax(degree[:half].mean(axis=0))==np.argmax(degree[half:].mean(axis=0))),
               hb_first=float(degree[:half].mean()),hb_second=float(degree[half:].mean()),temperature=float(temp.mean()),
               temperature_first=float(temp[:half].mean()),temperature_second=float(temp[half:].mean()),
               potential_first=float(ene[:half,0].mean()/n),potential_second=float(ene[half:,0].mean()/n),
               wall_fraction=float((np.linalg.norm(x[:,:,0],axis=-1)>meta['radius_nm']).mean()))
    mol=[]
    for i in range(n):
        remote=hb.copy();remote[:,i,:]=False;remote[:,:,i]=False
        mol.append(dict(stem=path.stem,site=i,hb_degree=float(degree[:,i].mean()),coordination=float(coord[:,i].mean()),nn_exchange=float(switch[i]),centrality=float(central[i]),
                        topshare=float(topshare[i]),remote_hb_degree=float(remote.sum()/(nf*(n-1)))))
    resp=None
    if 'response_positions_nm' in z:
        a,sc,mats,score=response_metrics(z);row.update(a)
        for i in range(n):mol[i]['response_out']=float(sc[i])
        resp=(path.stem,mats.tolist(),score.tolist())
    # Duration diagnostic uses nested prefix from the very same trajectory (never independent samples).
    row['nn_first'] = float((nn[1:half]!=nn[:half-1]).sum()/((half-1)*dt*n))
    row['nn_stride2'] = float((nn[2::2]!=nn[:-2:2]).sum()/(duration*n))
    row['hb_rmst_stride2'] = lifetimes(hb[::2,upper[0],upper[1]],2*dt)['rmst_ps']
    # Label invariance: compare stored physical rows after restoring original indices.
    label=None
    if meta['label']:
        a=np.load(path.with_name(path.stem.replace('_label','_A')+'.npz'))
        err=float(np.max(np.abs(x[:,np.argsort(z['permutation'])]-a['positions_nm'])))
        label=dict(stem=path.stem,max_position_error_nm=err)
    return row,mol,resp,label

def write_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def summarize(values, baseline=None):
    values=np.array(values,float);values=values[np.isfinite(values)];n=len(values)
    if n<2:return dict(n=n,mean=float(values.mean()) if n else None)
    rng=np.random.default_rng(733102);boots=values[rng.integers(n,size=(10000,n))].mean(axis=1)
    if n<=16:
        signs=2*((np.arange(2**n)[:,None]>>np.arange(n))&1)-1
        p=float(np.mean(np.abs((signs*values).mean(axis=1))>=abs(values.mean())-1e-14))
    else:p=None
    sd=values.std(ddof=1)
    r=dict(n=n,mean=float(values.mean()),ci95=np.percentile(boots,[2.5,97.5]).tolist(),sd=float(sd),dz=float(values.mean()/sd) if sd>0 else None,signflip_p=p,cluster_values=values.tolist())
    if baseline is not None:r['baseline']=float(np.mean(baseline));r['percent']=float(100*values.mean()/np.mean(baseline)) if np.mean(baseline) else None
    return r

def inference(rows,molecules):
    mol={(x['stem'],x['site']):x for x in molecules};out={};primary=['hb_degree','nn_exchange','response_offdiag']
    keys=['hb_degree','coordination','degree_sd','nn_exchange','hb_rmst','residence_rmst','hb_survival_05','response_offdiag','response_002','response_010','response_participation','response_far_share','rank_stability','nn_stride2','hb_rmst_stride2']
    for suite in sorted(set(x['suite'] for x in rows)):
        subset=[r for r in rows if r['suite']==suite and not r['label']];cfgs=sorted(set(r['config'] for r in subset));effects={};interaction={};local={}
        for metric in keys:
            if metric not in subset[0]:continue
            vals=[];bases=[];ivs=[]
            for c in cfgs:
                r=[x for x in subset if x['config']==c];a={x['velocity']:x for x in r if x['heavy']<0};ds=[x for x in r if x['heavy']>=0]
                vals.append(np.mean([x[metric]-a[x['velocity']][metric] for x in ds]));bases.append(np.mean([x[metric] for x in a.values()]))
                inner=[x[metric] for x in ds if x['heavy']==0];outer=[x[metric] for x in ds if x['heavy']==x['n']-1];ivs.append(np.mean(inner)-np.mean(outer))
            effects[metric]=summarize(vals,bases);interaction[metric]=summarize(ivs)
        for metric in ['hb_degree','coordination','nn_exchange','centrality','topshare','remote_hb_degree','response_out']:
            vals=[];ints=[]
            if metric not in mol[(subset[0]['stem'],0)]:continue
            for c in cfgs:
                r=[x for x in subset if x['config']==c];a={x['velocity']:x for x in r if x['heavy']<0};ds=[x for x in r if x['heavy']>=0]
                deltas=[(x['heavy'],mol[(x['stem'],x['heavy'])][metric]-mol[(a[x['velocity']]['stem'],x['heavy'])][metric]) for x in ds]
                vals.append(np.mean([v for _,v in deltas]));ints.append(np.mean([v for h,v in deltas if h==0])-np.mean([v for h,v in deltas if h==r[0]['n']-1]))
            local[metric]=dict(mass=summarize(vals),mass_position=summarize(ints))
        out[suite]=dict(mass=effects,initial_position_interaction=interaction,tagged=local)
    pvals=sorted([(out['main']['mass'][k]['signflip_p'],k) for k in primary]);last=0
    for rank,(p,k) in enumerate(pvals):
        last=max(last,min(1,p*(len(primary)-rank)));out['main']['mass'][k]['holm_p']=last
    for suite in [s for s in out if s!='main']:
        out[suite]['paired_sensitivity_shift']={}
        for metric in ['hb_degree','nn_exchange','hb_rmst','residence_rmst','response_offdiag']:
            if metric not in out[suite]['mass']:continue
            changes=[]
            cfgs=sorted(set(r['config'] for r in rows if r['suite']==suite))
            for c in cfgs:
                effects=[]
                for which in ['main',suite]:
                    r=[x for x in rows if x['suite']==which and x['config']==c and not x['label'] and x['heavy'] in [-1,0,x['n']-1]]
                    a={x['velocity']:x for x in r if x['heavy']<0}
                    effects.append(np.mean([x[metric]-a[x['velocity']][metric] for x in r if x['heavy']>=0]))
                changes.append(effects[1]-effects[0])
            out[suite]['paired_sensitivity_shift'][metric]=summarize(changes)
    # Absolute stationarity shifts (A only); sensitivity comparisons remain descriptively paired by config.
    out['diagnostics']={}
    for suite in sorted(set(x['suite'] for x in rows)):
        aa=[x for x in rows if x['suite']==suite and x['heavy']<0 and not x['label']];cfg=sorted(set(x['config'] for x in aa))
        out['diagnostics'][suite]={k:summarize([np.mean([x[end]-x[start] for x in aa if x['config']==c]) for c in cfg]) for k,start,end in [('hb_second_minus_first','hb_first','hb_second'),('potential_second_minus_first','potential_first','potential_second'),('temperature_second_minus_first','temperature_first','temperature_second')]}
        out['diagnostics'][suite]['temperature_K']=summarize([np.mean([x['temperature'] for x in aa if x['config']==c]) for c in cfg])
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=Path('data'));p.add_argument('--out',type=Path,default=Path('analysis'));p.add_argument('--workers',type=int,default=4);args=p.parse_args();args.out.mkdir(exist_ok=True)
    paths=sorted(args.data.glob('*_c*_v*_*.json'));rows=[];mol=[];responses={};labels=[]
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        for i,(r,m,re,l) in enumerate(pool.map(one,paths)):
            rows.append(r);mol.extend(m)
            if re:responses[re[0]]=dict(matrices=re[1],scores=re[2])
            if l:labels.append(l)
            if i%50==0:print(f'analyzed {i+1}/{len(paths)}',flush=True)
    write_csv(args.out/'trajectories.csv',rows);write_csv(args.out/'molecules.csv',mol)
    dump(args.out/'response_matrices.json',responses);dump(args.out/'label_controls.json',labels)
    dump(args.out/'statistics.json',inference(rows,mol))
    dump(args.out/'inventory.json',dict(completed=len(rows),counts={s:sum(r['suite']==s for r in rows) for s in sorted(set(r['suite'] for r in rows))},label_controls=len(labels),max_label_error_nm=max(x['max_position_error_nm'] for x in labels),frames=sum(round(r['production_ps']/.05)+1 for r in rows)))

if __name__=='__main__':main()
