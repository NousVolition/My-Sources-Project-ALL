"""Measure which passive labels have the fastest mean neighbor separation."""
from pathlib import Path
import os,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def wrap(x):return (x+3)%6-3
def scores(p,v,edges):
    r=wrap(p[edges[:,1]]-p[edges[:,0]])
    assert np.all(np.linalg.norm(r,axis=1)>1e-14)
    assert np.all(abs(abs(r)-3)>1e-12), 'Periodic cut locus: derivative undefined'
    rates=np.sum(r*(v[edges[:,1]]-v[edges[:,0]]),axis=1)/np.linalg.norm(r,axis=1)
    out=np.zeros(250);degree=np.zeros(250)
    for col in (0,1):np.add.at(out,edges[:,col],rates);np.add.at(degree,edges[:,col],1)
    return out/degree
def leading(s,offset=0):
    peak=float(s.max());tol=1e-10*max(1.,float(abs(s).max()))
    return [int(i+offset) for i in np.flatnonzero(peak-s<=tol)] if peak>tol else []
def main():
    data=json.loads((ROOT/'results.json').read_text());edges=np.array(data['edges']);runs={};maxerr=0.
    for name,record in data['runs'].items():
        path=ROOT/(name+'.npz');assert hashlib.sha256(path.read_bytes()).hexdigest()==record['trajectory_sha256']
        z=np.load(path,allow_pickle=False);rows=[[],[]]
        for k,t in enumerate(z['time']):
            p=z['positions'][k,:250];v=z['velocities'][k,:250];s=scores(p,v,edges)
            # Independent scalar evaluation of the node-wise means.
            check=[]
            for i in range(250):
                neighbors=[int(b if a==i else a) for a,b in edges if a==i or b==i]
                check.append(np.mean([np.dot(wrap(p[j]-p[i]),v[j]-v[i])/np.linalg.norm(wrap(p[j]-p[i])) for j in neighbors]))
            err=float(np.max(abs(s-check)));assert err<1e-10;maxerr=max(maxerr,err)
            for g in (0,1):
                values=s[g*125:(g+1)*125];ordered=np.sort(values)
                rows[g].append({'time':float(t),'labels':leading(values,g*125),'maximum_rate':float(ordered[-1]),'top_two_gap':float(ordered[-1]-ordered[-2])})
        summary=[]
        for group in rows:
            labels=sorted({i for row in group for i in row['labels']});longest=0;holder=None
            for i in labels:
                streak=0
                for row in group:
                    streak=streak+1 if i in row['labels'] else 0
                    if streak>longest:longest=streak;holder=i
            summary.append({'distinct_selected_labels':len(labels),'selected_labels':labels,'changes_between_recorded_sets':sum(a['labels']!=b['labels'] for a,b in zip(group,group[1:])),'longest_consecutive_samples':longest,'longest_streak_label':holder,'samples_with_ties':sum(len(r['labels'])>1 for r in group),'samples_without_positive_separation':sum(not r['labels'] for r in group)})
        runs[name]={'rows':rows,'summary':summary,'source_sha256':record['trajectory_sha256']}
    p=np.array(data['initial_positions']);assert np.all(scores(p,np.ones_like(p),edges)==0)
    offset=np.array([.17,.08,-.11]);v=np.sin(p);assert np.allclose(scores(p,v,edges),scores(p+offset,v,edges),atol=1e-12)
    central=p[:125];dummy=np.vstack([central,central+[1,0,0]])
    assert leading(scores(dummy,-dummy,edges))==[]
    controls={}
    for a,b in [('tracking_medium','reference'),('grid128','grid128_refined'),('reference','fluid_half'),('reference','grid128_refined')]:
        controls[a+' vs '+b]=[]
        for g in (0,1):
            ra=runs[a]['rows'][g];rb=runs[b]['rows'][g];assert [r['time'] for r in ra]==[r['time'] for r in rb]
            controls[a+' vs '+b].append({'same_selected_set':sum(x['labels']==y['labels'] for x,y in zip(ra,rb)),'samples':len(ra),'overlap_nonempty':sum(bool(set(x['labels'])&set(y['labels'])) for x,y in zip(ra,rb))})
    out={'definition':'Mean instantaneous increase in shortest periodic distance to fixed original lattice neighbors. All 125 labels in each cloud are eligible; positive maximum selected, near ties retained. No new fluid force or feedback.','time_window':[0,.3],'tie_tolerance':'1e-10 * max(1, maximum absolute score in cloud)','units':'Model length per model time','groups':['central','periodic_join'],'runs':runs,'controls':controls,'verification':{'status':'passed','independent_score_max_error':maxerr,'uniform_translation_zero_score':True,'translation_invariance':True,'contracting_control_no_positive_role':True},'limits':['Saved times only; switches between samples can be missed.','Consecutive selected samples do not prove continuous residence.','All labels are eligible; not all are observed to lead.','Original grid and initial-data limitations remain.','This is a relative motion statistic; no causal influence or Hénon escape classification is inferred.']}
    (ROOT/'separation-roles.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for g,ax in enumerate(axes):
        rows=runs['reference']['rows'][g];labels=runs['reference']['summary'][g]['selected_labels'];raster=np.zeros((len(labels),len(rows)))
        for k,row in enumerate(rows):
            for i in row['labels']:raster[labels.index(i),k]=1
        ax.imshow(raster,origin='lower',aspect='auto',cmap='Blues',vmin=0,vmax=1,interpolation='nearest');ax.set_yticks(range(len(labels)),labels);ax.set_xticks([0,10,20,30],['0','0.10','0.20','0.30']);ax.set(title=['Central cloud','Periodic-join cloud'][g],xlabel='Model time (saved samples)',ylabel='Fixed marker ID');ax.tick_params(axis='y',labelsize=8)
    fig.suptitle('Fastest mean separation changes between labels | 64-grid reconstruction');fig.savefig(ROOT/'separation-roles.png',dpi=140);plt.close(fig)
    text=['# Which marker separates from its neighbors fastest?','','Each of the 125 labels in a cloud is eligible at every saved time. A fixed identity is kept throughout the track. The selected role depends on the current positions and velocities.','','For marker i, calculate the mean of `(X_j-X_i) dot (u_j-u_i) / |X_j-X_i|` over its original lattice neighbors, using shortest periodic displacements. This is the instantaneous rate of increase in its mean distance to those neighbors. A positive maximum receives the role; near ties are retained. If every score is nonpositive within tolerance, no marker receives it. The fixed lattice gives each label three to six neighbors; their mean gives each label one comparable score.','','![Selected marker IDs](separation-roles.png)','','## Observed changes through time 0.30','','| Grid | Cloud | Different labels selected | Changes between recorded selected sets | Longest consecutive selected samples |','| --- | --- | ---: | ---: | ---: |']
    for run,grid in [('reference',64),('grid128_refined',128)]:
        for g,label in enumerate(['Central','Periodic join']):
            s=runs[run]['summary'][g];text.append(f"| {grid} | {label} | {s['distinct_selected_labels']} | {s['changes_between_recorded_sets']} | {s['longest_consecutive_samples']} |")
    text+=['','There are 31 observations separated by 0.01 model time. Consecutive observations do not establish uninterrupted residence between them. Every marker is eligible, but this finite run does not show that every marker eventually takes the role.','','## Do the controls select the same IDs?','','| Comparison | Central matches / 31 | Join matches / 31 |','| --- | ---: | ---: |']
    for label,c in controls.items():text.append(f"| {label} | {c[0]['same_selected_set']} | {c[1]['same_selected_set']} |")
    text+=['','Matching selected sets is a stringent rank diagnostic and may be sensitive to nearly tied scores. The full records retain the top-two score gap, ties and the score itself. Tracking-step agreement does not establish spatial resolution; grid comparisons also include the known differences in the sampled starting fields.','','This measures relative motion. It does not show that a passive label makes other labels move, that it has a permanent trait, or that it follows an unbounded Hénon orbit. No feedback or force has been added.','','[Full scores and checks](separation-roles.json) · [Rebuild this analysis](separation_roles.py) · [Marker stress test](README.md)','']
    (ROOT/'separation-roles.md').write_text('\n'.join(text),encoding='utf-8');print(json.dumps({'summary':{n:runs[n]['summary'] for n in ['reference','grid128_refined']},'controls':controls,'verification':out['verification']},indent=2))
if __name__=='__main__':main()
