"""Reporting only: remeasure saved arrays, make charts and retain every comparison."""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
p=json.loads((ROOT/'protocol.json').read_text())
v=json.loads((ROOT/'results/verification.json').read_text())
for name,sha in json.loads((ROOT/'freeze.json').read_text())['files'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
cases=json.loads((ROOT/'results/cases.json').read_text())
z=np.load(ROOT/'results/trajectories.npz');t=z['time'];y=z['state']
assert y.shape==(601,540,4) and np.isfinite(y).all()
lookup={(c['pair'],c['mode'],c['weight'],c['disturbance']):c['id'] for c in cases}
records=[]
for c in cases:
    i=c['id'];base=lookup[c['pair'],c['mode'],c['weight'],'none']
    delta=y[:,i]-y[:,base]
    if c['disturbance']=='none':continue
    donor=0 if c['disturbance'].startswith('A') else 1;receiver=1-donor
    a=np.abs(delta[:,receiver]);peak=int(a.argmax())
    records.append(dict(id=i,pair=c['pair'],seed=c['seed'],relabel=c['swapped'],mode=c['mode'],weight=c['weight'],
        pulse=c['disturbance'],direction='A_to_B' if donor==0 else 'B_to_A',
        sign='positive' if c['disturbance'].endswith('positive') else 'negative',
        receiver_peak=float(a[peak]),receiver_peak_time=float(t[peak]),
        receiver_area=float(np.trapezoid(a,t)),receiver_end=float(delta[-1,receiver]),
        own_peak=float(np.abs(delta[:,donor]).max()),guard_response=float(np.abs(delta[:,2:]).max()),
        total_guard_excursion=float(np.abs(y[:,i,2:]-.5).max())))
with (ROOT/'results/measurements.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
def stats(a):
    a=np.array(a,dtype=float)
    return dict(mean=float(a.mean()),minimum=float(a.min()),maximum=float(a.max()),count=len(a))
group=[]
for mode in p['guard_modes']:
    for weight in p['weights']:
        for direction in ('A_to_B','B_to_A'):
            r=[x for x in records if x['mode']==mode and x['weight']==weight and x['direction']==direction]
            group.append(dict(mode=mode,weight=weight,direction=direction,peak=stats([x['receiver_peak'] for x in r]),
                         positive=stats([x['receiver_peak'] for x in r if x['sign']=='positive']),
                         negative=stats([x['receiver_peak'] for x in r if x['sign']=='negative'])))
measure={x['id']:x for x in records}
contrasts=[]
for r in records:
    if r['mode']=='frozen':continue
    frozen=measure[lookup[r['pair'],'frozen',r['weight'],r['pulse']]]
    balanced=measure[lookup[r['pair'],r['mode'],'equal',r['pulse']]]
    contrasts.append(dict(id=r['id'],mode=r['mode'],weight=r['weight'],direction=r['direction'],
                          adaptation_difference=r['receiver_peak']-frozen['receiver_peak'],
                          adaptation_relative=r['receiver_peak']/frozen['receiver_peak']-1,
                          unequal_difference=r['receiver_peak']-balanced['receiver_peak'],
                          unequal_relative=r['receiver_peak']/balanced['receiver_peak']-1))
with (ROOT/'results/contrasts.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(contrasts[0]));w.writeheader();w.writerows(contrasts)
contrasts_summary=[]
for mode in ('shared','directional'):
    for weight in p['weights']:
        for direction in ('A_to_B','B_to_A'):
            r=[x for x in contrasts if (x['mode'],x['weight'],x['direction'])==(mode,weight,direction)]
            contrasts_summary.append(dict(mode=mode,weight=weight,direction=direction,
                adaptation_relative=stats([x['adaptation_relative'] for x in r]),
                unequal_relative=stats([x['unequal_relative'] for x in r])))
gain_stats={}
for mode in ('shared','directional'):
    ids=[c['id'] for c in cases if c['mode']==mode and c['weight']=='equal']
    diff=np.abs(y[:,ids,2]-y[:,ids,3])
    gain_stats[mode]=dict(maximum_equal_setting_gain_difference=float(diff.max()),
        active_guard_excursion=stats([float(np.abs(y[:,i,2:]-.5).max()) for i in ids]))
ratios=[]
for mode in ('shared','directional'):
    a=[]
    for pair in range(12):
        for sign in ('positive','negative'):
            ab=measure[lookup[pair,mode,'equal','A_'+sign]]['receiver_peak']
            ba=measure[lookup[pair,mode,'equal','B_'+sign]]['receiver_peak']
            a.append(ab/ba)
            ratios.append(dict(mode=mode,pair=pair,sign=sign,A_to_B_over_B_to_A=ab/ba))
gap={}
for mode in ('shared','directional'):
    for weight in p['weights']:
        ids=[c['id'] for c in cases if c['mode']==mode and c['weight']==weight and c['disturbance']=='none']
        gap[mode+'/'+weight]=stats(np.abs(y[-1,ids,0]-y[-1,ids,1]))
summary=dict(status='completed',numerical_checks_passed=v['passed'],counts=p['expected'],groups=group,
    contrasts=contrasts_summary,gains=gain_stats,equal_setting_transfer_ratios=ratios,final_undisturbed_gaps=gap,
    active_guard_response_to_pulse=stats([r['guard_response'] for r in records if r['mode']!='frozen']),
    max_end_receiver_response=max(abs(r['receiver_end']) for r in records),
    max_time_refinement_difference=max(c['maximum'] for c in v['checks'] if c['name'].startswith('step')),
    independent_solver_maximum=next(c['maximum'] for c in v['checks'] if c['name']=='independent_DOP853'))
(ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'figure.facecolor':'#fafbf9','axes.facecolor':'#fafbf9','axes.titleweight':'bold'})
figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
colors=['#246b8e','#c95b36'];labels=['Equal settings','B → A stronger','A → B stronger']
fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
for ax,mode,title in zip(axes,('shared','directional'),('Shared evolving guard','Separate directional guards')):
    for j,direction in enumerate(('A_to_B','B_to_A')):
        g=[next(x for x in group if (x['mode'],x['weight'],x['direction'])==(mode,w,direction)) for w in p['weights']]
        means=np.array([x['peak']['mean'] for x in g])
        low=means-np.array([x['peak']['minimum'] for x in g]);high=np.array([x['peak']['maximum'] for x in g])-means
        ax.bar(np.arange(3)+(j-.5)*.32,means,.3,color=colors[j],label=direction.replace('_to_',' → '),alpha=.92)
        ax.errorbar(np.arange(3)+(j-.5)*.32,means,yerr=[low,high],fmt='none',ecolor='#303b41',capsize=3,lw=.8)
    ax.set_xticks(range(3),labels,rotation=12);ax.set_title(title);ax.legend(frameon=False)
axes[0].set_ylabel("Receiver's peak disturbance response")
fig.suptitle('Unequal settings change how much of a disturbance reaches the other SIM')
fig.text(.5,.01,'Means and full ranges over 12 matched starts and ± pulses. Mirrored cases are controls, not independent samples.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.06,1,.94])
for ext in ('png','svg'):fig.savefig(figdir/f'comparison.{ext}',dpi=160)
plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,6.5),sharex=True)
for j,mode in enumerate(('shared','directional')):
    for weight,style in zip(p['weights'],('-','--',':')):
        i=lookup[0,mode,weight,'B_positive'];base=lookup[0,mode,weight,'none']
        axes[0,j].plot(t,y[:,i,0]-y[:,base,0],style,label=labels[list(p['weights']).index(weight)],lw=1.8)
    i=lookup[0,mode,'equal','B_positive']
    axes[1,j].plot(t,y[:,i,2],label='A receives B',color=colors[0])
    axes[1,j].plot(t,y[:,i,3],'--',label='B receives A',color=colors[1])
    axes[0,j].set_title(['Shared evolving guard','Separate directional guards'][j])
    for ax in axes[:,j]:ax.axvspan(2,2.5,color='#c5d3bb',alpha=.4);ax.grid(alpha=.15)
    axes[0,j].legend(frameon=False,fontsize=8);axes[1,j].legend(frameon=False,fontsize=8)
    axes[1,j].set_xlabel('Model time')
axes[0,0].set_ylabel('A response to a pulse on B');axes[1,0].set_ylabel('Guard transmission, equal settings')
fig.suptitle('The guard keeps evolving during and after the disturbance')
fig.text(.5,.01,'Illustration: first seeded pair; all other pairs and both pulse signs are retained in the saved arrays.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,.94])
for ext in ('png','svg'):fig.savefig(figdir/f'timeline.{ext}',dpi=160)
plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4))
for j,mode in enumerate(('shared','directional')):
    r=[x['A_to_B_over_B_to_A'] for x in ratios if x['mode']==mode]
    ax.scatter(np.arange(len(r)),r,color=colors[j],label=['Shared evolving guard','Separate directional guards'][j],s=30,alpha=.75)
ax.axhline(1,color='#424a52',ls='--',lw=1);ax.set_yscale('log');ax.set_xlabel('Matched start / pulse-sign comparison (including relabelings)')
ax.set_ylabel('A → B response / B → A response');ax.set_title('Equal settings do not guarantee equal disturbance responses');ax.legend(frameon=False)
fig.tight_layout()
for ext in ('png','svg'):fig.savefig(figdir/f'within-case.{ext}',dpi=160)
plt.close(fig)
print(json.dumps({'counts':summary['counts'],'numerical_passed':v['passed'],'groups':group,'contrasts':contrasts_summary,
                  'gains':gain_stats,'response_ratio_ranges':{m:stats([r['A_to_B_over_B_to_A'] for r in ratios if r['mode']==m]) for m in ('shared','directional')},
                 'guard_response':summary['active_guard_response_to_pulse'],'refinement':summary['max_time_refinement_difference'],
                 'independent':summary['independent_solver_maximum']},indent=2))

