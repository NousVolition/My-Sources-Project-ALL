"""Reporting only; retain all numerical screens whether they pass or fail."""
from pathlib import Path
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import json,hashlib,csv
from datetime import datetime, timezone
import numpy as np
from scipy.fft import fft2,ifft2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
def read(n):return json.loads((R/n).read_text(encoding='utf-8'))
def write(n,x):(R/n).write_text(json.dumps(x,indent=2,default=lambda y:y.tolist() if isinstance(y,np.ndarray) else y.item()),encoding='utf-8')
handoff=json.loads((R/'handoff.json').read_text(encoding='utf-8-sig'))
if handoff.get('analysis'):
    handoff.setdefault('analysis_attempts',[]).append(handoff['analysis'])
handoff['analysis']=dict(pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),stdout='analysis.log',stderr='analysis.err')
handoff['status']='analysis_running'
write('handoff.json',handoff)
p=read('protocol.json');cs=read('results/cases.json');run_summary=read('results/summary.json')
for name,sha in read('freeze.json')['files'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==sha,name
labels=run_summary['batch_labels'];data={a:read('results/'+a+'.json') for a in labels}
fields={a:np.load(R/'results'/f'{a}.npz')['fields'] for a in labels}
look={(c['mode'],c['weight'],c['pulse']):c['id'] for c in cs}
fine='n48-dt0.005';f=fields[fine];records=data[fine]['records'];time=np.array([r['t'] for r in records])
checks=[]
def check(name,value,limit,details=None):
    item=dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit))
    if details is not None:item['details']=details
    checks.append(item)
def embed(f,n):
    old=f.shape[-1];modes=np.rint(np.fft.fftfreq(old)*old).astype(int)%n
    out=np.zeros(f.shape[:-2]+(n,n),dtype=complex);out[...,modes[:,None],modes[None,:]]=fft2(f,norm='forward')
    return ifft2(out,norm='forward').real
field_checks=[]
for a,b,limit in [('n24-dt0.01','n32-dt0.01',p['screens']['state_grid_rms']),('n32-dt0.01','n48-dt0.01',p['screens']['state_grid_rms']),('n48-dt0.01',fine,p['screens']['state_time_rms'])]:
    fa=fields[a];fb=fields[b]
    diff=embed(fa,fb.shape[-1])-fb if fa.shape[-1]!=fb.shape[-1] else fa-fb
    rms=np.sqrt(np.mean(diff**2,axis=(-2,-1)))
    detail=dict(pair=[a,b],component_maxima=rms.max(axis=(0,1)),per_case_component_maxima=rms.max(axis=0))
    check('full_band_fields_'+a+'_vs_'+b,rms.max(),limit,detail);field_checks.append(detail)
    delta=np.abs(np.array([x['response_rms'] for x in data[a]['records']])-np.array([x['response_rms'] for x in data[b]['records']]))
    denom=np.maximum(np.max(np.array([x['response_rms'] for x in data[b]['records']]),axis=0),p['screens']['response_absolute_floor'])
    relative=delta.max(axis=0)/denom
    rlimit=p['screens']['response_grid_relative'] if fa.shape[-1]!=fb.shape[-1] else p['screens']['response_time_relative']
    check('paired_response_curves_'+a+'_vs_'+b,relative.max(),rlimit,dict(per_case_component_relative=relative))
for label,d in data.items():
    b=d['bounds']
    check('energy_balance_'+label,b['max_energy_residual'],p['hard_limits']['energy_relative'])
    check('incompressibility_'+label,b['max_divergence'],p['hard_limits']['divergence'])
    check('momentum_'+label,b['max_mean_drift'],p['hard_limits']['mean_drift'])
    check('CFL_'+label,b['max_CFL'],p['hard_limits']['CFL'])
    check('guard_bounds_'+label,max(0,-b['guard_min'],b['guard_max']-1),p['hard_limits']['guard_tolerance'])
swapped=np.load(R/'results/relabel.npz')['fields'];ids=[c['id'] for c in cs if c['mode']=='directional']
reference=fields['n32-dt0.01'][:,ids][:,:,[0,1,3,2,5,4]]
check('participant_relabel',np.max(np.abs(reference-swapped)),p['screens']['relabel_max'])
check('independent_curl_DOP853',max(x['max_field_error'] for x in read('results/independent.json')),p['independent']['limit'])
shared=[c['id'] for c in cs if c['mode']=='shared']
check('shared_guard_equal',np.max(np.abs(f[:,shared,4]-f[:,shared,5])),1e-11)
shared_flow=[]
for pulse in p['pulses']:
    equal=look['shared','equal',pulse]
    shared_flow.extend(float(np.max(np.abs(f[:,look['shared',w,pulse],:2]-f[:,equal,:2]))) for w in p['weights'])
check('shared_weight_flow_null',max(shared_flow),p['screens']['shared_weight_flow_null'])
off=[c['id'] for c in cs if c['mode']=='no_feedback']
check('feedback_disabled_flow_null',np.max(np.abs(f[:,off,:2]-f[:,off[0]:off[0]+1,:2])),p['screens']['feedback_disabled_flow_null'])

metrics=[]
def values(label,key):return np.array([r[key] for r in data[label]['records']])
for c in cs:
    if c['pulse']=='none':continue
    i=c['id'];receiver=3 if c['pulse']=='A' else 2;direction='A_to_B' if c['pulse']=='A' else 'B_to_A'
    rr=values(fine,'response_rms')[:,i,receiver];vr=values(fine,'flow_response_rms')[:,i]
    gr=values(fine,'response_rms')[:,i,4:6]
    get_peak=lambda label,key,component=None: float((values(label,key)[:,i] if component is None else values(label,key)[:,i,component]).max())
    peak=float(rr.max());flow=float(vr.max())
    metrics.append(dict(id=i,mode=c['mode'],weight=c['weight'],pulse=c['pulse'],direction=direction,
        receiver_peak_rms=peak,receiver_peak_time=float(time[rr.argmax()]),flow_peak_rms=flow,flow_peak_time=float(time[vr.argmax()]),
        guard_peak_rms=float(gr.max()),receiver_end_rms=float(rr[-1]),flow_end_rms=float(vr[-1]),
        energy_end=records[-1]['energy'][i],gap_end=records[-1]['gap_rms'][i],
        receiver_grid_peak_difference=abs(get_peak('n48-dt0.01','response_rms',receiver)-get_peak('n32-dt0.01','response_rms',receiver)),
        receiver_time_peak_difference=abs(peak-get_peak('n48-dt0.01','response_rms',receiver)),
        flow_grid_peak_difference=abs(get_peak('n48-dt0.01','flow_response_rms')-get_peak('n32-dt0.01','flow_response_rms')),
        flow_time_peak_difference=abs(flow-get_peak('n48-dt0.01','flow_response_rms'))))
with (R/'results/measurements.csv').open('w',newline='',encoding='utf-8') as out:
    w=csv.DictWriter(out,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
metric={r['id']:r for r in metrics};contrasts=[]
for r in metrics:
    if r['mode'] not in ('shared','directional'):continue
    equal=metric[look[r['mode'],'equal',r['pulse']]];frozen=metric[look['frozen',r['weight'],r['pulse']]]
    contrasts.append(dict(id=r['id'],mode=r['mode'],weight=r['weight'],direction=r['direction'],
        receiver_weight_change=r['receiver_peak_rms']-equal['receiver_peak_rms'],
        receiver_weight_change_relative=r['receiver_peak_rms']/equal['receiver_peak_rms']-1,
        flow_weight_change=r['flow_peak_rms']-equal['flow_peak_rms'],
        receiver_adaptation_relative=r['receiver_peak_rms']/frozen['receiver_peak_rms']-1,
        receiver_numerical_sensitivity=max(r['receiver_grid_peak_difference'],r['receiver_time_peak_difference'],equal['receiver_grid_peak_difference'],equal['receiver_time_peak_difference']),
        flow_numerical_sensitivity=max(r['flow_grid_peak_difference'],r['flow_time_peak_difference'],equal['flow_grid_peak_difference'],equal['flow_time_peak_difference'])))
ratios={mode:metric[look[mode,'equal','A']]['receiver_peak_rms']/metric[look[mode,'equal','B']]['receiver_peak_rms'] for mode in ('shared','directional')}
controls=[]
for r in metrics:
    if r['mode']!='directional':continue
    for mode in ('no_feedback','no_strain','frozen'):
        other=metric[look[mode,r['weight'],r['pulse']]]
        controls.append(dict(weight=r['weight'],direction=r['direction'],control=mode,
            receiver_difference=r['receiver_peak_rms']-other['receiver_peak_rms'],
            flow_difference=r['flow_peak_rms']-other['flow_peak_rms'],
            guard_difference=r['guard_peak_rms']-other['guard_peak_rms']))
verification=dict(passed=all(x['passed'] for x in checks),checks=checks,preflight=read('preflight.json'),frozen_hashes_matched=True)
write('results/verification.json',verification)
summary=dict(status='analyzed',verification_passed=verification['passed'],counts=p['expected'],metrics=metrics,contrasts=contrasts,
             equal_weight_response_ratios=ratios,control_contrasts=controls,
             elapsed_seconds=run_summary['elapsed_seconds'],initial_energy=records[0]['energy'][0],
             numerical_failures=[x for x in checks if not x['passed']],bounds={k:v['bounds'] for k,v in data.items()})
write('results/analysis.json',summary)
F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fafbf9','axes.facecolor':'#fafbf9'})
colors=['#246b8e','#c95b36'];weights=list(p['weights']);weight_labels=['Equal settings','B → A stronger','A → B stronger']
fig,axes=plt.subplots(2,2,figsize=(11,7),sharey='row')
for col,mode in enumerate(('shared','directional')):
    for j,pulse in enumerate(('A','B')):
        rr=[metric[look[mode,w,pulse]] for w in weights]
        for row,key in enumerate(('receiver_peak_rms','flow_peak_rms')):
            axes[row,col].bar(np.arange(3)+(j-.5)*.32,[r[key] for r in rr],width=.3,color=colors[j],label=['A → B pulse','B → A pulse'][j])
    axes[0,col].set_title(['Shared evolving guard','Draft directional guards'][col],weight='bold')
    for ax in axes[:,col]:ax.set_xticks(range(3),weight_labels,rotation=10);ax.grid(axis='y',alpha=.15);ax.legend(frameon=False,fontsize=8)
axes[0,0].set_ylabel('Other SIM peak RMS response');axes[1,0].set_ylabel('Flow peak RMS response')
fig.suptitle('A disturbance changes the SIMS and the flow through the active guard')
fig.text(.5,.01,'One smooth spatial starting pattern; each response subtracts its own matched undisturbed trajectory.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,.95])
for ext in ('png','svg'):fig.savefig(F/f'comparison.{ext}',dpi=155)
plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,6.5))
for w,ls in zip(weights,('-','--',':')):
    i=look['directional',w,'B'];label=weight_labels[weights.index(w)]
    axes[0,0].plot(time,values(fine,'response_rms')[:,i,2],ls,label=label)
    axes[0,1].plot(time,values(fine,'flow_response_rms')[:,i],ls,label=label)
i=look['directional','equal','B'];base=look['directional','equal','none']
for j in range(2):
    axes[1,0].plot(time,values(fine,'mean_guard')[:,i,j],label=['A receives B','B receives A'][j],color=colors[j])
    axes[1,0].plot(time,values(fine,'mean_guard')[:,base,j],':',color=colors[j],alpha=.7)
for mode in ('directional','no_feedback','no_strain','frozen'):
    axes[1,1].plot(time,values(fine,'energy')[:,look[mode,'equal','B']],label=mode.replace('_',' '))
for ax in axes.ravel():ax.axvspan(.5,1,color='#cad5bf',alpha=.4);ax.set_xlabel('Model time');ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
axes[0,0].set_ylabel('A response to pulse on B');axes[0,1].set_ylabel('Flow response to pulse on B');axes[1,0].set_ylabel('Spatial mean guard transmission');axes[1,1].set_ylabel('Kinetic energy')
fig.suptitle('The guard remains in the fluid interaction before, during and after the pulse')
fig.text(.5,.01,'Directional model, pulse on B. Dotted guard curves are the matched undisturbed paths. Energy curves include different controls.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,.95])
for ext in ('png','svg'):fig.savefig(F/f'timeline.{ext}',dpi=155)
plt.close(fig)
# Show spatial changes at a stated saved time, with shared scales across settings.
snapshot=2;selected=[look['directional',w,'B'] for w in weights];baseids=[look['directional',w,'none'] for w in weights]
delta=f[snapshot,selected]-f[snapshot,baseids]
arrays=[delta[:,2],np.sqrt(delta[:,0]**2+delta[:,1]**2),delta[:,4]]
fig,axes=plt.subplots(3,3,figsize=(10,8),layout='constrained')
for row,arr in enumerate(arrays):
    lim=float(np.max(np.abs(arr)))
    for col in range(3):
        im=axes[row,col].imshow(arr[col].T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],cmap='viridis' if row==1 else 'RdBu_r',vmin=0 if row==1 else -lim,vmax=lim)
        axes[row,col].set_xticks([]);axes[row,col].set_yticks([])
        if row==0:axes[row,col].set_title(weight_labels[col],fontsize=10)
    axes[row,0].set_ylabel(['A state change','Flow speed of difference','A receives B guard change'][row]);fig.colorbar(im,ax=axes[row,:],shrink=.75)
fig.suptitle('Pulse on B at model time 1.5: changes from each matched undisturbed run\nIllustrative spatial fields; the detailed guard-field refinement screen failed.',fontsize=11)
for ext in ('png','svg'):fig.savefig(F/f'spatial.{ext}',dpi=155)
plt.close(fig)
print(json.dumps(dict(passed=verification['passed'],failures=[{k:v for k,v in c.items() if k!='details'} for c in summary['numerical_failures']],ratios=ratios),indent=2))
handoff['status']='analysis_completed_pending_review'
handoff['analysis']['completed_utc']=datetime.now(timezone.utc).isoformat()
write('handoff.json',handoff)
