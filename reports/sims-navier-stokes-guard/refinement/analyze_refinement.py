"""Compare saved refinement fields and independently remeasure diagnostics."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,csv
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.fft import fft2,ifft2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def serial(v):return v.tolist() if isinstance(v,np.ndarray) else v.item()
def write(name,data):(R/name).write_text(json.dumps(data,indent=2,default=serial),encoding='utf-8')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
p=read(R/'protocol.json');B=Path(p['refinement']['parent_directory']);h=read(R/'handoff.json')
assert read(R/'results/summary.json')['status']=='completed'
for name,digest in read(R/'freeze.json')['files'].items():assert sha(R/name)==digest,name
for name,digest in p['refinement']['parent_file_hashes'].items():assert sha(B/name)==digest,name
h['analysis']=dict(pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),stdout='analysis.log',stderr='analysis.err');h['status']='analysis_running';write('handoff.json',h)
cs=read(R/'results/cases.json');lookup={(c['mode'],c['weight'],c['pulse']):c['id'] for c in cs}
assert cs==read(B/'results/cases.json'),'Case definitions differ from the original comparison'
components=['u','v','x_A','x_B','g_AB','g_BA'];base=np.array([lookup[c['mode'],c['weight'],'none'] for c in cs])
checks=[]
def check(name,value,limit,details=None):
    row=dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit))
    if details is not None:row['details']=details
    checks.append(row)
def embed(f,n):
    old=f.shape[-1];modes=np.rint(np.fft.fftfreq(old)*old).astype(int)%n
    out=np.zeros(f.shape[:-2]+(n,n),dtype=complex)
    out[...,modes[:,None],modes[None,:]]=fft2(f,axes=(-2,-1),norm='forward',workers=1)
    return ifft2(out,axes=(-2,-1),norm='forward',workers=1).real
old=read(B/'results/verification.json');historical=[]
for c in old['checks']:
    if c['name'].startswith('full_band_fields_n24') or c['name'].startswith('full_band_fields_n32'):
        historical.append(dict(c,scope='Original frozen comparison; remains unchanged'))
labels=['n48-dt0.005','n64-dt0.005','n64-dt0.0025'];data={};fields={};times={}
for label in labels:
    directory=B if label.startswith('n48') else R
    data[label]=read(directory/'results'/f'{label}.json')
    assert data[label]['dt']==float(label.split('-dt')[1])
    assert np.allclose([row['t'] for row in data[label]['records']],np.arange(81)*.05,rtol=0,atol=1e-12)
    z=np.load(directory/'results'/f'{label}.npz');fields[label]=z['fields'];times[label]=z['time']
    assert np.array_equal(times[label],p['time']['field_times'])
    assert fields[label].shape[:3]==(5,45,6)
    assert np.isfinite(fields[label]).all()
spatial_details=None
for a,b,kind,limit in [('n48-dt0.005','n64-dt0.005','spatial',p['screens']['state_grid_rms']),('n64-dt0.005','n64-dt0.0025','temporal',p['screens']['state_time_rms'])]:
    fa,fb=fields[a],fields[b]
    difference=(embed(fa,fb.shape[-1]) if fa.shape[-1]!=fb.shape[-1] else fa)-fb
    rms=np.sqrt(np.mean(difference**2,axis=(-2,-1)))
    ix=np.unravel_index(np.argmax(rms),rms.shape)
    detail=dict(pair=[a,b],component_maxima=rms.max(axis=(0,1)),per_case_component_maxima=rms.max(axis=0),per_saved_time_maxima=rms.max(axis=(1,2)),worst=dict(time=float(times[b][ix[0]]),case=cs[ix[1]],component=components[ix[2]],value=float(rms[ix])))
    check(kind+'_full_band_fields',rms.max(),limit,detail)
    if kind=='spatial':spatial_details=detail
    for key in ('response_rms','flow_response_rms'):
        va=np.array([d[key] for d in data[a]['records']]);vb=np.array([d[key] for d in data[b]['records']])
        denominator=np.maximum(vb.max(axis=0),p['screens']['response_absolute_floor'])
        check(kind+'_'+key+'_curves',np.max(np.max(abs(va-vb),axis=0)/denominator),p['screens']['response_grid_relative' if kind=='spatial' else 'response_time_relative'])
    del difference

# Re-measure every saved new field without calling model diagnostics.
saved_checks=[]
for label in labels[1:]:
    d=data[label];f=fields[label];n=d['n'];ztime=times[label]
    assert d['n']==64 and len(d['records'])==81
    rec=[next(row for row in d['records'] if abs(row['t']-t)<1e-12) for t in ztime]
    x,y=np.meshgrid(np.arange(n)*2*np.pi/n,np.arange(n)*2*np.pi/n,indexing='ij')
    initial=np.stack([np.sin(x)*np.cos(y)-.25*np.sin(2*x+y)-.3*np.cos(x-2*y),-np.cos(x)*np.sin(y)+.5*np.sin(2*x+y)-.15*np.cos(x-2*y),-.55+.15*np.cos(x)+.1*np.sin(y),.35+.1*np.sin(x)-.15*np.cos(y),np.full_like(x,.5),np.full_like(x,.5)])
    check(label+'_analytical_initial_fields',np.max(abs(f[0]-initial)),1e-12)
    delta=f-f[:,base]
    actual=dict(energy=.5*np.mean(f[:,:,0]**2+f[:,:,1]**2,axis=(-2,-1)),mean_velocity=f[:,:,:2].mean(axis=(-2,-1)),mean_guard=f[:,:,4:6].mean(axis=(-2,-1)),gap_rms=np.sqrt(np.mean((f[:,:,3]-f[:,:,2])**2,axis=(-2,-1))),response_rms=np.sqrt(np.mean(delta**2,axis=(-2,-1))),flow_response_rms=np.sqrt(np.mean(np.sum(delta[:,:,:2]**2,axis=2),axis=(-2,-1))))
    k=np.rint(np.fft.fftfreq(n)*n).astype(int);hat=fft2(f,axes=(-2,-1),norm='forward',workers=1)
    dx=ifft2(1j*k[:,None]*hat,axes=(-2,-1),norm='forward',workers=1).real
    dy=ifft2(1j*k[None,:]*hat,axes=(-2,-1),norm='forward',workers=1).real
    strain=dx[:,:,0]**2+dy[:,:,1]**2+.5*(dy[:,:,0]+dx[:,:,1])**2
    nu=.02+.04*(f[:,:,4]+f[:,:,5])
    for c in cs:
        if c['mode']=='no_feedback':nu[:,c['id']]=.06
    actual.update(dissipation=np.mean(2*nu*strain,axis=(-2,-1)),enstrophy=.5*np.mean((dx[:,:,1]-dy[:,:,0])**2,axis=(-2,-1)),mean_viscosity=np.mean(nu,axis=(-2,-1)))
    for key,value in actual.items():
        error=float(np.max(abs(value-np.array([row[key] for row in rec]))));check(label+'_saved_'+key,error,1e-12);saved_checks.append(dict(batch=label,quantity=key,max_abs_error=error))
    check(label+'_saved_divergence',np.max(abs(dx[:,:,0]+dy[:,:,1])),p['hard_limits']['divergence'])
    b=d['bounds']
    for key,limitkey in [('max_CFL','CFL'),('max_divergence','divergence'),('max_mean_drift','mean_drift'),('max_energy_residual','energy_relative')]:check(label+'_'+key,b[key],p['hard_limits'][limitkey])
    check(label+'_guard_range',max(0,-b['guard_min'],b['guard_max']-1),p['hard_limits']['guard_tolerance'])
    energy=np.array([row['energy'] for row in d['records']]);loss=np.array([row['loss'] for row in d['records']])
    check(label+'_energy_monotonic',max(0,np.diff(energy,axis=0).max()),1e-12)
    check(label+'_budget_remeasurement',np.max(abs((energy-energy[0]+loss)/energy[0])),p['hard_limits']['energy_relative'])
    check(label+'_batch_seconds',d['elapsed_seconds'],p['resources']['max_batch_seconds'])
    del delta,hat,dx,dy,strain,actual

fine=labels[-1];f=fields[fine];records=data[fine]['records'];t=np.array([row['t'] for row in records])
def values(label,key):return np.array([row[key] for row in data[label]['records']])
metrics=[]
for c in cs:
    if c['pulse']=='none':continue
    i=c['id'];receiver=3 if c['pulse']=='A' else 2
    rr=values(fine,'response_rms')[:,i,receiver];vr=values(fine,'flow_response_rms')[:,i]
    def peak(label,key,component=None):
        v=values(label,key)[:,i];return float((v if component is None else v[:,component]).max())
    metrics.append(dict(id=i,mode=c['mode'],weight=c['weight'],pulse=c['pulse'],direction='A_to_B' if c['pulse']=='A' else 'B_to_A',receiver_peak_rms=float(rr.max()),receiver_peak_time=float(t[rr.argmax()]),flow_peak_rms=float(vr.max()),flow_peak_time=float(t[vr.argmax()]),guard_peak_rms=float(values(fine,'response_rms')[:,i,4:6].max()),receiver_end_rms=float(rr[-1]),flow_end_rms=float(vr[-1]),
        receiver_spatial_sensitivity=abs(peak(labels[1],'response_rms',receiver)-peak(labels[0],'response_rms',receiver)),receiver_temporal_sensitivity=abs(peak(labels[2],'response_rms',receiver)-peak(labels[1],'response_rms',receiver)),flow_spatial_sensitivity=abs(peak(labels[1],'flow_response_rms')-peak(labels[0],'flow_response_rms')),flow_temporal_sensitivity=abs(peak(labels[2],'flow_response_rms')-peak(labels[1],'flow_response_rms'))))
with (R/'results/measurements.csv').open('w',encoding='utf-8',newline='') as out:
    w=csv.DictWriter(out,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
byid={m['id']:m for m in metrics};contrasts=[]
for m in metrics:
    if m['mode'] not in ('shared','directional'):continue
    equal=byid[lookup[m['mode'],'equal',m['pulse']]];frozen=byid[lookup['frozen',m['weight'],m['pulse']]]
    contrasts.append(dict(id=m['id'],mode=m['mode'],weight=m['weight'],direction=m['direction'],receiver_relative_to_equal=m['receiver_peak_rms']/equal['receiver_peak_rms']-1,flow_relative_to_equal=m['flow_peak_rms']/equal['flow_peak_rms']-1,receiver_relative_to_frozen=m['receiver_peak_rms']/frozen['receiver_peak_rms']-1))
shared=[c['id'] for c in cs if c['mode']=='shared'];off=[c['id'] for c in cs if c['mode']=='no_feedback']
check('shared_guards_equal',np.max(abs(f[:,shared,4]-f[:,shared,5])),1e-11)
check('shared_weight_fluid_null',max(float(np.max(abs(f[:,lookup['shared',w,pulse],:2]-f[:,lookup['shared','equal',pulse],:2]))) for w in p['weights'] for pulse in p['pulses']),p['screens']['shared_weight_flow_null'])
check('feedback_disabled_fluid_null',np.max(abs(f[:,off,:2]-f[:,off[0]:off[0]+1,:2])),p['screens']['feedback_disabled_flow_null'])
check('total_runtime',read(R/'results/summary.json')['elapsed_seconds'],p['resources']['max_total_seconds'])
check('total_storage_GiB',sum(z.stat().st_size for z in R.rglob('*') if z.is_file())/2**30,p['resources']['max_storage_GiB'])
verification=dict(new_checks_passed=all(c['passed'] for c in checks),new_checks=checks,historical_checks=historical,source_and_parent_hashes_matched=True,saved_remeasurements=saved_checks,reused_independent_checks=read(B/'results/independent.json'),reused_relabel_check=next(c for c in old['checks'] if c['name']=='participant_relabel'))
write('results/verification.json',verification)
summary=dict(status='analyzed',new_checks_passed=verification['new_checks_passed'],new_failures=[c for c in checks if not c['passed']],historical_failures_preserved=historical,spatial=spatial_details,spatial_reduction_vs_previous_pair=1-spatial_details['worst']['value']/historical[-1]['value'],metrics=metrics,contrasts=contrasts,directional_equal_response_ratio_B_to_A_over_A_to_B=byid[lookup['directional','equal','B']]['receiver_peak_rms']/byid[lookup['directional','equal','A']]['receiver_peak_rms'],bounds={label:data[label]['bounds'] for label in labels[1:]})
write('results/analysis.json',summary)

F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fafbf9','axes.facecolor':'#fafbf9'})
matrix=np.array([c['details']['component_maxima'] for c in historical]+[spatial_details['component_maxima']])
fig,ax=plt.subplots(figsize=(9,5.6))
colors=['#78899b','#8c9b79','#285e81','#994e24','#577d45','#a22f49']
for j,c in enumerate(components):ax.plot(range(3),matrix[:,j],'-o',label=c,color=colors[j],linewidth=1.7)
ax.axhline(p['screens']['state_grid_rms'],ls='--',color='#222',label='Unchanged tolerance: 0.002')
ax.set_yscale('log');ax.set_xticks(range(3),['24² → 32²\nstep 0.01','32² → 48²\nstep 0.01','48² → 64²\nstep 0.005']);ax.set_ylabel('Largest spatial RMS field difference');ax.grid(alpha=.16,which='both');ax.legend(frameon=False,ncol=3,fontsize=9)
ax.set_title('Does a finer grid reduce the spatial differences?',weight='bold')
fig.text(.5,.012,'All 45 conditions and five saved times; finer-grid modes are retained in each difference.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.05,1,1])
for ext in ('png','svg'):fig.savefig(F/f'refinement.{ext}',dpi=160)
plt.close(fig)

# Inspect the worst new spatial comparison at a common observation grid.
worst=spatial_details['worst'];it=p['time']['field_times'].index(worst['time']);case=worst['case']['id'];component=components.index(worst['component'])
a=embed(fields[labels[0]][it,case,component],64);b=fields[labels[1]][it,case,component];difference=b-a
fig,axes=plt.subplots(1,3,figsize=(11,4),layout='constrained')
lo=min(a.min(),b.min());hi=max(a.max(),b.max())
for ax,arr,title in zip(axes[:2],[a,b],['48² field on a 64² observation grid','64² field']):
    im=ax.imshow(arr.T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],vmin=lo,vmax=hi,cmap='viridis');ax.set_title(title,fontsize=10);ax.set_xlabel('x');ax.set_ylabel('y')
fig.colorbar(im,ax=axes[:2],shrink=.8)
im=axes[2].imshow(difference.T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],cmap='RdBu_r',vmin=-abs(difference).max(),vmax=abs(difference).max());axes[2].set_title('64² minus 48²',fontsize=10);axes[2].set_xlabel('x');axes[2].set_ylabel('y');fig.colorbar(im,ax=axes[2],shrink=.8)
fig.suptitle(f"Largest new field difference: {worst['component']} at model time {worst['time']:g}\n{worst['case']['mode']}, {worst['case']['weight']}, pulse {worst['case']['pulse']}",fontsize=11)
for ext in ('png','svg'):fig.savefig(F/f'worst-field.{ext}',dpi=160)
plt.close(fig)
h['status']='analysis_completed_pending_review';h['analysis']['completed_utc']=datetime.now(timezone.utc).isoformat();write('handoff.json',h)
print(json.dumps(dict(new_checks_passed=verification['new_checks_passed'],checks=len(checks),new_spatial_max=spatial_details['worst']['value'],reduction=summary['spatial_reduction_vs_previous_pair'],failed_names=[c['name'] for c in checks if not c['passed']],historical_failures_retained=len(historical)),indent=2))
