"""Verify retained raw fields; compare grids/timesteps and finite perturbations."""
from pathlib import Path
import hashlib, json
import numpy as np
from scipy import fft
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from solver import Solver, restrict

ROOT=Path(__file__).resolve().parent
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def field(run,t):return np.load(run['folder']/f'field-{t:.2f}.npz')['h']
def vals(r,k):return np.array([x[k] for x in r['rows']])
def key(r):
    s=r['settings'];return (s['family'],s['nu'],s['epsilon'],s['seed'],s['dt'],s['n'])
def prefix(rows,field='pass'):
    last=None
    for row in rows:
        if not row[field]:break
        last=row['t']
    return last
def relative(a,b,s):
    if a.shape==b.shape:return float(np.sqrt(s.inner(a-b,a-b)/s.inner(b,b)))
    n=s.n; nf=b.shape[1];c=restrict(b,n)*s.keep
    m=fft.fftfreq(nf)*nf;p=fft.rfftfreq(nf)*nf
    high=np.maximum.reduce(np.broadcast_arrays(abs(m[:,None,None]),abs(m[None,:,None]),p[None,None,:]))>=n/3
    weight=np.full((1,1,nf//2+1),2.);weight[...,0]=weight[...,-1]=1.
    norm=float(np.sum(abs(b)**2*weight)/nf**6)
    tail=float(np.sum(abs(b*high)**2*weight)/nf**6)
    return float(np.sqrt((s.inner(a-c,a-c)+tail)/norm))

def compare(a,b,kind):
    sa=a['settings'];s=Solver(sa['n'],sa['nu']);points=[]
    b_rows={round(r['t'],8):r for r in b['rows']};af={r['t']:r for r in a['fields']};bf={r['t']:r for r in b['fields']}
    for t in sorted(set(af)&set(bf)):
        aa=field(a,t);bb=field(b,t);rel=relative(aa,bb,s)
        ar=next(r for r in a['rows'] if abs(r['t']-t)<1e-8);br=b_rows[round(t,8)]
        werr=abs(ar['Wmax']-br['Wmax'])/br['Wmax']
        threshold=(.01,.05) if kind=='grid' else (.002,.01)
        checks=bool(rel<threshold[0] and werr<threshold[1])
        screen=bool(all(x['local_screen_pass'] for rr in [a,b] for x in rr['rows'] if x['t']<=t+1e-9))
        points.append(dict(t=t,velocity_relative_L2=rel,W_relative_difference=werr,comparison_pass=checks,local_pass=screen,pass_=checks and screen))
    q=dict(kind=kind,coarse=a['folder'].name,fine=b['folder'].name,points=points)
    q['last_contiguous_pass_time']=prefix(points,'pass_')
    return q

def main():
    runs=[];bad=[];verified=0;verrors=[]
    data=ROOT/'data-gpu' if (ROOT/'data-gpu').exists() else ROOT/'data'
    for p in sorted(data.glob('*/result.json')):
        r=json.loads(p.read_text());r['folder']=p.parent;runs.append(r)
        s=Solver(r['settings']['n'],r['settings']['nu'])
        for f in r['fields']:
            fp=p.parent/f['file'];h=np.load(fp)['h'];verified+=1
            row=next(x for x in r['rows'] if abs(x['t']-f['t'])<1e-8)
            E=.5*s.inner(h,h);Z=.5*s.inner(s.curl(h),s.curl(h))
            if digest(fp)!=f['sha256'] or not np.isfinite(h).all() or abs(E-row['energy'])>1e-12 or abs(Z-row['enstrophy'])>1e-10:verrors.append(str(fp))
        if digest(p.parent/'diagnostics.npz')!=r['diagnostics_sha256']:verrors.append(str(p.parent/'diagnostics.npz'))
        for fn,sha in r['source_sha256'].items():
            if digest(ROOT/fn)!=sha:verrors.append('Source changed: '+fn)
        if r['status']!='complete':bad.append(dict(run=p.parent.name,status=r['status'],error=r['error']))
    idx={key(r):r for r in runs}; comparisons=[]; perturbations=[];summary=[]
    for r in runs:
        f,nu,eps,seed,dt,n=key(r)
        if dt==.01:
            bigger=[rr for rr in runs if key(rr)[:5]==key(r)[:5] and rr['settings']['n']>n]
            if bigger:comparisons.append(compare(r,min(bigger,key=lambda x:x['settings']['n']),'grid'))
            smaller=[rr for rr in runs if key(rr)[:4]==key(r)[:4] and rr['settings']['n']==n and rr['settings']['dt']<dt]
            for half in smaller:comparisons.append(compare(r,half,'time'))
        if eps and dt==.01:
            base=idx.get((f,nu,0,1,dt,n))
            if base:
                s=Solver(n,nu);initial=field(r,0)-field(base,0);norm0=s.inner(initial,initial);ps=[]
                for ff in r['fields']:
                    t=ff['t'];delta=field(r,t)-field(base,t);u=s.real(delta);power=np.sum(u*u,axis=0)
                    ps.append(dict(t=t,amplification=float(np.sqrt(s.inner(delta,delta)/norm0)),
                        relative_to_baseline=float(np.sqrt(s.inner(delta,delta)/s.inner(field(base,t),field(base,t)))),
                        participation_fraction=float(np.mean(power)**2/np.mean(power**2)),
                        baseline_and_perturbed_local_pass=bool(next(x for x in r['rows'] if x['t']==t)['local_screen_pass'] and next(x for x in base['rows'] if x['t']==t)['local_screen_pass'])))
                perturbations.append(dict(run=r['folder'].name,n=n,epsilon=eps,seed=seed,points=ps))
        summary.append(dict(run=r['folder'].name,settings=r['settings'],status=r['status'],seconds=r['seconds'],steps=r['steps'],last_t=r['rows'][-1]['t'],
                            peak_W=float(vals(r,'Wmax').max()),peak_W_time=float(vals(r,'t')[vals(r,'Wmax').argmax()]),initial_W=r['rows'][0]['Wmax'],
                            W_amplification=float(vals(r,'Wmax').max()/r['rows'][0]['Wmax']),energy_remaining_fraction=r['rows'][-1]['energy']/r['rows'][0]['energy'],
                            max_tail_Z=float(vals(r,'tail_Z').max()),max_energy_residual=float(abs(vals(r,'energy_residual')).max()),
                            max_divergence=float(vals(r,'divergence_rms').max()),last_local_prefix=prefix(r['rows'],'local_screen_pass')))
    response_comparisons=[]
    for p in perturbations:
        candidates=[q for q in perturbations if q['epsilon']==p['epsilon'] and q['seed']==p['seed'] and q['n']>p['n']]
        if not candidates:continue
        q=min(candidates,key=lambda q:q['n']);n=p['n'];nf=q['n'];eps=p['epsilon'];seed=p['seed']
        a=idx[('kida',.005,eps,seed,.01,n)];b=idx[('kida',.005,eps,seed,.01,nf)]
        ba=idx[('kida',.005,0,1,.01,n)];bb=idx[('kida',.005,0,1,.01,nf)];s=Solver(n,.005);points=[]
        for j,pp in enumerate(p['points']):
            t=pp['t'];da=field(a,t)-field(ba,t);db=field(b,t)-field(bb,t)
            rel=relative(da,db,s);gain_error=abs(pp['amplification']-q['points'][j]['amplification'])/q['points'][j]['amplification']
            local=all(x['local_screen_pass'] for r in [a,b,ba,bb] for x in r['rows'] if x['t']<=t+1e-9)
            points.append(dict(t=t,response_field_relative_L2=rel,gain_relative_difference=gain_error,pass_=bool(rel<.01 and gain_error<.01 and local)))
        response_comparisons.append(dict(coarse=p['run'],fine=q['run'],points=points,last_contiguous_pass_time=prefix(points,'pass_')))
    planned=json.loads((ROOT/'planned_tasks.json').read_text())
    followup=json.loads((ROOT/'refinement-plan.json').read_text())['cases'] if (ROOT/'refinement-plan.json').exists() else []
    backend=[];cpu_verified=0
    if data.name=='data-gpu':
        for r in runs:
            cp_path=ROOT/'data'/r['folder'].name/'result.json'
            if not cp_path.exists():continue
            cpu=json.loads(cp_path.read_text());cpu['folder']=cp_path.parent
            if cpu['status']!='complete':continue
            cmp=compare(cpu,r,'backend')
            for f in cpu['fields']:
                fp=cpu['folder']/f['file'];h=np.load(fp)['h'];cpu_verified+=1
                if digest(fp)!=f['sha256'] or not np.isfinite(h).all():verrors.append(str(fp))
            mx=max(p['velocity_relative_L2'] for p in cmp['points'])
            backend.append(dict(run=r['folder'].name,max_checkpoint_relative_L2=mx,passed=mx<1e-9))
            if mx>=1e-9:verrors.append('CPU/GPU production comparison: '+r['folder'].name)
    dump(ROOT/'backend-comparison.json',dict(comparisons=backend,cpu_fields_verified=cpu_verified,passed=all(x['passed'] for x in backend)))
    dump(ROOT/'verification.json',dict(primary_data=data.name,verified_fields=verified,cpu_control_fields_verified=cpu_verified,errors=verrors,passed=not verrors,completed_runs=sum(r['status']=='complete' for r in runs),original_planned_runs=len(planned),followup_planned_runs=len(followup),planned_runs=len(planned)+len(followup),failed_runs=bad))
    seed_spread=[]
    seeded=[p for p in perturbations if p['n']==48 and p['epsilon']==.001]
    if len(seeded)==4:
        for j in range(len(seeded[0]['points'])):
            a=np.array([p['points'][j]['amplification'] for p in seeded]);rng=np.random.default_rng(20261009)
            boot=a[rng.integers(len(a),size=(10000,len(a)))].mean(axis=1)
            seed_spread.append(dict(t=seeded[0]['points'][j]['t'],seeds=4,mean=float(a.mean()),min=float(a.min()),max=float(a.max()),bootstrap_ci95=np.percentile(boot,[2.5,97.5]).tolist(),warning='Four random directions around one deterministic start; no population claim; late values can fail spatial resolution.'))
    dump(ROOT/'summary.json',dict(runs=summary,comparisons=comparisons,perturbations=perturbations,response_comparisons=response_comparisons,seed_spread=seed_spread,
        interpretation='Screen pass is descriptive convergence evidence at saved times, not a proof of continuum regularity. Any subsequent failure censors the reliable prefix. Broad initial perturbations do not measure causal propagation from a single particle. Response comparisons normalize by the perturbation difference itself; agreement of the much larger baseline flow cannot hide an inaccurate small response.'))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':145})
    colors={32:'#9ca3af',48:'#b45309',64:'#2563eb',96:'#059669',128:'#7c3aed',192:'#e11d48',256:'#111827'}
    fig,axes=plt.subplots(2,3,figsize=(13,7),sharex=True)
    for j,nu in enumerate([.02,.005,.001]):
        for r in runs:
            f,v,e,seed,dt,n=key(r)
            if f=='kida' and v==nu and e==0 and dt==.01:
                tt=vals(r,'t');ww=vals(r,'Wmax');bad_indices=np.flatnonzero(~vals(r,'local_screen_pass'))
                last=max(0,int(bad_indices[0])-1) if len(bad_indices) else len(tt)-1
                axes[0,j].plot(tt[:last+1],ww[:last+1],color=colors[n])
                if last<len(tt)-1:axes[0,j].plot(tt[last:],ww[last:],color=colors[n],ls='--')
                axes[1,j].semilogy(vals(r,'t'),np.maximum(vals(r,'tail_Z'),1e-12),color=colors[n])
        axes[0,j].set_title(f'Kida–Pelz, viscosity {nu:g}');axes[0,j].set_ylabel('Maximum vorticity')
        axes[1,j].axhline(.01,color='#dc2626',ls='--',label='1% screen');axes[1,j].set_ylim(1e-10,1);axes[1,j].set_ylabel('Enstrophy near cutoff');axes[1,j].set_xlabel('Model time')
    present=sorted({r['settings']['n'] for r in runs if r['settings']['family']=='kida'})
    handles=[plt.Line2D([0],[0],color=colors[n],label=f'{n}³') for n in present]
    handles.append(plt.Line2D([0],[0],color='#444444',ls='--',label='after local screen failure'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.955),ncol=4,fontsize=9)
    fig.suptitle('Strong vortex growth versus loss of spatial resolution',fontsize=15);fig.tight_layout(rect=(0,0,1,.88));fig.savefig(ROOT/'stress.png');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    for nu in [.02,.005,.001]:
        r=idx.get(('kida',nu,0,1,.01,96))
        if r:
            axes[0].plot(vals(r,'t'),vals(r,'energy')/r['rows'][0]['energy'],label=f'ν={nu:g}')
            axes[1].plot(vals(r,'t'),vals(r,'stretch_production')-vals(r,'enstrophy_dissipation'),label=f'ν={nu:g}')
            axes[2].semilogy(vals(r,'t'),np.maximum(abs(vals(r,'energy_residual')),1e-14),label=f'ν={nu:g}')
    axes[0].set_ylabel('Energy / initial energy');axes[1].set_ylabel('Stretching − viscous enstrophy loss');axes[2].set_ylabel('|Energy budget residual|');axes[2].axhline(.001,color='r',ls='--')
    for ax in axes:ax.set_xlabel('Model time');ax.legend(fontsize=8)
    fig.suptitle('96³ diagnostics: energy balance alone cannot establish resolution');fig.tight_layout();fig.savefig(ROOT/'budgets.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for c in comparisons:
        if c['kind']=='grid' and c['coarse'].startswith('kida-') and '-e0-' in c['coarse']:
            label=c['coarse'].replace('kida-','').split('-e')[0]+' → '+c['fine'].split('-')[1]
            axes[0].semilogy([p['t'] for p in c['points']],np.maximum([p['velocity_relative_L2'] for p in c['points']],1e-12),label=label)
    axes[0].axhline(.01,color='r',ls='--');axes[0].set_ylabel('Whole-field spatial relative error');axes[0].legend(fontsize=6,ncol=2)
    for c in comparisons:
        if c['kind']=='time':axes[1].semilogy([p['t'] for p in c['points']],np.maximum([p['velocity_relative_L2'] for p in c['points']],1e-12),label=c['coarse'].replace('kida-','').replace('-s1-dt0.01',''))
    axes[1].axhline(.002,color='r',ls='--');axes[1].set_ylabel('Whole-field timestep relative error');axes[1].legend(fontsize=8)
    for ax in axes:ax.set_xlabel('Model time')
    fig.suptitle('Independent checks of spatial and temporal accuracy');fig.tight_layout();fig.savefig(ROOT/'convergence.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for p in perturbations:
        label=f"{p['n']}³, ε={p['epsilon']:g}, seed {p['seed']}";x=[v['t'] for v in p['points']]
        axes[0].plot(x,[v['amplification'] for v in p['points']],label=label)
        axes[1].plot(x,[v['participation_fraction'] for v in p['points']],label=label)
    axes[0].set_ylabel('Velocity difference / starting difference');axes[1].set_ylabel('Difference participation / domain volume');axes[1].set_ylim(0,1)
    for ax in axes:ax.set_xlabel('Model time')
    axes[0].legend(fontsize=6,ncol=2);fig.suptitle('Finite-time sensitivity to energy-matched flow perturbations\nLater values may be spatially unresolved; these are not molecular tests',fontsize=12);fig.tight_layout();fig.savefig(ROOT/'perturbations.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for r in runs:
        f,nu,eps,seed,dt,n=key(r)
        if f=='taylor_green':
            axes[0].plot(vals(r,'t'),vals(r,'Wmax'),label=f'{n}³, ν={nu:g}')
            axes[1].semilogy(vals(r,'t'),np.maximum(vals(r,'tail_Z'),1e-12),label=f'{n}³, ν={nu:g}')
    axes[0].set_ylabel('Maximum vorticity');axes[1].set_ylabel('Enstrophy near cutoff');axes[1].axhline(.01,color='r',ls='--')
    for ax in axes:ax.set_xlabel('Model time')
    axes[0].legend(fontsize=8);fig.suptitle('Second smooth starting field: three-dimensional Taylor–Green');fig.tight_layout();fig.savefig(ROOT/'taylor-green.png');plt.close(fig)
    from plot_slices import main as plot_slices
    plot_slices()
    print(json.dumps(dict(runs=len(runs),fields=verified,errors=verrors,comparisons=len(comparisons)),indent=2))

if __name__=='__main__':main()
